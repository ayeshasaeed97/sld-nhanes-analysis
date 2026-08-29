# ckm baseline
library(tidyverse)
library(flextable)
library(survey)
library(gtsummary)

ckm <- read_csv('clean_data/ckm.csv')


ckm_table1 <- ckm %>%
  mutate(
    
    # ----------------------------------------------------
    # Combined MASLD / CKM classification
    # ----------------------------------------------------
    
    MASLD_CKM_GROUP = factor(
      case_when(
        LIVER_GROUP4 == "Neither" & CKM %in% c("0", "1", "2") ~ "No MASLD: CKM 0-2",
        LIVER_GROUP4 == "Neither" & CKM %in% c("3", "4") ~ "No MASLD: CKM 3-4",
        LIVER_GROUP4 == "MASLD" & CKM %in% c("0", "1", "2") ~ "MASLD: CKM 0-2",
        LIVER_GROUP4 == "MASLD" & CKM %in% c("3", "4") ~ "MASLD: CKM 3-4",
        TRUE ~ NA_character_
      ),
      levels = c(
        "No MASLD: CKM 0-2",
        "No MASLD: CKM 3-4",
        "MASLD: CKM 0-2",
        "MASLD: CKM 3-4"
      )
    ),
    
    
    # ----------------------------------------------------
    # Baseline characteristics
    # ----------------------------------------------------
    
    # Sex
    sex = factor(
      HSSEX,
      levels = c(1, 2),
      labels = c("Men", "Women")
    ),
    
    
    # Race/ethnicity
    race_ethnicity = factor(
      case_when(
        DMARETHN == 1 ~ "Non-Hispanic white",
        DMARETHN == 2 ~ "Non-Hispanic black",
        DMARETHN == 3 ~ "Mexican-American",
        DMARETHN == 4 ~ "Other",
        TRUE ~ NA_character_
      ),
      levels = c(
        "Non-Hispanic white",
        "Non-Hispanic black",
        "Mexican-American",
        "Other"
      )
    ),
    
    
    # Education
    education_lt12 = factor(
      case_when(
        HFA8R < 12 ~ "Yes",
        HFA8R >= 12 ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    
    # BMI categories
    bmi_cat = factor(
      case_when(
        BMPBMI < 18.5 ~ "<18.5",
        BMPBMI >= 18.5 & BMPBMI < 25 ~ "18.5-24.9",
        BMPBMI >= 25 & BMPBMI < 30 ~ "25-29.9",
        BMPBMI >= 30 & BMPBMI < 35 ~ "30-34.9",
        BMPBMI >= 35 ~ "≥35",
        TRUE ~ NA_character_
      ),
      levels = c(
        "<18.5",
        "18.5-24.9",
        "25-29.9",
        "30-34.9",
        "≥35"
      )
    ),
    
    
    # High waist circumference
    high_waist = factor(
      case_when(
        HSSEX == 1 & BMPWAIST > 102 ~ "Yes",
        HSSEX == 2 & BMPWAIST > 88 ~ "Yes",
        !is.na(BMPWAIST) ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    
    # Diabetes
    # prevent_ckm already creates diabetes as TRUE/FALSE
    diabetes_table = factor(
      case_when(
        diabetes ~ "Yes",
        !diabetes ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    
    # Hypertension
    # prevent_ckm already creates hypertension as TRUE/FALSE
    hypertension_table = factor(
      case_when(
        hypertension ~ "Yes",
        !hypertension ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    
    # Sedentary
    sedentary = factor(
      case_when(
        HAT2 == 2 & HAT4 == 2 & HAT6 == 2 & HAT8 == 2 &
          HAT10 == 2 & HAT12 == 2 & HAT14 == 2 &
          HAT16 == 2 & HAT18 == 2 ~ "Yes",
        
        if_any(
          c(
            HAT2, HAT4, HAT6, HAT8,
            HAT10, HAT12, HAT14, HAT16, HAT18
          ),
          ~ .x == 1
        ) ~ "No",
        
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    
    # History of cardiovascular disease
    history_cvd = factor(
      case_when(
        HAF10 == 1 | HAC1D == 1 | HAC1C == 1 ~ "Yes",
        HAF10 == 2 & HAC1D == 2 & HAC1C == 2 ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    
    # History of cancer
    history_cancer = factor(
      case_when(
        HAC1N == 1 | HAC1O == 1 ~ "Yes",
        HAC1N == 2 & HAC1O == 2 ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    
    # Smoking
    smoking_status = factor(
      case_when(
        HAR1 == 2 ~ "Never",
        HAR1 == 1 & HAR3 == 1 ~ "Current",
        HAR1 == 1 & HAR3 == 2 ~ "Former",
        TRUE ~ NA_character_
      ),
      levels = c("Never", "Former", "Current")
    ),
    
    current_smoking = factor(
      case_when(
        smoking_status == "Current" ~ "Yes",
        smoking_status %in% c("Never", "Former") ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    
    # Alcohol
    alcohol_category = factor(
      case_when(
        MAPE1 == 2 ~ "Never",
        ALC == "Low-Moderate alcohol consumption" ~ "Low-moderate",
        TRUE ~ NA_character_
      ),
      levels = c("Never", "Low-moderate")
    )
  ) %>%
  
  # Remove people who do not belong to one of our 5 table columns
  filter(!is.na(MASLD_CKM_GROUP))


table(ckm_table1$MASLD_CKM_GROUP, useNA = "ifany")

options(survey.lonely.psu = "adjust")

nhanes_design_ckm <- svydesign(
  ids = ~SDPPSU6,
  strata = ~SDPSTRA6,
  weights = ~WTPFEX6,
  nest = TRUE,
  data = ckm_table1
)

baseline_ckm <- nhanes_design_ckm %>%
  tbl_svysummary(
    by = MASLD_CKM_GROUP,
    
    include = c(
      HSAGEIR,
      sex,
      race_ethnicity,
      education_lt12,
      bmi_cat,
      high_waist,
      diabetes_table,
      hypertension_table,
      hypercholesterolaemia,
      sedentary,
      history_cvd,
      history_cancer,
      current_smoking,
      alcohol_category,
      GHP,
      ATPSI,
      ASPSI,
      GGPSI
    ),
    
    value = list(
      education_lt12 ~ "Yes"
    ),
    
    label = list(
      HSAGEIR ~ "Mean age (years)",
      sex ~ "Sex",
      race_ethnicity ~ "Race or ethnicity",
      education_lt12 ~ "<12 years education",
      bmi_cat ~ "Body mass index",
      high_waist ~ "High waist circumference",
      diabetes_table ~ "Diabetes",
      hypertension_table ~ "Hypertension",
      hypercholesterolaemia ~ "Hypercholesterolaemia",
      sedentary ~ "Sedentary",
      history_cvd ~ "History of cardiovascular disease",
      history_cancer ~ "History of cancer",
      current_smoking ~ "Current smoking",
      alcohol_category ~ "Alcohol consumption",
      GHP ~ "Mean (SE) glycated haemoglobin (%)",
      ATPSI ~ "Alanine aminotransferase (U/L)",
      ASPSI ~ "Aspartate aminotransferase (U/L)",
      GGPSI ~ "γ-glutamyltransferase (U/L)"
    ),
    
    statistic = list(
      all_continuous() ~ "{mean} ({mean.std.error})",
      all_categorical() ~ "{p}%"
    ),
    
    digits = list(
      all_continuous() ~ c(1, 1),
      all_categorical() ~ 1
    ),
    
    missing = "no"
  ) %>%
  
  add_p() %>%
  
  modify_header(
    label ~ "**Characteristics**",
    stat_1 ~ "**No MASLD: CKM 0-2**<br>N = {n}",
    stat_2 ~ "**No MASLD: CKM 3-4**<br>N = {n}",
    stat_3 ~ "**MASLD: CKM 0-2**<br>N = {n}",
    stat_4 ~ "**MMASLD: CKM 3-4**<br>N = {n}",
    p.value ~ "**p-value**"
  ) %>%
  
  bold_labels()



baseline_ckm #%>%
  #as_flex_table() %>%
  #save_as_docx(
   # path = "baseline_characteristics_masld_ckm.docx"
#  )


#write_csv(
  ckm_table1,
  "clean_data/ckm_table1.csv"
)


table(ckm_table1$MASLD_CKM_GROUP)
