# loading packages
library(tidyverse)
library(flextable)
library(survey)
library(gtsummary)

# loading working data

nhanes_work <- read.csv("clean_data/nhanes_work.csv")

nhanes_table1 <- nhanes_work %>%
  mutate(
    # Liver group labels
    LIVER_GROUP = factor(
      LIVER_GROUP,
      levels = c("Neither", "NAFLD", "NASH"),
      labels = c(
        "No hepatic steatosis",
        "Non-alcoholic fatty liver disease",
        "Non-alcoholic steatohepatitis"
      )
    ),
    
    # Sex
    sex = factor(
      HSSEX,
      levels = c(1, 2),
      labels = c("Men", "Women")
    ),
    
    men = factor(
      if_else(HSSEX == 1, "Yes", "No"),
      levels = c("No", "Yes")
    ),
    
    # Race/ethnicity
    # Prefer DMARETHN if available because it is already race-ethnicity.
    race_ethnicity = factor(
      case_when(
        DMARETHN == 1 ~ "Non-Hispanic white",
        DMARETHN ==2 ~ "Non-Hispanic black",
        DMARETHN == 3 ~ "Mexican-American",
        DMARETHN == 4 ~ "Other"
      ),
      levels = c("Non-Hispanic white", "Non-Hispanic black", "Mexican-American", "Other")
    ),
    
    # Education: <12 years
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
      levels = c("<18.5", "18.5-24.9", "25-29.9", "30-34.9", "≥35")
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
    
    # Diabetes: self-reported only here unless you add glucose/medication data
    # Paper also used diabetes meds and glucose.
    diabetes = diabetes,
    
    # Hypertension: self-report OR measured BP high
    # Paper also used medication.
    hypertension = factor(
      case_when(
        HAE2 == 1 | PEPMNK1R >= 140 | PEPMNK5R >= 90 ~ "Yes",
        HAE2 == 2 & PEPMNK1R < 140 & PEPMNK5R < 90 ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    # Hypercholesterolaemia: self-report only here unless total cholesterol/meds added
    hypercholesterolaemia = hypercholesterolaemia,
    
    # Sedentary: no to all listed activities
    sedentary = factor(
      case_when(
        HAT2 == 2 & HAT4 == 2 & HAT6 == 2 & HAT8 == 2 &
          HAT10 == 2 & HAT12 == 2 & HAT14 == 2 &
          HAT16 == 2 & HAT18 == 2 ~ "Yes",
        if_any(
          c(HAT2, HAT4, HAT6, HAT8, HAT10, HAT12, HAT14, HAT16, HAT18),
          ~ .x == 1
        ) ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    # History of cardiovascular disease:
    # heart attack, stroke, or congestive heart failure
    history_cvd = factor(
      case_when(
        HAF10 == 1 | HAC1D == 1 | HAC1C == 1 ~ "Yes",
        HAF10 == 2 & HAC1D == 2 & HAC1C == 2 ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    # History of cancer:
    # skin cancer or other cancer
    history_cancer = factor(
      case_when(
        HAC1N == 1 | HAC1O == 1 ~ "Yes",
        HAC1N == 2 & HAC1O == 2 ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    # Smoking: never/former/current
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
  )

nhanes_table1 <- nhanes_table1 %>%
  filter(!is.na(LIVER_GROUP))

options(survey.lonely.psu = "adjust")

nhanes_design_t1 <- svydesign(
  ids = ~SDPPSU6,
  strata = ~SDPSTRA6,
  weights = ~WTPFEX6,
  nest = TRUE,
  data = nhanes_table1
)


baseline_table <- nhanes_design_t1 %>%
  tbl_svysummary(
    by = LIVER_GROUP,
    include = c(
      HSAGEIR,
      sex,
      race_ethnicity,
      education_lt12,
      bmi_cat,
      high_waist,
      diabetes,
      hypertension,
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
      diabetes ~ "Diabetes",
      hypertension ~ "Hypertension",
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
      all_categorical() ~ c(1)
    ),
    missing = "no"
  ) %>%
  add_p() %>%
  modify_header(label ~ "**Characteristics**") %>%
  bold_labels()


baseline_table <- baseline_table %>%
  modify_header(
    stat_1 ~ "**No hepatic steatosis**<br>N = 9,169",
    stat_2 ~ "**Non-alcoholic fatty liver disease**<br>N = 2,179",
    stat_3 ~ "**Non-alcoholic steatohepatitis**<br>N = 416",
    p.value ~ "**p-value**"
  )

baseline_table

# exporting baseline_table

#baseline_table %>%
#   as_flex_table() %>%
#   save_as_docx(path = "/Users/ayeshasaeed/Documents/MSAB/Practical Training/NHANES/baseline_characteristics_nhanes.docx")


# exporting table1 data

write_csv(nhanes_table1, "clean_data/nhanes_table1.csv")
