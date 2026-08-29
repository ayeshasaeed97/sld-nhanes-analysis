# Dataset creation for CKM analysis----
############################################

setwd("~/Documents/FL_NHANES")

# loading packages ----
library(tidyverse)
library(flextable)
library(haven) # for reading xpt files
library(readr)

# loading masld data with prevent score data
prevent <- read_csv("clean_data/masld_prevent.csv")
prevent <- prevent %>%
  filter(over_years == 10)

# need to add additional variables due to ckm classification:
# 

# loading data -------------------------------
adult <- read.csv('raw_data/adult_output.csv')
lab <- read.csv('raw_data/lab_output.csv')
exam <- read.csv('raw_data/exam_output.csv')

# Variable Selection ----------------------------------

# lab
lab_vars <- c(
  "SEQN", "TGP", "HDP", "UBP", "URP", "G1P", "WTPFHSD6", "G2P", "G1PCODE"
)

# G1P = first plasma glucose
# G2P = second plasma glucose
# WTPFHSD6 = weights for fasting subsample
# G1PCODE = incomplete glucose test code

lab_value_vars <- c("TGP", "HDP", "UBP", "URP", "G1P", "G2P")

lab_clean <- lab %>%
  select(all_of(lab_vars)) %>%
  mutate(
    across(
      all_of(lab_value_vars),
      ~ na_if(.x, 8888)
    ),
    across(
      all_of(lab_value_vars),
      ~ na_if(.x, 888)
    )
  )

adult_clean <- adult %>%
  select(SEQN, HAD10, HAD6)


# merging data by SEQN number
prevent <- prevent %>%
  left_join(adult_clean, by = "SEQN") %>%
  left_join(lab_clean,  by = "SEQN")


# ckm------
# ckm classification

prevent <- prevent %>%
  mutate(
  LIVER_GROUP4 = factor(
  case_when(
    GUPHSPFR == 2 & LIVER_ENZ == "Normal" & MASLD_CRITERIA_COUNT >= 1 ~ "MASLD",
    
    GUPHSPFR == 2 &
      LIVER_ENZ == "Raised" &
      HBP == 2 &
      HCP == 2 &
      IRON_OVERLOAD == "No iron overload" ~ "MASLD",
    
    GUPHSPFR == 1 ~ "Neither",
    
    TRUE ~ NA_character_
  ),
  levels = c("MASLD", "Neither")
  )
  )

table(prevent$LIVER_GROUP2, prevent$MASLD_CRITERIA_COUNT)
table(prevent$LIVER_GROUP4, prevent$MASLD_CRITERIA_COUNT)

# ckm

prevent_ckm <- prevent %>% 
  mutate(
    male = HSSEX == 1,
    female = HSSEX == 2,
    
    # identifying those who were AM participants with fasting
    morning_fasting = !is.na(WTPFHSD6) & WTPFHSD6 > 0,

    # Stage 0 / Stage 1 adiposity criteria
    
    normal_bmi = BMPBMI < 25,
    elevated_bmi = BMPBMI >= 25,
    
    # waist circumference
    normal_wc = case_when(
      male ~ BMPWAIST < 102,
      female ~ BMPWAIST < 88,
      TRUE ~ NA
    ),
    
    elevated_wc = case_when(
      male ~ BMPWAIST >= 102,
      female ~ BMPWAIST >= 88,
      TRUE ~ NA
    ),
    
    
    # Glucose/OGTT validity indicators
    complete_pgtt = is.na(G1PCODE),
    
    valid_fasting_glucose =
      morning_fasting & !is.na(G1P),
    
    valid_2hr_glucose =
      morning_fasting & complete_pgtt & !is.na(G2P),
    
    # Diabetes
    diabetes =
      coalesce(GHP >= 6.5, FALSE) |
      coalesce(HAD10 == 1, FALSE) |
      coalesce(HAD6 == 1, FALSE) |
      coalesce(valid_fasting_glucose & G1P >= 126, FALSE) |
      coalesce(valid_2hr_glucose & G2P >= 200, FALSE),
    # add groups of people who don't fall into fasting glucose
    
    # Prediabetes
    prediabetes =
      !diabetes &
      (
        coalesce(GHP >= 5.7 & GHP < 6.5, FALSE) |
          coalesce(valid_fasting_glucose & G1P >= 100 & G1P < 126, FALSE) |
          coalesce(valid_2hr_glucose & G2P >= 140 & G2P < 200, FALSE)
      ),
         
         
    # Blood pressure / hypertension
    elevated_bp =
      PEPMNK1R >= 130 | #SBP
      PEPMNK5R >= 80 |
      HAE5A == 1,
    
    hypertension =
      elevated_bp |
      HAE2 == 1, # hypertension diagnosis

    # Triglycerides and HDL
    
    high_tg_135 = TGP >= 135,
    high_tg_150 = TGP >= 150,
    
    low_hdl = case_when(
      male ~ HDP < 40,
      female ~ HDP < 50,
      TRUE ~ NA
    ),
    

    # Metabolic syndrome
    # >=3 of:
    # elevated waist, low HDL, TG >=150, elevated BP, prediabetes
    
    metabolic_syndrome_count = rowSums(
      cbind(
        elevated_wc,
        low_hdl,
        high_tg_150,
        elevated_bp,
        prediabetes
      ),
      na.rm = TRUE
    ),
    
    metabolic_syndrome = metabolic_syndrome_count >= 3,
    
    UACR = (UBP/URP)*100,
    
    # CKD KDIGO risk category
    # Replace egfr and UACR with your variables
    # egfr units: mL/min/1.73m2
    # UACR units: mg/g

    
    kdigo_risk = case_when(
      is.na(egfr) | is.na(UACR) ~ NA_character_,
      
      egfr >= 90 & UACR < 30 ~ "low",
      egfr >= 90 & UACR >= 30 & UACR < 300 ~ "moderate",
      egfr >= 90 & UACR >= 300 ~ "high",
      
      egfr >= 60 & egfr < 90 & UACR < 30 ~ "low",
      egfr >= 60 & egfr < 90 & UACR >= 30 & UACR < 300 ~ "moderate",
      egfr >= 60 & egfr < 90 & UACR >= 300 ~ "high",
      
      egfr >= 45 & egfr < 60 & UACR < 30 ~ "moderate",
      egfr >= 45 & egfr < 60 & UACR >= 30 & UACR < 300 ~ "high",
      egfr >= 45 & egfr < 60 & UACR >= 300 ~ "very_high",
      
      egfr >= 30 & egfr < 45 & UACR < 30 ~ "high",
      egfr >= 30 & egfr < 45 & UACR >= 30 ~ "very_high",
      
      egfr < 30 ~ "very_high",
      
      TRUE ~ NA_character_
    ),
    
    moderate_high_ckd = kdigo_risk %in% c("moderate", "high"),
    very_high_ckd = kdigo_risk == "very_high",
    

    # PREVENT 10-year CVD risk
    # Replace PREVENT_10YR_CVD with your variable
    # If risk is stored as 22.5, converts to 0.225
    # If risk is already 0.225, leaves as 0.225

    
    prevent10_risk = if_else(
      total_cvd > 1,
      total_cvd / 100,
      total_cvd
    ),
    
    high_prevent_risk = prevent10_risk >= 0.20,
    

    # Established cardiovascular disease
    # Replace these with your CVD variables
    # CKM stage 4 includes:
    # coronary heart disease, angina, heart attack, heart failure, stroke
    # Assuming 1 = yes, 2 = no

    
    established_cvd =
      HAC1C == 1 |
      HAF10 == 1 |
      HAC1D == 1,
    

    # CKM stage indicators
    
    ckm_stage0 = normal_bmi & normal_wc,
    
    ckm_stage1 =
      elevated_bmi |
      elevated_wc |
      prediabetes,
    
    ckm_stage2 =
      high_tg_135 |
      hypertension |
      diabetes |
      metabolic_syndrome |
      moderate_high_ckd,
    
    ckm_stage3 =
      very_high_ckd |
      high_prevent_risk,
    
    ckm_stage4 =
      established_cvd,
    

    # Final CKM classification
    
    CKM = case_when(
      ckm_stage4 ~ "4",
      ckm_stage3 ~ "3",
      ckm_stage2 ~ "2",
      ckm_stage1 ~ "1",
      ckm_stage0 ~ "0",
      TRUE ~ NA_character_
    ),
    
    CKM = factor(CKM, levels = c("0", "1", "2", "3", "4"))
  )


table(prevent_ckm$LIVER_GROUP4, prevent_ckm$CKM, useNA = "ifany")


write_csv(prevent_ckm, "clean_data/ckm.csv")

