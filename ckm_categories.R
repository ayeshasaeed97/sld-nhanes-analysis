# ckm categories
setwd("~/Documents/FL_NHANES")

library(tidyverse)

# loading masld data with prevent score data
prevent <- read_csv("clean_data/masld_prevent.csv")

# 1582, ref group, no HS diagnosis, no CMRF
# 0-1 is ref group, 2, 3, 4, 5
# additional analysis 0-2 reference, and 3+
# ckm classification


# ckm stage only use 10-year risk, so will filter prevent score dataset to only have observations with 10-year risk

prevent <- prevent %>%
  filter(over_years == 10) %>%
  mutate(
    LIVER_GROUP3 = factor(
      case_when(
        LIVER_GROUP2 == "MASLD" | 
          LIVER_GROUP2 == "MASH" ~ "MASLD",
        LIVER_GROUP2 == "Neither" ~ "None",
        TRUE ~ NA_character_
      ),
      levels = c("None", "MASLD")
    ), 
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
  filter(over_years == 10) %>%
  mutate(
    
    # -------------------------
    # Sex
    # Assuming HSSEX: 1 = male, 2 = female
    # -------------------------
    
    male = HSSEX == 1,
    female = HSSEX == 2,
    
    # -------------------------
    # Stage 0 / Stage 1 adiposity criteria
    # No Asian race adjustment
    # -------------------------
    
    normal_bmi = BMPBMI < 25,
    elevated_bmi = BMPBMI >= 25,
    
    normal_waist = case_when(
      male ~ BMPWAIST < 102,
      female ~ BMPWAIST < 88,
      TRUE ~ NA
    ),
    
    elevated_waist = case_when(
      male ~ BMPWAIST >= 102,
      female ~ BMPWAIST >= 88,
      TRUE ~ NA
    ),
    
    # -------------------------
    # Prediabetes
    # Replace A1C and FASTING_GLUCOSE with your variables
    # -------------------------
    
    prediabetes = 
      (A1C >= 5.7 & A1C < 6.5) |
      (FASTING_GLUCOSE >= 100 & FASTING_GLUCOSE < 126),
    
    # -------------------------
    # Diabetes
    # Replace diabetes diagnosis/medication variables as needed
    # Assuming 1 = yes, 2 = no for diagnosis variables
    # -------------------------
    
    diabetes =
      A1C >= 6.5 |
      FASTING_GLUCOSE >= 126 |
      DIABETES_DX == 1 |
      DIABETES_MED == 1,
    
    # -------------------------
    # Blood pressure / hypertension
    # Replace SBP, DBP, BP_MED, HTN_DX with your variables
    # -------------------------
    
    elevated_bp =
      SBP >= 130 |
      DBP >= 80 |
      BP_MED == 1,
    
    hypertension =
      elevated_bp |
      HTN_DX == 1,
    
    # -------------------------
    # Triglycerides and HDL
    # Replace TRIGLYCERIDES and HDL with your variables
    # -------------------------
    
    high_tg_135 = TRIGLYCERIDES >= 135,
    high_tg_150 = TRIGLYCERIDES >= 150,
    
    low_hdl = case_when(
      male ~ HDL < 40,
      female ~ HDL < 50,
      TRUE ~ NA
    ),
    
    # -------------------------
    # Metabolic syndrome
    # >=3 of:
    # elevated waist, low HDL, TG >=150, elevated BP, prediabetes
    # -------------------------
    
    metabolic_syndrome_count = rowSums(
      cbind(
        elevated_waist,
        low_hdl,
        high_tg_150,
        elevated_bp,
        prediabetes
      ),
      na.rm = TRUE
    ),
    
    metabolic_syndrome = metabolic_syndrome_count >= 3,
    
    # -------------------------
    # CKD KDIGO risk category
    # Replace EGFR and UACR with your variables
    # EGFR units: mL/min/1.73m2
    # UACR units: mg/g
    # -------------------------
    
    kdigo_risk = case_when(
      is.na(EGFR) | is.na(UACR) ~ NA_character_,
      
      EGFR >= 90 & UACR < 30 ~ "low",
      EGFR >= 90 & UACR >= 30 & UACR < 300 ~ "moderate",
      EGFR >= 90 & UACR >= 300 ~ "high",
      
      EGFR >= 60 & EGFR < 90 & UACR < 30 ~ "low",
      EGFR >= 60 & EGFR < 90 & UACR >= 30 & UACR < 300 ~ "moderate",
      EGFR >= 60 & EGFR < 90 & UACR >= 300 ~ "high",
      
      EGFR >= 45 & EGFR < 60 & UACR < 30 ~ "moderate",
      EGFR >= 45 & EGFR < 60 & UACR >= 30 & UACR < 300 ~ "high",
      EGFR >= 45 & EGFR < 60 & UACR >= 300 ~ "very_high",
      
      EGFR >= 30 & EGFR < 45 & UACR < 30 ~ "high",
      EGFR >= 30 & EGFR < 45 & UACR >= 30 ~ "very_high",
      
      EGFR < 30 ~ "very_high",
      
      TRUE ~ NA_character_
    ),
    
    moderate_high_ckd = kdigo_risk %in% c("moderate", "high"),
    very_high_ckd = kdigo_risk == "very_high",
    
    # -------------------------
    # PREVENT 10-year CVD risk
    # Replace PREVENT_10YR_CVD with your variable
    # If risk is stored as 22.5, converts to 0.225
    # If risk is already 0.225, leaves as 0.225
    # -------------------------
    
    prevent10_risk = if_else(
      PREVENT_10YR_CVD > 1,
      PREVENT_10YR_CVD / 100,
      PREVENT_10YR_CVD
    ),
    
    high_prevent_risk = prevent10_risk >= 0.20,
    
    # -------------------------
    # Established cardiovascular disease
    # Replace these with your CVD variables
    # CKM stage 4 includes:
    # coronary heart disease, angina, heart attack, heart failure, stroke
    # Assuming 1 = yes, 2 = no
    # -------------------------
    
    established_cvd =
      HEART_FAILURE == 1 |
      CORONARY_HEART_DISEASE == 1 |
      ANGINA == 1 |
      HEART_ATTACK == 1 |
      STROKE == 1,
    
    # -------------------------
    # CKM stage indicators
    # -------------------------
    
    ckm_stage0 = normal_bmi & normal_waist,
    
    ckm_stage1 =
      elevated_bmi |
      elevated_waist |
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
    
    # -------------------------
    # Final CKM classification
    # Highest stage wins
    # -------------------------
    
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

