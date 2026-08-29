# survival analysis with ckm groups

############################################################
# Survey-weighted Cox regression by joint MASLD + CKM stage
#
# Exposure:
#   No MASLD          <- reference
#   MASLD: CKM 1
#   MASLD: CKM 2
#   MASLD: CKM 3
#   MASLD: CKM 4
#
# Model 1:
#   CKM exposure + sex + race/ethnicity
#
# Model 2:
#   Model 1 + education + smoking + alcohol + sedentary
#
# IMPORTANT:
# CKM stage itself contains cardiometabolic information.
# Therefore do NOT additionally adjust for BMI, waist circumference,
# hypertension, diabetes, triglycerides, HDL, CKD, CVD, etc.
############################################################

setwd("~/Documents/FL_NHANES")

library(tidyverse)
library(survey)
library(survival)


# ---------------------------------------------------------
# Load CKM baseline dataset
# ---------------------------------------------------------

nhanes_table1 <- read_csv(
  "clean_data/ckm_table1.csv",
  show_col_types = FALSE
)

EXPOSURE <- "MASLD_CKM_GROUP"


# ---------------------------------------------------------
# Check exposure exists
# ---------------------------------------------------------

if (!EXPOSURE %in% names(nhanes_table1)) {
  stop(
    paste0(
      EXPOSURE,
      " is missing. Re-run the CKM classification and baseline table code."
    )
  )
}


# ---------------------------------------------------------
# Set exposure factor levels
#
# First level = Cox model reference category
# ---------------------------------------------------------

nhanes_table1$masld

nhanes_table1 <- nhanes_table1 %>%
  mutate(
    MASLD_CKM_GROUP = factor(
      MASLD_CKM_GROUP,
      levels = c(
        "No MASLD",
        "MASLD: CKM 1",
        "MASLD: CKM 2",
        "MASLD: CKM 3",
        "MASLD: CKM 4"
      )
    )
  )


# Remove people without a valid exposure category
nhanes_table1 <- nhanes_table1 %>%
  filter(!is.na(.data[[EXPOSURE]]))


# Check sample sizes
table(nhanes_table1$MASLD_CKM_GROUP, useNA = "ifany")


# ---------------------------------------------------------
# Mortality outcomes
# ---------------------------------------------------------

nhanes_surv <- nhanes_table1 %>%
  mutate(
    
    # Age as the survival time scale
    age_entry = HSAGEIR,
    age_exit = HSAGEIR + permth_exm / 12,
    
    
    # -----------------------------------------------------
    # All-cause mortality
    # -----------------------------------------------------
    
    all_cause_death = if_else(
      mortstat == 1,
      1,
      0
    ),
    
    
    # -----------------------------------------------------
    # Cardiometabolic mortality
    #
    # 1 = heart disease
    # 5 = cerebrovascular disease
    # 7 = diabetes
    # -----------------------------------------------------
    
    cardiometabolic_death = case_when(
      mortstat == 1 & ucod_leading %in% c(1, 5, 7) ~ 1,
      mortstat == 1 & !(ucod_leading %in% c(1, 5, 7)) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    
    # -----------------------------------------------------
    # Respiratory mortality
    #
    # 3 = chronic lower respiratory disease
    # 8 = influenza / pneumonia
    # -----------------------------------------------------
    
    respiratory_death = case_when(
      mortstat == 1 & ucod_leading %in% c(3, 8) ~ 1,
      mortstat == 1 & !(ucod_leading %in% c(3, 8)) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    
    # -----------------------------------------------------
    # Cancer mortality
    # -----------------------------------------------------
    
    cancer_death = case_when(
      mortstat == 1 & ucod_leading == 2 ~ 1,
      mortstat == 1 & ucod_leading != 2 ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    
    # -----------------------------------------------------
    # Accident mortality
    # -----------------------------------------------------
    
    accident_death = case_when(
      mortstat == 1 & ucod_leading == 4 ~ 1,
      mortstat == 1 & ucod_leading != 4 ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    
    # -----------------------------------------------------
    # Alzheimer disease mortality
    # -----------------------------------------------------
    
    alz_death = case_when(
      mortstat == 1 & ucod_leading == 6 ~ 1,
      mortstat == 1 & ucod_leading != 6 ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    
    # -----------------------------------------------------
    # Kidney disease mortality
    # -----------------------------------------------------
    
    neph_death = case_when(
      mortstat == 1 & ucod_leading == 9 ~ 1,
      mortstat == 1 & ucod_leading != 9 ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    
    # -----------------------------------------------------
    # Other mortality
    # -----------------------------------------------------
    
    other_death = case_when(
      mortstat == 1 & ucod_leading == 10 ~ 1,
      mortstat == 1 & ucod_leading != 10 ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  
  filter(
    !is.na(age_entry),
    !is.na(age_exit),
    !is.na(all_cause_death),
    age_exit > age_entry
  )

# ---------------------------------------------------------
# NHANES survey design
# ---------------------------------------------------------

options(survey.lonely.psu = "adjust")

nhanes_surv_design <- svydesign(
  ids = ~SDPPSU6,
  strata = ~SDPSTRA6,
  weights = ~WTPFEX6,
  nest = TRUE,
  data = nhanes_surv
)

# ---------------------------------------------------------
# Sample size and event counts
# ---------------------------------------------------------

event_counts <- nhanes_surv %>%
  group_by(.data[[EXPOSURE]]) %>%
  summarise(
    
    n = n(),
    
    all_cause =
      sum(all_cause_death == 1, na.rm = TRUE),
    
    cardiometabolic =
      sum(cardiometabolic_death == 1, na.rm = TRUE),
    
    respiratory =
      sum(respiratory_death == 1, na.rm = TRUE),
    
    cancer =
      sum(cancer_death == 1, na.rm = TRUE),
    
    accident =
      sum(accident_death == 1, na.rm = TRUE),
    
    alzheimer =
      sum(alz_death == 1, na.rm = TRUE),
    
    kidney =
      sum(neph_death == 1, na.rm = TRUE),
    
    other =
      sum(other_death == 1, na.rm = TRUE),
    
    .groups = "drop"
  )

print(event_counts, n = Inf)

# ---------------------------------------------------------
# Helper function:
# Fit Model 1 and Model 2 for one mortality outcome
# ---------------------------------------------------------

fit_svy_cox <- function(outcome, design, exposure) {
  
  # -------------------------------------------------------
  # Model 1
  # -------------------------------------------------------
  
  form_m1 <- as.formula(
    paste0(
      "Surv(age_entry, age_exit, ", outcome, ") ~ ",
      exposure,
      " + sex + race_ethnicity"
    )
  )
  
  
  # -------------------------------------------------------
  # Model 2
  # -------------------------------------------------------
  
  form_m2 <- as.formula(
    paste0(
      "Surv(age_entry, age_exit, ", outcome, ") ~ ",
      exposure,
      " + sex + race_ethnicity + education_lt12 + ",
      "smoking_status + DRINKS2 + sedentary"
    )
  )
  
  
  # -------------------------------------------------------
  # Fit models
  # -------------------------------------------------------
  
  list(
    
    model1 = svycoxph(
      form_m1,
      design = design
    ),
    
    model2 = svycoxph(
      form_m2,
      design = design
    )
  )
}



# ---------------------------------------------------------
# Outcomes
# ---------------------------------------------------------

outcomes <- c(
  "all_cause_death",
  "cardiometabolic_death",
  "respiratory_death",
  "cancer_death",
  "accident_death",
  "alz_death",
  "neph_death",
  "other_death"
)



# ---------------------------------------------------------
# Fit Cox models
# ---------------------------------------------------------

cox_models <- setNames(
  lapply(
    outcomes,
    fit_svy_cox,
    design = nhanes_surv_design,
    exposure = EXPOSURE
  ),
  outcomes
)




# ---------------------------------------------------------
# Model summaries
# ---------------------------------------------------------

summary(cox_models$all_cause_death$model1)
summary(cox_models$all_cause_death$model2)


summary(cox_models$cardiometabolic_death$model1)
summary(cox_models$cardiometabolic_death$model2)


summary(cox_models$respiratory_death$model1)
summary(cox_models$respiratory_death$model2)


summary(cox_models$cancer_death$model1)
summary(cox_models$cancer_death$model2)


summary(cox_models$accident_death$model1)
summary(cox_models$accident_death$model2)


summary(cox_models$alz_death$model1)
summary(cox_models$alz_death$model2)


summary(cox_models$neph_death$model1)
summary(cox_models$neph_death$model2)


summary(cox_models$other_death$model1)
summary(cox_models$other_death$model2)




