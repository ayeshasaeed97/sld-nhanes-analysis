# survival analysis using Cox regression

# competing risk (especially between cmd and alz)

# loading packages and data
library(tidyverse)
library(survey)
library(survival)

nhanes_table1 <- read_csv("clean_data/masld_table1.csv")

nhanes_table1 <- nhanes_table1 %>%
  mutate(
    LIVER_GROUP3 = factor(
      case_when(
      LIVER_GROUP2 == "Metabolic dysfunction–associated steatotic liver disease" | 
        LIVER_GROUP2 == "Metabolic dysfunction–associated steatohepatitis" ~ "MASLD",
      LIVER_GROUP2 == "No hepatic steatosis" ~ "None",
      TRUE ~ NA_character_
    ),
    levels = c("None", "MASLD")
  )
  ) %>%
  filter(egfr > 15)

table(nhanes_table1$LIVER_GROUP3)
table(nhanes_table1$LIVER_GROUP2)

# Cox regression/Survival ----

nhanes_surv <- nhanes_table1 %>%
  mutate(
    age_entry = HSAGEIR,
    age_exit = HSAGEIR + permth_exm / 12,
    all_cause_death = if_else(mortstat == 1, 1, 0)) %>%
  
  mutate(
    cancer_death = case_when(
      mortstat == 1 & ucod_leading == 2 ~ 1,
      mortstat == 1 & !(ucod_leading == 2) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  
  mutate(
    cvd_death = case_when(
      mortstat == 1 & ucod_leading == 1 ~ 1,
      mortstat == 1 & !(ucod_leading == 1) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  
  mutate(
    resp_death = case_when(
      mortstat == 1 & ucod_leading == 3 ~ 1,
      mortstat ==1 & !(ucod_leading == 3) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  
  mutate(
    acdnt_death = case_when(
      mortstat == 1 & ucod_leading == 4 ~ 1,
      mortstat == 1 & !(ucod_leading == 4) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  
  mutate(
    cereb_death = case_when(
      mortstat == 1 & ucod_leading == 5 ~ 1,
      mortstat == 1 & !(ucod_leading == 5) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  
  mutate(
    alz_death = case_when(
      mortstat == 1 & ucod_leading == 6 ~ 1,
      mortstat == 1 & !(ucod_leading == 6) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    ) 
  ) %>%
  
  mutate(
    diab_death = case_when(
      mortstat == 1 & ucod_leading == 7 ~ 1,
      mortstat == 1 & !(ucod_leading == 7) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    ) 
  ) %>%
  
  mutate(
    flu_death = case_when(
      mortstat == 1 & ucod_leading == 8 ~ 1,
      mortstat == 1 & !(ucod_leading == 8) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    ) 
  ) %>%
  
  mutate(
    neph_death = case_when(
      mortstat == 1 & ucod_leading == 9 ~ 1,
      mortstat == 1 & !(ucod_leading == 9) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    other_death = case_when(
      mortstat == 1 & ucod_leading == 10 ~ 1,
      mortstat == 1 & !(ucod_leading == 10) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(LIVER_GROUP3),
    !is.na(age_entry),
    !is.na(age_exit),
    !is.na(all_cause_death)
  )

# nhanes_surv recoding mortality groups-----

nhanes_surv2 <- nhanes_table1 %>%
  mutate(
    age_entry = HSAGEIR,
    age_exit = HSAGEIR + permth_exm / 12,
    all_cause_death = if_else(mortstat == 1, 1, 0)) %>%
  
  # Heart disease + cerebrovascular + diabetes -> cardiometabolic mortality
  mutate(
    cardiometabolic_death = case_when(
      mortstat == 1 & ucod_leading %in% c(1, 5, 7) ~ 1,
      mortstat == 1 & !(ucod_leading %in% c(1, 5, 7)) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  
  # chronic lower respiratory disease + influenza/pneumonia -> respiratory mortality
  mutate(
    respiratory_death = case_when(
      mortstat == 1 & ucod_leading %in% c(3, 8) ~ 1,
      mortstat == 1 & !(ucod_leading %in% c(3, 8)) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    other_death = case_when(
      mortstat == 1 & ucod_leading == 10 ~ 1,
      mortstat == 1 & !(ucod_leading == 10) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  
  mutate(
    neph_death = case_when(
      mortstat == 1 & ucod_leading == 9 ~ 1,
      mortstat == 1 & !(ucod_leading == 9) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  
  mutate(
    acdnt_death = case_when(
      mortstat == 1 & ucod_leading == 4 ~ 1,
      mortstat == 1 & !(ucod_leading == 4) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>% mutate(
    cancer_death = case_when(
      mortstat == 1 & ucod_leading == 2 ~ 1,
      mortstat == 1 & !(ucod_leading == 2) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  
  mutate(
    alz_death = case_when(
      mortstat == 1 & ucod_leading == 6 ~ 1,
      mortstat == 1 & !(ucod_leading == 6) ~ 0,
      mortstat == 0 ~ 0,
      TRUE ~ NA_real_
    ) 
  ) %>%
  
  filter(
    !is.na(LIVER_GROUP3),
    !is.na(age_entry),
    !is.na(age_exit),
    !is.na(all_cause_death),
    age_exit > age_entry
  )

options(survey.lonely.psu = "adjust") # to center the stratum at the population mean rather than the stratum mean,

nhanes_surv_design <- svydesign(
  ids = ~SDPPSU6,
  strata = ~SDPSTRA6,
  weights = ~WTPFEX6,
  nest = TRUE,
  data = nhanes_surv2
)






## All cause mortality ----
nhanes_surv2 %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    all_events = sum(all_cause_death == 1, na.rm = TRUE),
    all_event_percent = mean(all_cause_death == 1, na.rm = TRUE) * 100
  )

cox_all_m1 <- svycoxph(
  Surv(age_entry, age_exit, all_cause_death) ~
    LIVER_GROUP3 + sex + race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_all_m1)


cox_all_m2 <- svycoxph(
  Surv(age_entry, age_exit, all_cause_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 + # replace with grams per day
    sedentary, # sedentary lifestyle is on the causal pathway
  
  design = nhanes_surv_design
)

summary(cox_all_m2)

cox_all_m3 <- svycoxph(
  Surv(age_entry, age_exit, all_cause_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    diabetes +
    egfr,
  
  design = nhanes_surv_design
)

summary(cox_all_m3)


# Cardiometabolic mortality

nhanes_surv2 %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    cardiometabolic_events = sum(cardiometabolic_death == 1, na.rm = TRUE),
    cardiometabolic_event_percent = mean(cardiometabolic_death == 1, na.rm = TRUE) * 100
  )

cox_cmd_m1 <- svycoxph(
  Surv(age_entry, age_exit, cardiometabolic_death) ~
    LIVER_GROUP3 + sex + race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_cmd_m1)

cox_cmd_m2 <- svycoxph(
  Surv(age_entry, age_exit, cardiometabolic_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_cmd_m2)


cox_cmd_m3 <- svycoxph(
  Surv(age_entry, age_exit, cardiometabolic_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia,
  design = nhanes_surv_design
)

summary(cox_cmd_m3)


## respiratory mortality ----


nhanes_surv2 %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    respiratory_death_events = sum(respiratory_death == 1, na.rm = TRUE),
    respiratory_death_event_percent = mean(respiratory_death == 1, na.rm = TRUE) * 100
  )

cox_cmd_m1 <- svycoxph(
  Surv(age_entry, age_exit, respiratory_death) ~
    LIVER_GROUP3 + sex + race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_cmd_m1)


## Cardiovascular disease mortality ----

nhanes_surv2 %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    cvd_events = sum(cvd_death == 1, na.rm = TRUE),
    cvd_event_percent = mean(cvd_death == 1, na.rm = TRUE) * 100
  )

cox_cvd_m1 <- svycoxph(
  Surv(age_entry, age_exit, cvd_death) ~
    LIVER_GROUP3 + sex + race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_cvd_m1)

cox_cvd_m2 <- svycoxph(
  Surv(age_entry, age_exit, cvd_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_cvd_m2)




## Cancer Death ----

nhanes_surv %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    cancer_events = sum(cancer_death == 1, na.rm = TRUE),
    cancer_event_percent = mean(cancer_death == 1, na.rm = TRUE) * 100
  )


cox_cancer_m1 <- svycoxph(
  Surv(age_entry, age_exit, cancer_death) ~
    LIVER_GROUP3 + sex + race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_cancer_m1)

cox_cancer_m2 <- svycoxph(
  Surv(age_entry, age_exit, cancer_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_cancer_m2)



## Chronic lower respiratory diseases ----

nhanes_surv %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    resp_events = sum(resp_death == 1, na.rm = TRUE),
    resp_event_percent = mean(resp_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_resp_m1 <- svycoxph(
  Surv(age_entry, age_exit, resp_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_resp_m1)

cox_resp_m2 <- svycoxph(
  Surv(age_entry, age_exit, resp_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_resp_m2)


## Accident mortality ----

nhanes_surv %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    acdnt_events = sum(acdnt_death == 1, na.rm = TRUE),
    acdnt_event_percent = mean(acdnt_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_acdnt_m1 <- svycoxph(
  Surv(age_entry, age_exit, acdnt_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_acdnt_m1)

cox_acdnt_m2 <- svycoxph(
  Surv(age_entry, age_exit, acdnt_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_acdnt_m2)


## Cerebrovascular disease mortality ----

nhanes_surv %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    cereb_events = sum(cereb_death == 1, na.rm = TRUE),
    cereb_event_percent = mean(cereb_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_cereb_m1 <- svycoxph(
  Surv(age_entry, age_exit, cereb_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_cereb_m1)

cox_cereb_m2 <- svycoxph(
  Surv(age_entry, age_exit, cereb_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_cereb_m2)


## Alzheimer's disease mortality ----

nhanes_surv2 <- nhanes_surv %>%
  mutate(
    LIVER_GROUP3 = factor(
      case_when(
        LIVER_GROUP3 == "Non-alcoholic steatohepatitis" ~ "NASH/NAFLD",
        LIVER_GROUP3 == "Non-alcoholic fatty liver disease" ~ "NASH/NAFLD",
        LIVER_GROUP3 == "No hepatic steatosis" ~ "Neither",
        TRUE ~ NA_character_
      ),
      levels = c("Neither", "NASH/NAFLD")
    )
  )

options(survey.lonely.psu = "adjust")

nhanes_surv_design2 <- svydesign(
  ids = ~SDPPSU6,
  strata = ~SDPSTRA6,
  weights = ~WTPFEX6,
  nest = TRUE,
  data = nhanes_surv2
)


nhanes_surv2 %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    alz_events = sum(alz_death == 1, na.rm = TRUE),
    alz_event_percent = mean(alz_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_alz_m1 <- svycoxph(
  Surv(age_entry, age_exit, alz_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity,
  design = nhanes_surv_design2
)

summary(cox_alz_m1)

cox_alz_m2 <- svycoxph(
  Surv(age_entry, age_exit, alz_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design2
)

summary(cox_alz_m2)


## Diabetes mortality ----

nhanes_surv %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    diab_events = sum(diab_death == 1, na.rm = TRUE),
    diab_event_percent = mean(diab_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_diab_m1 <- svycoxph(
  Surv(age_entry, age_exit, diab_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_diab_m1)

cox_diab_m2 <- svycoxph(
  Surv(age_entry, age_exit, diab_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_diab_m2)

cox_diab_m3 <- svycoxph(
  Surv(age_entry, age_exit, diab_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia,
  design = nhanes_surv_design
)

summary(cox_diab_m3)

## Influenza and pneumonia mortality ----

nhanes_surv %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    flu_events = sum(flu_death == 1, na.rm = TRUE),
    flu_event_percent = mean(flu_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_flu_m1 <- svycoxph(
  Surv(age_entry, age_exit, flu_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_flu_m1)

cox_flu_m2 <- svycoxph(
  Surv(age_entry, age_exit, flu_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_flu_m2)


## Nephritis, nephrotic syndrome, and nephrosis mortality ----

nhanes_surv %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    neph_events = sum(neph_death == 1, na.rm = TRUE),
    neph_event_percent = mean(neph_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_neph_m1 <- svycoxph(
  Surv(age_entry, age_exit, neph_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_neph_m1)

cox_neph_m2 <- svycoxph(
  Surv(age_entry, age_exit, neph_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_neph_m2)


## Other-cause mortality ----

nhanes_surv %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    other_events = sum(other_death == 1, na.rm = TRUE),
    other_event_percent = mean(other_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_other_m1 <- svycoxph(
  Surv(age_entry, age_exit, other_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_other_m1)

cox_other_m2 <- svycoxph(
  Surv(age_entry, age_exit, other_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_other_m2)



## PH Assumption ----


## Proportional-hazards test: Model 1
ph_all_m1 <- cox.zph(
  cox_all_m1,
  transform = "km",
  terms = TRUE,
  global = TRUE
)

print(ph_all_m1)


## Proportional-hazards test: Model 2
ph_all_m2 <- cox.zph(
  cox_all_m2,
  transform = "km",
  terms = TRUE,
  global = TRUE
)

print(ph_all_m2)

## schoenfeld
plot(ph_all_m1)

plot(ph_all_m2)
















nhanes_surv2 %>%
  group_by(LIVER_GROUP3) %>%
  summarise(
    n = n(),
    all_events = sum(other_death == 1, na.rm = TRUE),
    all_event_percent = mean(other_death == 1, na.rm = TRUE) * 100
  )

cox_all_m1 <- svycoxph(
  Surv(age_entry, age_exit, other_death) ~
    LIVER_GROUP3 + sex + race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_all_m1)


cox_all_m2 <- svycoxph(
  Surv(age_entry, age_exit, other_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 + # replace with grams per day
    sedentary, # sedentary lifestyle is on the causal pathway
  
  design = nhanes_surv_design
)

summary(cox_all_m2)

cox_all_m3 <- svycoxph(
  Surv(age_entry, age_exit, other_death) ~
    LIVER_GROUP3 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    DRINKS2 +
    sedentary +
    BMPBMI +
    hypertension +
    diabetes +
    egfr,
  
  design = nhanes_surv_design
)

summary(cox_all_m3)


