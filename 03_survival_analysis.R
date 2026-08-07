# survival analysis using Cox regression

# loading packages and data
library(tidyverse)
library(survey)
library(survival)

nhanes_table1 <- read.csv("clean_data/nhanes_table1.csv")
 
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
    !is.na(LIVER_GROUP),
    !is.na(age_entry),
    !is.na(age_exit),
    !is.na(all_cause_death)
  )



nhanes_surv <- nhanes_surv %>%
  filter(
    !is.na(age_entry),
    !is.na(age_exit),
    age_exit > age_entry
  )


options(survey.lonely.psu = "adjust")

nhanes_surv_design <- svydesign(
  ids = ~SDPPSU6,
  strata = ~SDPSTRA6,
  weights = ~WTPFEX6,
  nest = TRUE,
  data = nhanes_surv
)


## All cause mortality ----
cox_all_m1 <- svycoxph(
  Surv(age_entry, age_exit, all_cause_death) ~
    LIVER_GROUP + sex + race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_all_m1)

cox_all_m2 <- svycoxph(
  Surv(age_entry, age_exit, all_cause_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_all_m2)


nhanes_surv %>%
  group_by(LIVER_GROUP) %>%
  summarise(
    n = n(),
    events = sum(all_cause_death == 1, na.rm = TRUE),
    event_percent = mean(all_cause_death == 1, na.rm = TRUE) * 100
  )


cox_all_m3 <- svycoxph(
  Surv(age_entry, age_exit, all_cause_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
    sedentary +
    BMPBMI +
    hypertension +
    diabetes,
  design = nhanes_surv_design
)

summary(cox_all_m3)


## Cardiovascular disease mortality ----

nhanes_surv %>%
  group_by(LIVER_GROUP) %>%
  summarise(
    n = n(),
    cvd_events = sum(cvd_death == 1, na.rm = TRUE),
    cvd_event_percent = mean(cvd_death == 1, na.rm = TRUE) * 100
  )

cox_cvd_m1 <- svycoxph(
  Surv(age_entry, age_exit, cvd_death) ~
    LIVER_GROUP + sex + race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_cvd_m1)

cox_cvd_m2 <- svycoxph(
  Surv(age_entry, age_exit, cvd_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
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
  group_by(LIVER_GROUP) %>%
  summarise(
    n = n(),
    cancer_events = sum(cancer_death == 1, na.rm = TRUE),
    cancer_event_percent = mean(cancer_death == 1, na.rm = TRUE) * 100
  )


cox_cancer_m1 <- svycoxph(
  Surv(age_entry, age_exit, cancer_death) ~
    LIVER_GROUP + sex + race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_cancer_m1)

cox_cancer_m2 <- svycoxph(
  Surv(age_entry, age_exit, cancer_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
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
  group_by(LIVER_GROUP) %>%
  summarise(
    n = n(),
    resp_events = sum(resp_death == 1, na.rm = TRUE),
    resp_event_percent = mean(resp_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_resp_m1 <- svycoxph(
  Surv(age_entry, age_exit, resp_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_resp_m1)

cox_resp_m2 <- svycoxph(
  Surv(age_entry, age_exit, resp_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
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
  group_by(LIVER_GROUP) %>%
  summarise(
    n = n(),
    acdnt_events = sum(acdnt_death == 1, na.rm = TRUE),
    acdnt_event_percent = mean(acdnt_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_acdnt_m1 <- svycoxph(
  Surv(age_entry, age_exit, acdnt_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_acdnt_m1)

cox_acdnt_m2 <- svycoxph(
  Surv(age_entry, age_exit, acdnt_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
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
  group_by(LIVER_GROUP) %>%
  summarise(
    n = n(),
    cereb_events = sum(cereb_death == 1, na.rm = TRUE),
    cereb_event_percent = mean(cereb_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_cereb_m1 <- svycoxph(
  Surv(age_entry, age_exit, cereb_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_cereb_m1)

cox_cereb_m2 <- svycoxph(
  Surv(age_entry, age_exit, cereb_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
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
    liver_group2 = factor(
      case_when(
      LIVER_GROUP == "Non-alcoholic steatohepatitis" ~ "NASH/NAFLD",
      LIVER_GROUP == "Non-alcoholic fatty liver disease" ~ "NASH/NAFLD",
      LIVER_GROUP == "No hepatic steatosis" ~ "Neither",
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
  group_by(liver_group2) %>%
  summarise(
    n = n(),
    alz_events = sum(alz_death == 1, na.rm = TRUE),
    alz_event_percent = mean(alz_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_alz_m1 <- svycoxph(
  Surv(age_entry, age_exit, alz_death) ~
    liver_group2 +
    sex +
    race_ethnicity,
  design = nhanes_surv_design2
)

summary(cox_alz_m1)

cox_alz_m2 <- svycoxph(
  Surv(age_entry, age_exit, alz_death) ~
    liver_group2 +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
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
  group_by(LIVER_GROUP) %>%
  summarise(
    n = n(),
    diab_events = sum(diab_death == 1, na.rm = TRUE),
    diab_event_percent = mean(diab_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_diab_m1 <- svycoxph(
  Surv(age_entry, age_exit, diab_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_diab_m1)

cox_diab_m2 <- svycoxph(
  Surv(age_entry, age_exit, diab_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
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
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
    sedentary +
    BMPBMI +
    hypertension +
    hypercholesterolaemia,
  design = nhanes_surv_design
)

summary(cox_diab_m3)

## Influenza and pneumonia mortality ----

nhanes_surv %>%
  group_by(LIVER_GROUP) %>%
  summarise(
    n = n(),
    flu_events = sum(flu_death == 1, na.rm = TRUE),
    flu_event_percent = mean(flu_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_flu_m1 <- svycoxph(
  Surv(age_entry, age_exit, flu_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_flu_m1)

cox_flu_m2 <- svycoxph(
  Surv(age_entry, age_exit, flu_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
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
  group_by(LIVER_GROUP) %>%
  summarise(
    n = n(),
    neph_events = sum(neph_death == 1, na.rm = TRUE),
    neph_event_percent = mean(neph_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_neph_m1 <- svycoxph(
  Surv(age_entry, age_exit, neph_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_neph_m1)

cox_neph_m2 <- svycoxph(
  Surv(age_entry, age_exit, neph_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
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
  group_by(LIVER_GROUP) %>%
  summarise(
    n = n(),
    other_events = sum(other_death == 1, na.rm = TRUE),
    other_event_percent = mean(other_death == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cox_other_m1 <- svycoxph(
  Surv(age_entry, age_exit, other_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity,
  design = nhanes_surv_design
)

summary(cox_other_m1)

cox_other_m2 <- svycoxph(
  Surv(age_entry, age_exit, other_death) ~
    LIVER_GROUP +
    sex +
    race_ethnicity +
    education_lt12 +
    smoking_status +
    alcohol_category +
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