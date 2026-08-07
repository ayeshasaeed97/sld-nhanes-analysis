##############################################
# Prevent score calculation
##############################################

library(tidyverse)
library(pillar)
library(skimr)
library(preventr)
library(haven)

# Load Data----
masld <- read_csv("clean_data/nhanes_masld.csv")


masld2 <- masld %>%
  mutate(sex = tolower(trimws(SEX)),
         
         statin = trimws(STATIN_MED),
         
         statin = if_else(statin=="True", TRUE, FALSE, missing = FALSE),
         
         bp_tx = if_else(HBP_MED=="TRUE", TRUE, FALSE, missing = FALSE),
         
         # Cap age to the range 30–79
         age2 = pmin(pmax(HSAGEIR, 30), 79),
         
         # Cap BMI to the range 18.5–39.9
         bmi2 = pmin(pmax(BMPBMI, 18.5), 39.9),
         
         # Cap systolic blood pressure to the range 90–180
         sbp2 = pmin(pmax(PEPMNK1R, 90), 180),
         
         # Cap total cholesterol to the range 3.36–8.28 mmol/L
         total_c2 = pmin(pmax(TCPSI, 3.36), 8.28),
         
         # Cap hdl cholesterol to range 0.5202.59 mmol/L
         hdl_c2 = pmin(pmax(HDPSI, 2.59), 0.52),
         
         # Convert serum creatinine from µmol/L to mg/dL
         scr_mg_dl = CEPSI / 88.4,
         
         # Sex-specific CKD-EPI constants
         kappa = case_when(
           sex == "f" ~ 0.7,
           sex == "m" ~ 0.9,
           TRUE ~ NA_real_
         ),
         
         alpha = case_when(
           sex == "f" ~ -0.241,
           sex == "m" ~ -0.302,
           TRUE ~ NA_real_
         ),
         
         # creatinine equation
         egfr = 142 *
           pmin(scr_mg_dl / kappa, 1)^alpha *
           pmax(scr_mg_dl / kappa, 1)^(-1.200) *
           0.9938^HSAGEIR *
           if_else(sex == "f", 1.012, 1),
         
         egfr2 = case_when(
           
           # Female, SCr <= 0.7
           sex == "f" & scr_mg_dl <= 0.7 ~
             142 *
             (scr_mg_dl / 0.7)^(-0.241) *
             (0.9938^HSAGEIR) *
             1.012,
           
           # Female, SCr > 0.7
           sex == "f" & scr_mg_dl > 0.7 ~
             142 *
             (scr_mg_dl / 0.7)^(-1.200) *
             (0.9938^HSAGEIR) *
             1.012,
           
           # Male, SCr <= 0.9
           sex == "m" & scr_mg_dl <= 0.9 ~
             142 *
             (scr_mg_dl / 0.9)^(-0.302) *
             (0.9938^HSAGEIR),
           
           # Male, SCr > 0.9
           sex == "m" & scr_mg_dl > 0.9 ~
             142 *
             (scr_mg_dl / 0.9)^(-1.200) *
             (0.9938^HSAGEIR),
           
           # Missing/invalid combinations
           TRUE ~ NA_real_
           )
  )
  
hist(masld2$egfr)
summary(masld2$egfr)

# missing was assigned to false

# Score calculation

prevent_score <- masld2 %>%
  rowwise() %>%
  mutate(
    risk_score = list(
      estimate_risk(
        age = HSAGEIR,
        sex = sex,
        bp_tx = bp_tx,
        sbp = PEPMNK1R,
        total_c = TCPSI,
        hdl_c = HDPSI,
        bmi = BMPBMI,
        statin = statin,
        dm = DIABETES,
        smoking = SMOKING,
        
        # Use the current row's creatinine value
        egfr = call(
          "calc_egfr",
          cr = CEPSI,
          units = "umol"
        ),
        
        chol_unit = "mmol/L",
        time = "both",
        collapse = TRUE
      )
    )
  ) %>%
  ungroup()


# capped version of PREVENT score
masld3 <- masld2 %>%
  drop_na(HDPSI,
          TCPSI,
          DIABETES,
          SMOKING,
          sbp2)


prevent_score_cap <- masld3 %>%
  rowwise() %>%
  mutate(
    risk_score = list(
      estimate_risk(
        age = age2,
        sex = sex,
        bp_tx = bp_tx,
        sbp = sbp2,
        total_c = total_c2,
        hdl_c = hdl_c2,
        bmi = bmi2,
        statin = statin,
        dm = DIABETES,
        smoking = SMOKING,
        egfr = egfr,
        chol_unit = "mmol/L",
        time = "both",
        collapse = TRUE,
        quiet = TRUE
      )
    )
  ) %>%
  ungroup()


prevent_results <- prevent_score %>%
  unnest(risk_score)

prevent_results_cap <- prevent_score_cap %>%
  unnest(risk_score)

summary(prevent_results_cap)
