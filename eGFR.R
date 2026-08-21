# checking creatinine scores
library(tidyverse)

masld <- read_csv("clean_data/nhanes_masld.csv")

# egfr correction
masld2 <- masld %>%
  mutate(
    CEP = CEPSI/88.4,
    creatinine = 0.960*CEP - 0.184
  )


masld2 <- masld2 %>%
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
         
         
         # eGFR equations
         egfr = case_when(
           
           # Female, SCr <= 0.7
           sex == "f" & creatinine <= 0.7 ~
             142 *
             (creatinine / 0.7)^(-0.241) *
             (0.9938^HSAGEIR) *
             1.012,
           
           # Female, SCr > 0.7
           sex == "f" & creatinine > 0.7 ~
             142 *
             (creatinine / 0.7)^(-1.200) *
             (0.9938^HSAGEIR) *
             1.012,
           
           # Male, SCr <= 0.9
           sex == "m" & creatinine <= 0.9 ~
             142 *
             (creatinine / 0.9)^(-0.302) *
             (0.9938^HSAGEIR),
           
           # Male, SCr > 0.9
           sex == "m" & creatinine > 0.9 ~
             142 *
             (creatinine / 0.9)^(-1.200) *
             (0.9938^HSAGEIR),
           
           # Missing/invalid combinations
           TRUE ~ NA_real_
         )
  )

# 

# female with serum creatinine <= 0.7

female1 <- masld2 %>%
  filter(HSSEX == 2 & creatinine <= 0.7) ## 77 observations

egfr1 <- (142*(female1$creatinine/0.7)^-0.241)*((0.9938)^(female1$HSAGEIR))*1.012
summary(egfr1)

# female with serum creatinine > 0.7
female2 <- masld2 %>%
  filter(HSSEX == 2 & creatinine > 0.7) ## 6703 observations

egfr2 <- (142 * (female2$creatinine / 0.7)^-1.200) *((0.9938)^(female2$HSAGEIR))*1.012
summary(egfr2)


# male with serum creatinine <= 0.9
male1 <- masld2 %>%
  filter(HSSEX == 1 & creatinine <= 0.9) ## 84 observations

egfr3 <- (142 * (male1$creatinine / 0.9)^-0.302) * ((0.9938)^(male1$HSAGEIR))
summary(egfr3)

# male with serum creatinine > 0.9
male2 <- masld2 %>%
  filter(HSSEX == 1 & creatinine > 0.9) ## 5501 observations

egfr4 <- (142 * (male2$creatinine / 0.9)^-1.200) * ((0.9938)^(male2$HSAGEIR))
summary(egfr4)

egfr_wh <- c(egfr1, egfr2, egfr3, egfr4)
summary(egfr_wh)


summary(masld2$egfr)



