setwd("~/Documents/FL_NHANES")

# loading packages ----
library(tidyverse)
library(flextable)
library(haven) # for reading xpt files
library(readr)

# loading data -------------------------------
adult <- read.csv('raw_data/adult_output.csv')
lab <- read.csv('raw_data/lab_output.csv')
exam <- read.csv('raw_data/exam_output.csv')
hep <- read_xpt('raw_data/HGUHS.xpt')
meds <- read_xpt('raw_data/pupremed.xpt')

srvyin <- paste('raw_data/NHANES_III_MORT_2019_PUBLIC.dat')   # full .DAT name here
srvyout <- "<SURVEY>" 

# read in the fixed-width format ASCII file
dsn <- read_fwf(file=srvyin,
                col_types = "iiiiiiii",
                fwf_cols(SEQN = c(1,6),
                         eligstat = c(15,15),
                         mortstat = c(16,16),
                         ucod_leading = c(17,19),
                         diabetes = c(20,20),
                         hyperten = c(21,21),
                         permth_int = c(43,45),
                         permth_exm = c(46,48)
                ),
                na = c("", ".")
)

# Structure and contents of data
str(dsn) # mortality data


# Variable Selection ----------------------------------

# lab
lab_vars <- c(
  "ATPSI", "ASPSI", "GGPSI", "HCP", "HBP",
  "FEPSI", "FRP", "FRPSI", "TIPSI", "PXP",
  "GHP", "TGPSI", "HDPSI", "I1PSI",
  "TCPSI", "SGPSI", "PHPFAST", "PLP", "CEPSI"
) # added CEP july 31st for prevent score

lab_clean <- lab %>%
  select(SEQN, all_of(lab_vars)) %>%
  mutate(
    across(
      all_of(lab_vars),
      ~ na_if(.x, 8888)
    )
  ) %>%
  mutate(
    across(
      all_of(lab_vars),
      ~ na_if(.x, 888)
    )
  ) %>%
  mutate(
    TFR = if_else(
      !is.na(FEPSI) & !is.na(TIPSI) & TIPSI > 0,
      (FEPSI / TIPSI) * 100,
      NA_real_
    )
  )

# exam
exam_clean <- exam %>%
  select(
    SEQN, BMPHT, BMPWTLBS, BMPWAIST, MAPE1,
    PEPMNK1R, PEPMNK5R, BMPBMI, MAPE3S, MAPE4
  ) %>%
  mutate(
    MAPE3S = na_if(MAPE3S, 888),
    MAPE3S = na_if(MAPE3S, 999),
    MAPE4  = na_if(MAPE4, 888),
    MAPE4  = na_if(MAPE4, 999),
    BMPBMI = na_if(BMPBMI, 8888),
    DRINK = (MAPE3S*MAPE4)/365)


adult_clean <- adult %>%
  select(SEQN, HSSEX, HSAGEIR, DMARACER, DMAETHNR, DMARETHN, HAD1,  HAF10, HAC1D, HAC1C, HAC1N, 
         HAC1O, HAT2, HAT4, HAT6, HAT8, HAT10,
         HAT12, HAT14, HAT16, HAT18, HAE7, HAR3, HAR1, HFA8R, HAE2,
         WTPFEX6, SDPSTRA6, SDPPSU6) %>%
  mutate(
    HFA8R = na_if(HFA8R, 88),
    HFA8R = na_if(HFA8R, 99)
  )


hep_med_codes <- c(0212, 3060, 2306, 1770, 1319, 1032)

meds_clean <- meds %>%
  select(SEQN, HQRXCODE) %>%
  group_by(SEQN) %>%
  summarise(
    HEP_MED = factor(
      if_else(
        any(HQRXCODE %in% hep_med_codes, na.rm = TRUE),
        "Yes",
        "No"
      ),
      levels = c("No", "Yes")
    ),
    .groups = "drop"
  )

# merging data by SEQN number
nhanes_merged <- adult_clean %>%
  left_join(exam_clean, by = "SEQN") %>%
  left_join(lab_clean,  by = "SEQN") %>% 
  left_join(hep, by = "SEQN") %>%
  left_join(dsn, by = "SEQN") %>%
  left_join(meds_clean, by = "SEQN") %>%
  mutate(
    HEP_MED = as.character(HEP_MED),
    HEP_MED = if_else(is.na(HEP_MED), "No", HEP_MED),
    HEP_MED = factor(HEP_MED, levels = c("No", "Yes"))
  )

# prior to checking missingness, drop N/A values from essential variables
nhanes_merged <- nhanes_merged %>%
  drop_na(GUPHSPF,
          BMPBMI,
          HBP,
          HCP,
          ATPSI,
          ASPSI,
          MAPE1,
          PEPMNK1R,
          PEPMNK5R,
          HAE2,
          mortstat,
          FRPSI
  ) 


nhanes_merged <- nhanes_merged %>%
  drop_na(GUPHSPF)


## evaluating missingess --------------------
nhanes_missing <- nhanes_merged %>%
  summarise(
    across(
      everything(),
      list(
        n_missing = ~sum(is.na(.)),
        pct_missing = ~ mean(is.na(.)) * 100
      )
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_(n_missing|pct_missing)$"
  ) %>%
  arrange(desc(pct_missing))

print(nhanes_missing, n= 100)


# NHANES Working Data----
nhanes_work <- nhanes_merged %>%
  mutate(
    ALC = factor(
      case_when(
        MAPE1 == 2 ~ "No alcohol consumption",
        HSSEX == 2 & DRINK >= 1 ~ "High alcohol consumption",
        HSSEX == 1 & DRINK >= 2 ~ "High alcohol consumption",
        TRUE ~ "Low-Moderate alcohol consumption"
      ),
      levels = c(
        "No alcohol consumption",
        "Low-Moderate alcohol consumption",
        "High alcohol consumption"
      )
    )
  ) %>%
  mutate(
    LIVER_ENZ = factor(
      case_when(
        HSSEX == 1 & (ATPSI > 40 | ASPSI > 37) ~ "Raised",
        HSSEX == 2 & (ATPSI > 31 | ASPSI > 31) ~ "Raised",
        HSSEX == 1 & ATPSI <= 40 & ASPSI <= 37 ~ "Normal",
        HSSEX == 2 & ATPSI <= 31 & ASPSI <= 31 ~ "Normal",
        TRUE ~ NA_character_
      ),
      levels = c("Normal", "Raised")
    )
  ) %>%
  mutate(
    IRON_OVERLOAD = factor(
      case_when(
        HSSEX == 1 & TFR > 45 & FRPSI > 500 ~ "Iron overload",
        HSSEX == 2 & TFR > 45 & FRPSI > 400 ~ "Iron overload",
        !is.na(TFR) & !is.na(FRPSI) ~ "No iron overload",
        TRUE ~ NA_character_
      ),
      levels = c("No iron overload", "Iron overload")
    )
  )%>%
  mutate(
    SED = factor(
      if_else(
        if_all(c(HAT2, HAT4, HAT6, HAT8, HAT10, HAT12, HAT14, HAT16, HAT18), ~ .x == 2),
        "Sedentary",
        "Not sedentary"
      ),
      levels = c("Not sedentary", "Sedentary")
    )) %>% 
  mutate(
    SMOKE = factor(
      case_when(
        HAR1 == 2 ~ "Never",
        HAR1 == 1 & HAR3 == 1 ~ "Current",
        HAR1 == 1 & HAR3 == 2 ~ "Former",
        TRUE ~ NA_character_
      ),
      levels = c("Never", "Former", "Current")
    )
  ) %>%  
  filter(ALC != "High alcohol consumption") %>%
  filter(HEP_MED != "Yes") %>%
  filter(!(
    GUPHSPFR == 2 &
      LIVER_ENZ == "Raised" &
      (
        HBP == 1 |
          HCP == 1 |
          IRON_OVERLOAD == "Iron overload"
      )
  )) %>%
  mutate(
    LIVER_GROUP = factor(
      case_when(
        GUPHSPFR == 2 & LIVER_ENZ == "Normal" ~ "NAFLD",
        
        GUPHSPFR == 2 &
          LIVER_ENZ == "Raised" &
          HBP == 2 &
          HCP == 2 &
          IRON_OVERLOAD == "No iron overload" ~ "NASH",
        
        GUPHSPFR == 1 ~ "Neither",
        
        TRUE ~ NA_character_
      ),
      levels = c("Neither", "NAFLD", "NASH")
    ) 
  )%>%
  mutate(
    fasting_status = factor(
      case_when(
        PHPFAST >= 8 ~ "Fasting",
        PHPFAST < 8 ~ "Not fasting",
        TRUE ~ NA_character_
      ),
      levels = c("Fasting", "Not fasting")
    ),
    
    diabetes = factor(
      case_when(
        # Keep self-reported diabetes because the paper includes doctor diagnosis
        HAD1 == 1 ~ "Yes",
        
        fasting_status == "Fasting" & SGPSI >= 7.0 ~ "Yes",
        
        fasting_status == "Not fasting" & SGPSI >= 11.0 ~ "Yes",
        
        HAD1 == 2 &
          (
            (fasting_status == "Fasting" & SGPSI < 7.0) |
              (fasting_status == "Not fasting" & SGPSI < 11.0)
          ) ~ "No",
        
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  ) %>%
  mutate(
    hypercholesterolaemia = factor(
      case_when(
        HAE7 == 1 ~ "Yes",
        !is.na(TCPSI) & TCPSI > 6.2 ~ "Yes",
        
        HAE7 == 2 ~ "No",
        !is.na(TCPSI) & TCPSI <= 6.2 ~ "No",
        
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  )

# Data subset including High Alcohol Consumption ----------
nhanes_work2 <- nhanes_merged %>%
  mutate(
    ALC = factor(
      case_when(
        MAPE1 == 2 ~ "No alcohol consumption",
        HSSEX == 2 & DRINK >= 1 ~ "High alcohol consumption",
        HSSEX == 1 & DRINK >= 2 ~ "High alcohol consumption",
        TRUE ~ "Low-Moderate alcohol consumption"
      ),
      levels = c(
        "No alcohol consumption",
        "Low-Moderate alcohol consumption",
        "High alcohol consumption"
      )
    )
  ) %>%
  mutate(
    LIVER_ENZ = factor(
      case_when(
        HSSEX == 1 & (ATPSI > 40 | ASPSI > 37) ~ "Raised",
        HSSEX == 2 & (ATPSI > 31 | ASPSI > 31) ~ "Raised",
        HSSEX == 1 & ATPSI <= 40 & ASPSI <= 37 ~ "Normal",
        HSSEX == 2 & ATPSI <= 31 & ASPSI <= 31 ~ "Normal",
        TRUE ~ NA_character_
      ),
      levels = c("Normal", "Raised")
    )
  ) %>%
  mutate(
    IRON_OVERLOAD = factor(
      case_when(
        HSSEX == 1 & TFR > 45 & FRPSI > 500 ~ "Iron overload",
        HSSEX == 2 & TFR > 45 & FRPSI > 400 ~ "Iron overload",
        !is.na(TFR) & !is.na(FRPSI) ~ "No iron overload",
        TRUE ~ NA_character_
      ),
      levels = c("No iron overload", "Iron overload")
    )
  )%>%
  mutate(
    SED = factor(
      if_else(
        if_all(c(HAT2, HAT4, HAT6, HAT8, HAT10, HAT12, HAT14, HAT16, HAT18), ~ .x == 2),
        "Sedentary",
        "Not sedentary"
      ),
      levels = c("Not sedentary", "Sedentary")
    )) %>% 
  mutate(
    SMOKE = factor(
      case_when(
        HAR1 == 2 ~ "Never",
        HAR1 == 1 & HAR3 == 1 ~ "Current",
        HAR1 == 1 & HAR3 == 2 ~ "Former",
        TRUE ~ NA_character_
      ),
      levels = c("Never", "Former", "Current")
    )
  ) %>%  
  filter(HEP_MED != "Yes") %>%
  filter(!(
    GUPHSPFR == 2 &
      LIVER_ENZ == "Raised" &
      (
        HBP == 1 |
          HCP == 1 |
          IRON_OVERLOAD == "Iron overload"
      )
  )) %>%
  mutate(
    LIVER_GROUP = factor(
      case_when(
        GUPHSPFR == 2 & LIVER_ENZ == "Normal" ~ "NAFLD",
        
        GUPHSPFR == 2 &
          LIVER_ENZ == "Raised" &
          HBP == 2 &
          HCP == 2 &
          IRON_OVERLOAD == "No iron overload" ~ "NASH",
        
        GUPHSPFR == 1 ~ "Neither",
        
        TRUE ~ NA_character_
      ),
      levels = c("Neither", "NAFLD", "NASH")
    ) 
  )%>%
  mutate(
    fasting_status = factor(
      case_when(
        PHPFAST >= 8 ~ "Fasting",
        PHPFAST < 8 ~ "Not fasting",
        TRUE ~ NA_character_
      ),
      levels = c("Fasting", "Not fasting")
    ),
    
    diabetes = factor(
      case_when(
        # Keep self-reported diabetes because the paper includes doctor diagnosis
        HAD1 == 1 ~ "Yes",
        
        fasting_status == "Fasting" & SGPSI >= 7.0 ~ "Yes",
        
        fasting_status == "Not fasting" & SGPSI >= 11.0 ~ "Yes",
        
        HAD1 == 2 &
          (
            (fasting_status == "Fasting" & SGPSI < 7.0) |
              (fasting_status == "Not fasting" & SGPSI < 11.0)
          ) ~ "No",
        
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  ) %>%
  mutate(
    hypercholesterolaemia = factor(
      case_when(
        HAE7 == 1 ~ "Yes",
        !is.na(TCPSI) & TCPSI > 6.2 ~ "Yes",
        
        HAE7 == 2 ~ "No",
        !is.na(TCPSI) & TCPSI <= 6.2 ~ "No",
        
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  )



# Exporting working data-------------------------------------------------------

write_csv(nhanes_work, "clean_data/nhanes_work.csv")
write_csv(nhanes_work2, "clean_data/nhanes_work_high_alc.csv")


# Create MASLD criteria & new liver groups ----------------------------------------------------

nhanes_masld <- nhanes_work %>%
  mutate(
    BMPWAIST = na_if(BMPWAIST, 88888),
    SGPSI = na_if(SGPSI, 88888),
    GHP = na_if(GHP, 8888),
    PEPMNK1R = na_if(PEPMNK1R, 888),
    PEPMNK5R = na_if(PEPMNK5R, 888), 
    TGPSI = na_if(TGPSI, 88888),
    HDPSI = na_if(HDPSI, 8888),
    CEPSI = na_if(CEPSI, 888888),
    TCPSI = na_if(TCPSI, 88888),
    PLP = na_if(PLP, 88888),
    FIB4 = (HSAGEIR * ASPSI) / (PLP * sqrt(ATPSI)),
    ) %>%
  mutate(
    
    # Criterion 1:
    # BMI >= 25 OR waist > 94 cm for men or > 80 cm for women
    MASLD_CRIT_ADIPOSITY =
      BMPBMI >= 25 |
      (HSSEX == 1 & BMPWAIST > 94) |
      (HSSEX == 2 & BMPWAIST > 80),
    
    # Criterion 2:
    # Fasting glucose >= 5.6 mmol/L OR HbA1c >= 5.7%
    # OR doctor-diagnosed diabetes
    MASLD_CRIT_GLUCOSE =
      (fasting_status == "Fasting" & SGPSI >= 5.6) |
      GHP >= 5.7 |
      HAD1 == 1,
    
    # Criterion 3:
    # Systolic BP >= 130 OR diastolic BP >= 85 mmHg
    MASLD_CRIT_BP =
      PEPMNK1R >= 130 |
      PEPMNK5R >= 85,
    
    # Criterion 4:
    # Triglycerides >= 1.70 mmol/L
    MASLD_CRIT_TRIGLYCERIDES =
      TGPSI >= 1.70,
    
    # Criterion 5:
    # HDL <= 1.0 mmol/L for men or <= 1.3 mmol/L for women
    MASLD_CRIT_LOW_HDL =
      (HSSEX == 1 & HDPSI <= 1.0) |
      (HSSEX == 2 & HDPSI <= 1.3),
    
    # Count how many of the five criteria each participant meets
    MASLD_CRITERIA_COUNT = rowSums(
      across(
        c(
          MASLD_CRIT_ADIPOSITY,
          MASLD_CRIT_GLUCOSE,
          MASLD_CRIT_BP,
          MASLD_CRIT_TRIGLYCERIDES,
          MASLD_CRIT_LOW_HDL
        )
      ),
      na.rm = TRUE
    ),
  
    LIVER_GROUP2 = factor(
      case_when(
        GUPHSPFR == 2 & LIVER_ENZ == "Normal" & MASLD_CRITERIA_COUNT >= 1 ~ "MASLD",
        
        GUPHSPFR == 2 &
          LIVER_ENZ == "Raised" &
          HBP == 2 &
          HCP == 2 &
          IRON_OVERLOAD == "No iron overload" ~ "MASH",
        
        GUPHSPFR == 1 ~ "Neither",
      
      TRUE ~ NA_character_
      ),
      levels = c("MASLD", "MASH", "Neither")
    )
  )

table(nhanes_masld$LIVER_GROUP2)
table(nhanes_masld$LIVER_GROUP)

# for ID metALD
#nhanes_masld <- nhanes_masld %>%
#  mutate(
#    DRINKS2 = DRINK*14
#  ) %>%
#  mutate(
#    ALC_INTAKE = factor(
#      case_when(
#      DRINKS2 >=50 & HSSEX == 1 ~ "MASLD Pre-dominant",
#      DRINKS2 >= 50 & HSSEX == 2 ~ "MASLD Pre-dominant",
#      DRINKS2 > 50 & HSSEX == 2 ~ "ALD Pre-dominant",
#      DRINKS2 > 60 & HSSEX == 1 ~ "ALD Pre-dominant",
#      TRUE ~ NA_character_
#    ),
#    levels = c("MASLD Pre-dominant", "ALD Pre-dominant")
#  )
#  )



# Summaries---------------------------------------------------------------------
#table(nhanes_masld$ALC_INTAKE)

table(nhanes_masld$MASLD)
boxplot(nhanes_masld$FIB4)
summary(nhanes_masld$FIB4)
table(nhanes_masld$LIVER_GROUP2)
table(nhanes_masld$LIVER_GROUP)
table(nhanes_masld$MASLD_CRITERIA_COUNT)


# Figures-----------------------------------------------------------------------
ggplot(
  data = nhanes_masld, aes(x = factor(MASLD_CRITERIA_COUNT))) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(x = "Number of MASLD Criteria", y = "Count") 

criterion_summary <- nhanes_masld %>%
  summarise(
    `Adiposity` = sum(MASLD_CRIT_ADIPOSITY, na.rm = TRUE),
    `Glucose` = sum(MASLD_CRIT_GLUCOSE, na.rm = TRUE),
    `Blood Pressure` = sum(MASLD_CRIT_BP, na.rm = TRUE),
    `Triglycerides` = sum(MASLD_CRIT_TRIGLYCERIDES, na.rm = TRUE),
    `Low HDL` = sum(MASLD_CRIT_LOW_HDL, na.rm = TRUE)
    ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Criteria", values_to = "Count"
  )

print(criterion_summary)

ggplot(data = criterion_summary, aes(x = Criteria, y = Count)) +
  geom_col() +
  geom_text(aes(label = Count), vjust = -0.5)



# Adding blood pressure treatment variable and statin medication (used later for PREVENT SCORE)

HBP_TREAT <- adult %>%
  select(SEQN, HAE5A) %>%
  mutate(
    HAE5A = na_if(HAE5A, 8)
  ) %>%
  mutate(HBP_MED = if_else(HAE5A == 1, TRUE, FALSE))

# check if those taking HBP_med have blood pressure measurements

nhanes_masld <- nhanes_masld %>%
  mutate(SMOKING = if_else(SMOKE=="Current", TRUE, FALSE)) %>%
  mutate(DIABETES = if_else(diabetes=="Yes", TRUE, FALSE)) %>%
  mutate(
    SEX = factor(
      case_when(
        HSSEX == 1 ~ "m",
        HSSEX == 2 ~ "f",
        TRUE ~ NA_character_
      ),
      levels = c("m", "f")
    ))


statin_med_codes <- c(1664, 3064, 3088)

meds_clean <- meds %>%
  select(SEQN, HQRXCODE) %>%
  group_by(SEQN) %>%
  summarise(
    STATIN_MED = factor(
      if_else(
        any(HQRXCODE %in% statin_med_codes, na.rm = TRUE),
        "TRUE",
        "FALSE"
      ),
      levels = c("TRUE", "FALSE")
    ),
    .groups = "drop"
  )
table(meds_clean$STATIN_MED) # 189 yes


nhanes_masld <- nhanes_masld %>%
  left_join(HBP_TREAT, by = "SEQN") %>%
  left_join(meds_clean, by = "SEQN")



# Export the modified dataset ----------------------------------------------
write_csv(nhanes_masld, "clean_data/nhanes_masld.csv")
