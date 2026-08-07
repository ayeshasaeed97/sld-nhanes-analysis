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

table(nhanes_work$LIVER_GROUP)
table(nhanes_work$ALC)
table(nhanes_work$IRON_OVERLOAD)
table(nhanes_work$LIVER_ENZ)
table(nhanes_work$SMOKE)

prop.table(table(nhanes_work$LIVER_GROUP, nhanes_work$ALC), margin = 1)
prop.table(table(nhanes_work$LIVER_GROUP, nhanes_work$SMOKE), margin = 1)



# Dietary inclusion -----

#library(remotes)
#library(dietaryindex)
#library(gtsummary)
#cff <- read.csv("/Users/ayeshasaeed/Documents/MSAB/Practical Training/NHANES/raw data/cff_output.csv")
#nutr <- read.csv('/Users/ayeshasaeed/Documents/MSAB/Practical Training/NHANES/raw data/examdr_output.csv')
#iff <- read.csv('/Users/ayeshasaeed/Documents/MSAB/Practical Training/NHANES/raw data/iff_output.csv')
#psdb <- read.csv('/Users/ayeshasaeed/Documents/MSAB/Practical Training/NHANES/raw data/pyramid_output.csv')



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
