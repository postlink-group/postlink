## code to prepare `lifem` dataset

library(tidyverse)
library(readxl)
library(piggyback)

# https://www.openicpsr.org/openicpsr/project/155186/version/V5/view
# Data from January 2022
master <- readRDS("lifem_master.rds")
death <- readRDS("lifem_death.rds")

# Link master and death data sets within life-m database
# using unique "lifemid" identifier
df <- left_join(master, death, by = "lifemid")

# focus on life-m individuals w/ yob in [1883, 1906]
# without missing year of death (only focus on complete cases)
df <- df %>% filter((yob >= 1883 & yob <= 1906) & !is.na(dob) & !is.na(yod))

# add age at death variable
df$surv_age <- df$yod  - df$yob

# remove individuals with age > 95 (such individuals should be
# top-coded per latest life-m version)
df <- df %>% filter(surv_age <= 95)

# hand-linked" variable (used as the latent match indicator in this application)
df$hndlnk[(df$lnkd %in% c(0,2)) & (df$lnkm %in% c(0,2)) &
             (df$lnkc40 %in% c(0,2)) & (df$lnkc20 %in% c(0,2)) &
                         (df$lnkc10 %in% c(0,2)) & (df$lnkc00 %in% c(0,2)) &
                                     (df$lnkc80 %in% c(0,2))] <- 0
 df$hndlnk <- as.factor(df$hndlnk)
 levels(df$hndlnk) <- c("Purely Machine-Linked", "Hand-Linked At Some Level")
# state-link variable - subset individuals linked to Ohio
df$stlnk <- as.factor(df$stlnk) # LIFE-M state
levels(df$stlnk) <- c("North Carolina", "Ohio")
df <- df %>% filter(stlnk == "Ohio") # 156,924 life-m individuals

# commonness score of first and last name (commf and comml, respectively)
# focus on individuals with both scores available
lifem <- df %>% select(lifemid, yob, yod, surv_age, hndlnk, commf, comml)
lifem <- lifem %>% filter(comml <= 1)
lifem$uyob <- (lifem$yob - min(lifem$yob)) / (max(lifem$yob) - min(lifem$yob))

# Save the full data set to add to current release
lifem_data <- lifem
lifem_data$hndlnk <- ifelse(lifem_data$hndlnk == "Hand-Linked At Some Level",
                            TRUE, FALSE)
saveRDS(lifem_data, file = "data-raw/lifem.rds")
pb_upload("data-raw/lifem.rds", repo = "postlink-group/postlink")

### Demo Data (n = 3,238) ###
### 2:1 (HL-ML Ratio)
hl <- lifem %>% filter(hndlnk == "Hand-Linked At Some Level") #2159 life-m individuals
ml <- lifem %>% filter(hndlnk == "Purely Machine-Linked") #154294 life-m individuals

n = nrow(hl)*0.5

set.seed(123)
ml <- sample_n(ml, n)
lifem <- rbind(hl, ml)

# Save the cleaned demo data to the package's data/ directory
lifem_data <- lifem
lifem_data$hndlnk <- ifelse(lifem_data$hndlnk == "Hand-Linked At Some Level",
                            TRUE, FALSE)
usethis::use_data(lifem, overwrite = TRUE)
