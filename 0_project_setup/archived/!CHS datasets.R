# CHS Datasets for CC and Longitudinal Analysis
library(tidyverse)
library(labelled)
library(DescTools)
library(fs)
library(here)
library(Hmisc)

source(here::here("!directories.R"))

# ADA diabetes diagnostic criteria (https://www.diabetes.org/a1c/diagnosis)
# Prediabetes: 
#1) fasting glucose between 100 mg/dL and 125 mg/dL
#2) 2-hour glucose between 140 and 199 mg/dL
#3) HbA1c between 5.7 and 6.4%

#type 2 diabetes
#1) fasting glucose >= 126 mg/dL
#2) 2-hour glucose >= 200 mg/dL
#3) HbA1c >= 6.5%

## MetaCHEM ##

### merge asa24 data and metachem clinic data ###
diet <- readRDS("/Volumes/projects_eh/MetaChem/MetaChem Follow up Visit/Datasets/diet.rds") #most recent diet data
load("/Volumes/projects_eh/MetaChem/MetaChem Follow up Visit/Data Quality Reports/out.Rdata") #most recent metachem data
#make correct CHS ID variable
#diet$id <- as.numeric(sapply(strsplit(diet$UserName, split='MEC', fixed=TRUE), function(x) (x[2]))) 
#diet <- diet %>%
  #rename(HEI = HEI2015_TOTAL_SCORE,
        # MDS = MDS_score,
         #DASH = DASH_Score_Mellan,
         #DII = DII_score_g,)
metachem <- merge(out, diet, by = "id", all=T) %>%
  remove_labels()
metachem$id <- as.numeric(metachem$id)
metachem <- subset(metachem, !is.na(HEI)) %>%
  mutate(bodyfatpercent = ifelse(!is.na(bodyfatpercent), bodyfatpercent, bodyfatpercent2),
         bmicat = ifelse(bmi <= 25.0, "Normal Weight", ifelse(bmi <= 30, "Overweight", "Obese")),
         fbg = ifelse(is.na(og_glu_5), og_glu_5_bedside, og_glu_5),
         glu2hr = ifelse(is.na(og_glu120), og_glu120_bedside, og_glu120),
         glucauc1=0.58*(og_glu_5+og_glu30)/2,
         glucauc2=0.5*(og_glu30+og_glu60)/2,
         glucauc3=0.5*(og_glu60+og_glu90)/2,
         glucauc4=0.5*(og_glu90+og_glu120)/2) %>%
  mutate(gluc_auc= glucauc1+glucauc2+glucauc3+glucauc4)

# # Select GLU to calculate AUC
glu <- metachem %>%
  #filter(is.na(og_glu120)) %>% 
  dplyr::select(id, og_glu_5, og_glu30, og_glu60,og_glu90, og_glu120) %>%
  pivot_longer(., cols = og_glu_5:og_glu120, names_to = "time", values_to = "glu")
  
glu <- glu %>% 
  mutate(time = case_when(time == "og_glu_5" ~ 0,
                          time == "og_glu30" ~ 30,
                          time == "og_glu60" ~ 60,
                          time == "og_glu90" ~ 90,
                          time == "og_glu120" ~ 120))

isi <- glu %>%
  group_by(id) %>%
  arrange(time) %>%
  summarise(g0 = glu[1],
            gauc = AUC(time/60, glu),
            glu_inc_auc = gauc - g0,
            max_glu = max(glu)) %>%
  ungroup()

metachem <- full_join(metachem, isi) 

label(metachem$fbg) <- "Fasting Glucose"
label(metachem$glu2hr) <- "Two Hour Glucose"
label(metachem$gauc) <- "Glucose AUC"
label(metachem$bmicat) <- "BMI Category"
label(metachem$bodyfatpercent) <- "Body Fat Percent"
label(metachem$diabetes) <- "Diabetes"
label(metachem$cage) <- "Age"
label(metachem$male_factor) <- "Sex"
label(metachem$race) <- "Race"

metachem$id <- remove_labels(metachem$id) #final metachem dataset

#clean up
rm(diet, out, isi, glu)

#subset to variables of interest for longitudinal analyses
metachem2 <- dplyr::select(metachem, id, male_factor, cage, race2, fbg, glu2hr, gauc, hba1c, bodyfatpercent, bmi,
                           HEI, MDS, DASH, DII, KCAL_mean, fatmass, fatfreemass, ffmi, fatmass_height_ratio,
                           android_gynoid_ratio, trunk_leg_ratio, trunk_limb_ratio, dexa_vat_vol, 
                           diabetes, diabetesbinary)

## MetaAIR ##
# (dataset from Will)
metaair <- read.csv(fs::path(dir_data_local, "Fullindices.csv")) %>%
  mutate(race2 = as.factor(race_eth) %>%
           fct_recode("White" = "0", "Hispanic" = "1", "Other" = "2"),
         diabetes = ifelse((!is.na(gluc_fasting) & gluc_fasting >= 126)|(!is.na(gluc_120min) & gluc_120min >= 200)|(!is.na(hba1c) & hba1c >= 6.5), "T2D",
                           ifelse((!is.na(gluc_fasting) & gluc_fasting >= 100)|(!is.na(gluc_120min) & gluc_120min >= 140)|(!is.na(hba1c) & hba1c >= 5.7), "prediabetes",
                                  ifelse(!is.na(gluc_fasting)|!is.na(gluc_120min)|!is.na(hba1c), "no diabetes", NA))),
         diabetesbinary = ifelse((!is.na(gluc_fasting) & gluc_fasting >= 100)|(!is.na(gluc_120min) & gluc_120min >= 140) | (!is.na(hba1c) & hba1c >= 5.7), "T2D or prediabetes",
                                 ifelse(!is.na(gluc_fasting)|!is.na(gluc_120min)|!is.na(hba1c), "no diabetes", NA)),
         sex = as.factor(male) %>%
           fct_recode("female" = "0", "male" = "1"))

temp <- read.csv(fs::path(dir_data_local, "DEXA2.csv")) %>%
  mutate(bmd = ifelse(is.na(bmd), bmd2, bmd),
         fatmass = ifelse(is.na(fatmass), fatmass2, fatmass),
         fatfreemass = ifelse(is.na(fatfreemass), fatfreemass2, fatfreemass),
         dexa_totalmass = ifelse(is.na(dexa_totalmass), dexa_totalmass2, dexa_totalmass),
         ffmi = (fatfreemass*0.001)/((heightcm*0.01)**2)) 

#calculate AUC
glu <- temp %>%
  #filter(is.na(og_glu120)) %>% 
  dplyr::select(id, gluc_fasting, gluc_30min, gluc_60min, gluc_90min, gluc_120min) %>%
  pivot_longer(., cols = gluc_fasting:gluc_120min, names_to = "time", values_to = "glu")

glu <- glu %>% 
  mutate(time = case_when(time == "gluc_fasting" ~ 0,
                          time == "gluc_30min" ~ 30,
                          time == "gluc_60min" ~ 60,
                          time == "gluc_90min" ~ 90,
                          time == "gluc_120min" ~ 120))

isi <- glu %>%
  group_by(id) %>%
  arrange(time) %>%
  summarise(g0 = glu[1],
            gauc = AUC(time/60, glu),
            glu_inc_auc = gauc - g0,
            max_glu = max(glu)) %>%
  ungroup()

metaair <- full_join(metaair, isi) 

metaair <- merge(metaair, temp)
rm(temp, glu, isi)

#subset to variables of interest for longitudinal analyses
#colnames(metaair)
metaair2 <- dplyr::select(metaair, id, sex, cage, race2, educ, DASH, HEI, DII, MDS, kcal, gluc_fasting, gluc_120min, 
                          hba1c, bodyfatpercent_convert, bmi_avg, fatmass, fatfreemass, ffmi, fatmass_height_ratio,
                          android_gynoid_ratio, trunk_leg_ratio, trunk_limb_ratio, dexa_vat_vol, gauc, 
                          diabetes, diabetesbinary)

#subset metaair to just those ids with follow up
metaair2 <- metaair2 %>%
  filter(id %in% metachem2$id)

metachem <- full_join(metachem, metaair2[c(1,5)]) 

#make variable names the same
metachem2 <- metachem2 %>%
  rename(sex = male_factor,
         kcal = KCAL_mean)
metaair2 <- metaair2 %>%
  rename(bodyfatpercent = bodyfatpercent_convert,
         bmi = bmi_avg,
         fbg = gluc_fasting,
         glu2hr = gluc_120min)

#merge and create visit #
metaair2$visit <- 1
metachem2$visit <- 2

metachem2 <- remove_labels(metachem2)
metachem2$sex <- as.numeric(metachem2$sex)
metaair2$sex <- as.numeric(metaair2$sex)
metanutrition_long <- full_join(metaair2, metachem2)

#create change in diet index, change in outcome variables
#colnames(metanutrition_long) #check this if dataset has changed 
nutrition <- metanutrition_long[,c(1,3,6:24)] %>%
  group_by(id) %>%
  summarise(across(everything(), ~diff(.x, na.rm=T))) %>%
  ungroup()
temp <- metaair2 %>% filter(id %in% nutrition$id)
temp2 <- dplyr::select(metachem2, id, diabetes, diabetesbinary)
#colnames(temp)
colnames(temp)[c(3,6:26)] = paste("baseline", names(temp)[c(3,6:26)], sep="_")
colnames(nutrition)[c(2:21)] = paste("change", names(nutrition)[c(2:21)], sep = "_")
metanutrition_wide <- left_join(left_join(temp, nutrition), temp2)

rm(nutrition, temp, temp2, metaair2, metachem2)
saveRDS(metaair, "/Volumes/projects_eh/MetaChem/MetaChem Follow up Visit/Datasets/metaair.rds")
saveRDS(metanutrition_long, "/Volumes/projects_eh/MetaChem/MetaChem Follow up Visit/Datasets/metanutrition_long.rds")
saveRDS(metanutrition_wide, "/Volumes/projects_eh/MetaChem/MetaChem Follow up Visit/Datasets/metanutrition_wide.rds")
