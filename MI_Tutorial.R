
############ Produce MAR covariate from complete data ############
##################################################################

### Code copied from How to generate missing values?, Teresa Alves de Sousa, 
### Imke Mayer

## This part of the code is to explore what package does what

# Load the required packages

library(mice)
library(MASS) # for mvrnorm function to generate complete dataset
library(dplyr)
library(tidyr)
library(ggplot2)
library(devtools)
library(gdata)
# library(amputation.R) # supposedly used for produce_NA function
library(VIM) # for matrixplot function
 
# Load the function to generate missing data

source_url('https://raw.githubusercontent.com/R-miss-tastic/website/master/static/how-to/generate/amputation.R') 
# perhaps this works instead of library(amputation.R)

set.seed(123)

# Sample data generation ------------------------------------------------------
# Generate complete data

mu.X <- c(1, 1)
Sigma.X <- matrix(c(1, 1, 1, 4), nrow = 2)
n <- 100
X.complete.cont <- mvrnorm(n, mu.X, Sigma.X)

lambda <- 0.5
X.complete.discr <- rpois(n, lambda)

n.cat <- 5
X.complete.cat <- rbinom(n, size = 5, prob = 0.5)

X.complete <- data.frame(cbind(X.complete.cont, X.complete.discr, X.complete.cat))
X.complete[,4] <- as.factor(X.complete[,4])
levels(X.complete[,4]) <- c("F", "E", "D", "C", "B", "A")


### Minimal set of arguments for produce_NA function 

# data: the initial data (can be complete or incomplete) as a matrix or data.frame

# one of "MCAR", "MAR", "MNAR" (default: "MCAR")

# perc.missing: the proportion of new missing values among the initially observed values (default: 0.5)


### Value - produce_NA returns a list containing three elements:

# data.init: the initial data

# data.incomp: the data with the newly generated missing values (and the initial missing values if applicable)

# idx_newNA: a matrix indexing only the newly generated missing values


### Example

## On complete data

# Minimal example for generating missing data ------------------------

X.miss <- produce_NA(X.complete, mechanism="MCAR", perc.missing = 0.2)

X.mcar <- X.miss$data.incomp
R.mcar <- X.miss$idx_newNA

matrixplot(X.mcar, cex.axis = 0.5, interactive = F)

## On incomplete data

# Minimal example for generating missing data on an incomplete data set ------------------------

X.miss <- produce_NA(rbind(X.complete[1:50,], X.mcar[51:100,]), mechanism="MCAR", perc.missing = 0.2)

X.mcar <- X.miss$data.incomp
R.mcar <- X.miss$idx_newNA

writeLines(paste0("Percentage of newly generated missing values: ", 100*sum(R.mcar)/prod(dim(R.mcar)), " %"))

matrixplot(X.mcar, cex.axis = 0.5, interactive = F)

## This part is the useful code that I need for my tutorial

# Follows are some useful arguments in produce_NA function and what they do

# 1. idx.incomplete: provides the position of the subset of variables where 
# missing values are to be generated.

# 2. idx.covariates: specifies which variables will be used in the missingness 
# model for MAR or MNAR mechanisms.

# 3. weights.covariates: specifies the corresponding weights of the variables  
# that will be used in the missingness model

# Note: 1 and 2 are complimentary. If idx.incomplete = c(1,1,0,0,...,0), then 
# idx.covariates must be of the form c(0,0,*,*,...,*)

# 4. X.miss$data.incomp to get the incomplete data

# 5. logit.model - specify which of the four logistic distribution functions 
# is to be implemented. Default is right-tailed ("RIGHT"), at which cases 
# with high weighted sum scores will have a higher probability to have missing 
# values, compared to cases with low weighted sum scores


### Now is the real data analysis task - creating missingness for the 
# warfarin dataset

# Read the warfarin dataset

setwd("C:/Users/u0164053/OneDrive - KU Leuven/PhD in KU Leuven/Projects/Multiple Imputation Tutorial in Pharmacometrics/Case study")
warfarin <- read.csv("warfarin_data.csv")
View(warfarin)

# Filter out dvid = 1 - focus on the PK part only. Also keep rows with dvid = "."

warfarin_pk <- warfarin %>% filter(dvid == 1 | dvid == ".")

## Make the dataset conform to NONMEM format

# First, capitalize the column names

colnames(warfarin_pk) <- toupper(colnames(warfarin_pk))

# Then, add the EVID column, which equals 0 if AMT equals ".", and 1 otherwise

warfarin_pk$EVID <- ifelse(warfarin_pk$AMT == ".", 0, 1)

# Next, add CMT column, which equals 1 if EVID equals 1, and 2 otherwise

warfarin_pk$CMT <- ifelse(warfarin_pk$EVID == 1, 1, 2)

# Then, add RATE column, which equals 0

warfarin_pk$RATE <- 0

# After that, add MDV, which equals EVID

warfarin_pk$MDV <- warfarin_pk$EVID

# Next, remove the DVID column

warfarin_pk <- warfarin_pk %>% select(-DVID)

# Finally, reorder the columns as follows: ID, TIME, DV, MDV, AMT, RATE, EVID, CMT, 
# WT, SEX, AGE

warfarin_pk <- warfarin_pk %>% select(ID, TIME, DV, MDV, AMT, RATE, EVID, CMT, 
                                      WT, SEX, AGE) 

## Then, perform data amputation

# Replace "." in DV and AMT with 0

warfarin_pk$DV <- ifelse(warfarin_pk$DV == ".", 0, warfarin_pk$DV)

warfarin_pk$AMT <- ifelse(warfarin_pk$AMT == ".", 0, warfarin_pk$AMT)

# Create a dataset with missing WT, as WT is the only covariate in the final 
# popPK model. Mechanism = "MAR". Covariates include: ID, TIME, DV, AMT, SEX, AGE 

# As WT is a baseline covariate, so only consider baseline rows in this procedure

warfarin_pk_baseline <- warfarin_pk %>% filter(TIME == 0)

# Change SEX into factor and AMT into numeric

warfarin_pk_baseline$SEX <- as.factor(warfarin_pk_baseline$SEX)

warfarin_pk_baseline$AMT <- as.numeric(warfarin_pk_baseline$AMT)

# Then, perform amputation with AMT, SEX, AGE as covariates, SEX & AGE have higher weights

warfarin_miss <- produce_NA(warfarin_pk_baseline, mechanism="MAR", perc.missing = 0.54, 
                            idx.incomplete = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0), 
                            idx.covariates = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1), 
                            weights.covariates = c(0, 0, 0, 0, 1/4, 0, 0, 0, 0, 3/8, 3/8)
                            ,logit.model = "RIGHT",seed=123) # Here I give SEX 
                            # and AGE a higher weight in the missingness model

# 250224 - I decided to remove AMT from the missing data model since it doesn't make sense
# in terms of MAR mechanism. Plus, AMT and WT in this particular dataset are highly correlated

warfarin_miss <- produce_NA(warfarin_pk_baseline, mechanism="MAR", perc.missing = 0.58, 
                            idx.incomplete = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0), 
                            idx.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1), 
                            weights.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1/2, 1/2)
                            ,logit.model = "RIGHT",seed=123)

warfarin_incomp <- warfarin_miss$data.incomp # extract the incomplete data

# count the number of missing WT cells 

sum(is.na(warfarin_incomp$WT)) # equals 15 - approximately 50% of the baseline
# exactly 46.9%


# 040324 - Generate datasets with 20, 50 and 80% WT missing -------------------- 

# 20% missing WT

warfarin_miss_20 <- produce_NA(warfarin_pk_baseline, mechanism="MAR", perc.missing = 0.20, 
                               idx.incomplete = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0), 
                               idx.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1), 
                               weights.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1/2, 1/2)
                               ,logit.model = "RIGHT",seed=123)

warfarin_incomp_20 <- warfarin_miss_20$data.incomp # extract the incomplete data

sum(is.na(warfarin_incomp_20$WT)) # equals 6 - 18.8% WT missing

# 30% missing WT

warfarin_miss_30 <- produce_NA(warfarin_pk_baseline, mechanism="MAR", perc.missing = 0.45, 
                               idx.incomplete = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0), 
                               idx.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1), 
                               weights.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1/2, 1/2)
                               ,logit.model = "RIGHT",seed=123)

warfarin_incomp_30 <- warfarin_miss_30$data.incomp # extract the incomplete data

sum(is.na(warfarin_incomp_30$WT)) # equals 10 - 31.2% WT missing

# 50% missing WT

warfarin_miss_50 <- produce_NA(warfarin_pk_baseline, mechanism="MAR", perc.missing = 0.583, 
                               idx.incomplete = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0), 
                               idx.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1), 
                               weights.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1/2, 1/2)
                               ,logit.model = "RIGHT",seed=123)

warfarin_incomp_50 <- warfarin_miss_50$data.incomp # extract the incomplete data

sum(is.na(warfarin_incomp_50$WT)) # equals 16 - exactly 50% WT missing

# 70% missing WT

warfarin_miss_70 <- produce_NA(warfarin_pk_baseline, mechanism="MAR", perc.missing = 0.6701, 
                               idx.incomplete = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0), 
                               idx.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1), 
                               weights.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1/3, 2/3)
                               ,logit.model = "RIGHT",seed=123)

warfarin_incomp_70 <- warfarin_miss_70$data.incomp # extract the incomplete data

sum(is.na(warfarin_incomp_70$WT)) # equals 16 - exactly 50% WT missing

# I will try 80% anyway, as 70% does not work!

# 80% missing WT

warfarin_miss_80 <- produce_NA(warfarin_pk_baseline, mechanism="MAR", perc.missing = 0.73, 
                               idx.incomplete = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0), 
                               idx.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1), 
                               weights.covariates = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1/2, 1/2)
                               ,logit.model = "RIGHT",seed=123)

warfarin_incomp_80 <- warfarin_miss_80$data.incomp # extract the incomplete data

sum(is.na(warfarin_incomp_80$WT)) # equals 25 - 78.1% WT missing

### Now perform data analysis --------------------------------------------------

## First, with the complete data -----------------------------------------------

# Change 0 in DV and AMT back to "."

warfarin_pk$DV <- ifelse(warfarin_pk$DV == 0, ".", warfarin_pk$DV)

warfarin_pk$AMT <- ifelse(warfarin_pk$AMT == 0, ".", warfarin_pk$AMT)

# Then, export it with the name warfarin_auth.csv

setwd("C:/Users/u0164053/OneDrive - KU Leuven/PhD in KU Leuven/Projects/Multiple Imputation Tutorial in Pharmacometrics/Case study/NONMEM datasets")
write.csv(warfarin_pk, "warfarin_auth.csv",quote = F, row.names = F)

# Check the median body weight

median(warfarin_pk_baseline$WT) # equals 71.7

summary(warfarin_pk_baseline$WT) # check the range as well

# 010324 - also check the mean (sd) of body weight

shapiro.test(warfarin_pk_baseline$WT) # p-value = 0.97891, normal

mean(warfarin_pk_baseline$WT) # equals 72.3

sd(warfarin_pk_baseline$WT) # equals 12.7


## Second, with the complete case analysis -------------------------------------


# Remove from warfarin_pk rows that has ID with missing WT in warfarin_incomp

warfarin_pk_cc <- warfarin_pk %>% filter(!(ID %in% warfarin_incomp$ID[is.na(warfarin_incomp$WT)]))

# 040324 - I do it for the 3 scenarios, 20, 50 and 80% missing WT

warfarin_pk_cc_20 <- warfarin_pk %>% filter(!(ID %in% warfarin_incomp_20$ID[is.na(warfarin_incomp_20$WT)]))

warfarin_pk_cc_30 <- warfarin_pk %>% filter(!(ID %in% warfarin_incomp_30$ID[is.na(warfarin_incomp_30$WT)]))

warfarin_pk_cc_50 <- warfarin_pk %>% filter(!(ID %in% warfarin_incomp_50$ID[is.na(warfarin_incomp_50$WT)]))

warfarin_pk_cc_80 <- warfarin_pk %>% filter(!(ID %in% warfarin_incomp_80$ID[is.na(warfarin_incomp_80$WT)]))


# Export it with the name warfarin_cc.csv

# So here I export as the beginning

write.csv(warfarin_pk_cc, "warfarin_cc.csv",quote = F, row.names = F)

# 040324 - Export different scenarios

write.csv(warfarin_pk_cc_20, "warfarin_cc_20.csv",quote = F, row.names = F)

write.csv(warfarin_pk_cc_30, "warfarin_cc_30.csv",quote = F, row.names = F)

write.csv(warfarin_pk_cc_50, "warfarin_cc_50.csv",quote = F, row.names = F)

write.csv(warfarin_pk_cc_80, "warfarin_cc_80.csv",quote = F, row.names = F)


## Third, median imputation ---------------------------------------------------

# Replace missing WT with the median of observed WT. But first replace WT in warfarin_pk dataset with WT from warfarin_incomp by ID

warfarin_pk$WT <- ifelse(warfarin_pk$ID %in% warfarin_incomp$ID, 
                         warfarin_incomp$WT[match(warfarin_pk$ID, warfarin_incomp$ID)], 
                         warfarin_pk$WT)

warfarin_pk$WT <- ifelse(is.na(warfarin_pk$WT), median(warfarin_pk_cc$WT), warfarin_pk$WT)

# Export it with the name warfarin_med.csv

write.csv(warfarin_pk, "warfarin_med.csv",quote = F, row.names = F)


# 040324 - Replace missing WT with the median of observed WT for the 3 scenarios and export them


# 20% missing WT

warfarin_pk_20 <- warfarin_pk

warfarin_pk_20$WT <- ifelse(warfarin_pk_20$ID %in% warfarin_incomp_20$ID, 
                         warfarin_incomp_20$WT[match(warfarin_pk_20$ID, warfarin_incomp_20$ID)], 
                         warfarin_pk_20$WT)

warfarin_pk_20$WT <- ifelse(is.na(warfarin_pk_20$WT), median(warfarin_pk_cc_20$WT), warfarin_pk_20$WT)

write.csv(warfarin_pk_20, "warfarin_med_20.csv",quote = F, row.names = F)

# 30% missing WT

warfarin_pk_30 <- warfarin_pk

warfarin_pk_30$WT <- ifelse(warfarin_pk_30$ID %in% warfarin_incomp_30$ID, 
                            warfarin_incomp_30$WT[match(warfarin_pk_30$ID, warfarin_incomp_30$ID)], 
                            warfarin_pk_30$WT)

warfarin_pk_30$WT <- ifelse(is.na(warfarin_pk_30$WT), median(warfarin_pk_cc_30$WT), warfarin_pk_30$WT)

write.csv(warfarin_pk_30, "warfarin_med_30.csv",quote = F, row.names = F)

# 50% missing WT

warfarin_pk_50 <- warfarin_pk

warfarin_pk_50$WT <- ifelse(warfarin_pk_50$ID %in% warfarin_incomp_50$ID, 
                         warfarin_incomp_50$WT[match(warfarin_pk_50$ID, warfarin_incomp_50$ID)], 
                         warfarin_pk_50$WT)

warfarin_pk_50$WT <- ifelse(is.na(warfarin_pk_50$WT), median(warfarin_pk_cc_50$WT), warfarin_pk_50$WT)

write.csv(warfarin_pk_50, "warfarin_med_50.csv",quote = F, row.names = F)

# 80% missing WT

warfarin_pk_80 <- warfarin_pk

warfarin_pk_80$WT <- ifelse(warfarin_pk_80$ID %in% warfarin_incomp_80$ID, 
                            warfarin_incomp_80$WT[match(warfarin_pk_80$ID, warfarin_incomp_80$ID)], 
                            warfarin_pk_80$WT)

warfarin_pk_80$WT <- ifelse(is.na(warfarin_pk_80$WT), median(warfarin_pk_cc_80$WT), warfarin_pk_80$WT)

write.csv(warfarin_pk_80, "warfarin_med_80.csv",quote = F, row.names = F)


## Fourth, multiple imputation ------------------------------------------------

# Produce 50 imputation datasets from warfarin_incomp as the % of WT missingness
# is approximately 50%

# Select variables to impute WT, including AMT, SEX & AGE

warfarin_impute <- warfarin_incomp %>% select(WT, AMT, SEX , AGE)

# Perform 50 imputations

imp0 <- mice(warfarin_impute, maxit=0)
meth <- imp0$method
meth["WT"] <- "norm"
maxit <- 20
nimp <- 50
warfarin_imputed <- mice::mice(data=warfarin_impute,remove.collinear = FALSE, method=meth,
                                         maxit=maxit, m=nimp, printFlag = FALSE,
                                         seed=123)

# 250224 - Remove AMT from the imputation model

warfarin_impute <- warfarin_incomp %>% select(WT, SEX , AGE)

# Perform 50 imputations

imp0 <- mice(warfarin_impute, maxit=0)
meth <- imp0$method
meth["WT"] <- "norm"
maxit <- 20
nimp <- 50
warfarin_imputed <- mice::mice(data=warfarin_impute,method=meth,
                               maxit=maxit, m=nimp, printFlag = FALSE,
                               seed=123)

# Increase to 60 imputations to account for 20% non-convergence rate in NONMEM

imp0 <- mice(warfarin_impute, maxit=0)
meth <- imp0$method
meth["WT"] <- "norm"
maxit <- 20
nimp <- 60
warfarin_imputed <- mice::mice(data=warfarin_impute,method=meth,
                               maxit=maxit, m=nimp, printFlag = FALSE,
                               seed=123)

# add imputed WT to warfarin_incomp

for (i in 1:60) {
  warfarin_incomp[[paste0("WT", sprintf("%02d", i))]] <- (complete(warfarin_imputed, action = i))$WT
}

z# Round up WT01 to WT60 in warfarin_incomp to 1 significant digit

for (i in 1:60) {
  warfarin_incomp[[paste0("WT", sprintf("%02d", i))]] <- round(warfarin_incomp[[paste0("WT", sprintf("%02d", i))]], 1)
}

# Copy columns WT01 to WT60 from warfarin_incomp to warfarin_pk, matching by ID
# It needs to be broken down into two parts

WT01_WT30 <- warfarin_incomp[, c(1, 12:41)]

WT31_WT60 <- warfarin_incomp[, c(1, 42:71)]

warfarin_MI01 <- merge(warfarin_pk, WT01_WT30, by = "ID", all.x = TRUE)

warfarin_MI02 <- merge(warfarin_pk, WT31_WT60, by = "ID", all.x = TRUE)

# Export it with the name warfarin_MI01.csv and warfarin_MI02.csv

write.csv(warfarin_MI01, "warfarin_MI01.csv",quote = F, row.names = F)

write.csv(warfarin_MI02, "warfarin_MI02.csv",quote = F, row.names = F)


## Pooling parameters ---------------------------------------------------------

# Read in the dataset
setwd("C:/Users/u0164053/OneDrive - KU Leuven/PhD in KU Leuven/Projects/Multiple Imputation Tutorial in Pharmacometrics/Case study/MI")
Pooled_War<-read.csv("Pooled_War.csv",sep=",")

# Creating Var (Variance) column 
Pooled_War$Var<-(Pooled_War$RSE*Pooled_War$Estimates/100)^2

# Then, pooling WT_CL

WT_CL<-mean(Pooled_War$Estimates[Pooled_War$Parameters == "WT_CL"]) #0.4794

W_WT_CL<-mean(Pooled_War$Var[Pooled_War$Parameters == "WT_CL"]) #0.0255248 #within imputation variance

B_WT_CL<-var(Pooled_War$Estimates[Pooled_War$Parameters == "WT_CL"]) #0.008850446 #between imputation variance

WT_CL
W_WT_CL
B_WT_CL


# Then, pooling WT_V

WT_V<-mean(Pooled_War$Estimates[Pooled_War$Parameters == "WT_V"]) #0.4794

W_WT_V<-mean(Pooled_War$Var[Pooled_War$Parameters == "WT_V"]) #0.0255248 #within imputation variance

B_WT_V<-var(Pooled_War$Estimates[Pooled_War$Parameters == "WT_V"]) #0.008850446 #between imputation variance

WT_V
W_WT_V
B_WT_V

## Note: Both WT_CL and WT_V are biased, probably because of the suboptimal imputation model.


## Finetuning the imputation model 280224 --------------------------------------

# MI scenario 03, 50% missing WT 060324 ----------------------------------------

# Variables for the imputation model

warfarin_impute_50 <- warfarin_pk %>% select(ID, TIME, DV, AMT, SEX, AGE, WT)

# Change AMT and DV to numeric

warfarin_impute_50$AMT <- as.numeric(warfarin_impute_50$AMT)

warfarin_impute_50$DV <- as.numeric(warfarin_impute_50$DV)

# Change SEX into factor

warfarin_impute_50$SEX <- as.numeric(warfarin_impute_50$SEX)

# Replace WT in warfarin_impute_50 dataset with WT from warfarin_incomp_50 by ID

warfarin_impute_50$WT <- ifelse(warfarin_impute_50$ID %in% warfarin_incomp_50$ID, 
                             warfarin_incomp_50$WT[match(warfarin_impute_50$ID, warfarin_incomp_50$ID)], 
                             warfarin_impute_50$WT)

# Leave out ID from warfarin_impute_50

warfarin_impute_50 <- warfarin_impute_50 %>% select(-ID)

# Perform 50 imputations

imp0 <- mice(warfarin_impute_50, maxit=0)
meth <- imp0$method
meth["WT"] <- "norm"
maxit <- 20
nimp <- 50
warfarin_imputed_50 <- mice::mice(data=warfarin_impute_50,method=meth,
                               maxit=maxit, m=nimp, printFlag = FALSE,
                               seed=123)

# add imputed WT to warfarin_pk_50

for (i in 1:50) {
  warfarin_pk_50[[paste0("WT", sprintf("%02d", i))]] <- (complete(warfarin_imputed_50, action = i))$WT
}

# Average WT01 to WT50 in warfarin_pk_50 by ID 

WT01_50 <- warfarin_pk_50 %>%
  group_by(ID) %>%
  summarise(across(WT01:WT50, mean, na.rm = TRUE))

View(WT01_50)

# Removing columns WT01 to WT50 from warfarin_pk_50

warfarin_pk_50 <- warfarin_pk_50 %>% select(-WT01:-WT50)

# Merge WT01_50 with warfarin_pk_50 by ID

warfarin_pk_50 <- merge(warfarin_pk_50, WT01_50, by = "ID", all.x = TRUE)

# Repeat the above steps from the rounding up one

# Round up WT01 to WT50 in warfarin_pk_50 to 1 significant digit

for (i in 1:50) {
  warfarin_pk_50[[paste0("WT", sprintf("%02d", i))]] <- round(warfarin_pk_50[[paste0("WT", sprintf("%02d", i))]], 1)
}

# Removing WT26 to WT50 columns in warfarin_pk_50 dataframe to make warfarin_MI_50_01 dataframe

warfarin_MI_50_01 <- warfarin_pk_50 %>% select(-WT26:-WT50)

# Similarly, warfarin_MI_50_02 is produced by removing WT01 to WT25 columns in warfarin_pk_50 dataframe

warfarin_MI_50_02 <- warfarin_pk_50 %>% select(-WT01:-WT25)

# Export it with the name warfarin_MI_50_01.csv and warfarin_MI_50_02.csv

write.csv(warfarin_MI_50_01, "warfarin_MI_50_01.csv",quote = F, row.names = F)

write.csv(warfarin_MI_50_02, "warfarin_MI_50_02.csv",quote = F, row.names = F)


# MI scenario 01, 20% missing WT -----------------------------------------------

# Select variables for the imputation model

warfarin_impute_20 <- warfarin_pk %>% select(ID, TIME, DV, AMT, SEX, AGE, WT)

# Change AMT and DV to numeric

warfarin_impute_20$AMT <- as.numeric(warfarin_impute_20$AMT)

warfarin_impute_20$DV <- as.numeric(warfarin_impute_20$DV)

# Change SEX into factor

warfarin_impute_20$SEX <- as.numeric(warfarin_impute_20$SEX)

# Replace WT in warfarin_impute_20 dataset with WT from warfarin_incomp_20 by ID

warfarin_impute_20$WT <- ifelse(warfarin_impute_20$ID %in% warfarin_incomp_20$ID, 
                             warfarin_incomp_20$WT[match(warfarin_impute_20$ID, warfarin_incomp_20$ID)], 
                             warfarin_impute_20$WT)

# Leave out ID from warfarin_impute_20

warfarin_impute_20 <- warfarin_impute_20 %>% select(-ID)

# Perform 20 imputations

imp0 <- mice(warfarin_impute_20, maxit=0)
meth <- imp0$method
meth["WT"] <- "norm"
maxit <- 20
nimp <- 20
warfarin_imputed_20 <- mice::mice(data=warfarin_impute_20,method=meth,
                               maxit=maxit, m=nimp, printFlag = FALSE,
                               seed=123)

# add imputed WT to warfarin_pk_20

for (i in 1:20) {
  warfarin_pk_20[[paste0("WT", sprintf("%02d", i))]] <- (complete(warfarin_imputed_20, action = i))$WT
}

# Average WT01 to WT20 in warfarin_pk_20 by ID 

WT01_20 <- warfarin_pk_20 %>%
  group_by(ID) %>%
  summarise(across(WT01:WT20, mean, na.rm = TRUE))

View(WT01_20)

# Removing columns WT01 to WT20 from warfarin_pk_20

warfarin_pk_20 <- warfarin_pk_20 %>% select(-WT01:-WT20)

# Merge WT01_20 with warfarin_pk_20 by ID

warfarin_pk_20 <- merge(warfarin_pk_20, WT01_20, by = "ID", all.x = TRUE)

# Round up WT01 to WT20 in warfarin_pk_20 to 1 significant digit

for (i in 1:20) {
  warfarin_pk_20[[paste0("WT", sprintf("%02d", i))]] <- round(warfarin_pk_20[[paste0("WT", sprintf("%02d", i))]], 1)
}

# Replace NA in WT with "."

warfarin_pk_20[is.na(warfarin_pk_20)] <- "."

# Export it with the name warfarin_MI_20.csv

write.csv(warfarin_pk_20, "warfarin_MI_20.csv",quote = F, row.names = F)


# MI scenario 02, 30% missing WT -----------------------------------------------

# Select variables for the imputation model

warfarin_impute_30 <- warfarin_pk %>% select(ID, TIME, DV, AMT, SEX, AGE, WT)

# Change AMT and DV to numeric

warfarin_impute_30$AMT <- as.numeric(warfarin_impute_30$AMT)

warfarin_impute_30$DV <- as.numeric(warfarin_impute_30$DV)

# Change SEX into factor

warfarin_impute_30$SEX <- as.numeric(warfarin_impute_30$SEX)

# Replace WT in warfarin_impute_30 dataset with WT from warfarin_incomp_30 by ID

warfarin_impute_30$WT <- ifelse(warfarin_impute_30$ID %in% warfarin_incomp_30$ID, 
                                warfarin_incomp_30$WT[match(warfarin_impute_30$ID, warfarin_incomp_30$ID)], 
                                warfarin_impute_30$WT)

# Leave out ID from warfarin_impute_30

warfarin_impute_30 <- warfarin_impute_30 %>% select(-ID)

# Perform 30 imputations

imp0 <- mice(warfarin_impute_30, maxit=0)
meth <- imp0$method
meth["WT"] <- "norm"
maxit <- 20
nimp <- 30
warfarin_imputed_30 <- mice::mice(data=warfarin_impute_30,method=meth,
                                  maxit=maxit, m=nimp, printFlag = FALSE,
                                  seed=123)

# add imputed WT to warfarin_pk_30

for (i in 1:30) {
  warfarin_pk_30[[paste0("WT", sprintf("%02d", i))]] <- (complete(warfarin_imputed_30, action = i))$WT
}

# Average WT01 to WT30 in warfarin_pk_30 by ID 

WT01_30 <- warfarin_pk_30 %>%
  group_by(ID) %>%
  summarise(across(WT01:WT30, mean, na.rm = TRUE))

View(WT01_30)

# Removing columns WT01 to WT30 from warfarin_pk_30

warfarin_pk_30 <- warfarin_pk_30 %>% select(-WT01:-WT30)

# Merge WT01_30 with warfarin_pk_30 by ID

warfarin_pk_30 <- merge(warfarin_pk_30, WT01_30, by = "ID", all.x = TRUE)

# Round up WT01 to WT30 in warfarin_pk_30 to 1 significant digit

for (i in 1:30) {
  warfarin_pk_30[[paste0("WT", sprintf("%02d", i))]] <- round(warfarin_pk_30[[paste0("WT", sprintf("%02d", i))]], 1)
}

# Replace NA in WT with "."

warfarin_pk_30[is.na(warfarin_pk_30)] <- "."

# Export it with the name warfarin_MI_30.csv

write.csv(warfarin_pk_30, "warfarin_MI_30.csv",quote = F, row.names = F)


# MI scenario 04, 80% missing WT -----------------------------------------------

# Variables for the imputation model

warfarin_impute_80 <- warfarin_pk %>% select(ID, TIME, DV, AMT, SEX, AGE, WT)

# Change AMT and DV to numeric

warfarin_impute_80$AMT <- as.numeric(warfarin_impute_80$AMT)

warfarin_impute_80$DV <- as.numeric(warfarin_impute_80$DV)

# Change SEX into factor

warfarin_impute_80$SEX <- as.numeric(warfarin_impute_80$SEX)

# Replace WT in warfarin_impute_80 dataset with WT from warfarin_incomp_80 by ID

warfarin_impute_80$WT <- ifelse(warfarin_impute_80$ID %in% warfarin_incomp_80$ID, 
                                warfarin_incomp_80$WT[match(warfarin_impute_80$ID, warfarin_incomp_80$ID)], 
                                warfarin_impute_80$WT)

# Leave out ID from warfarin_impute_80

warfarin_impute_80 <- warfarin_impute_80 %>% select(-ID)

# Perform 80 imputations

imp0 <- mice(warfarin_impute_80, maxit=0)
meth <- imp0$method
meth["WT"] <- "norm"
maxit <- 20
nimp <- 80
warfarin_imputed_80 <- mice::mice(data=warfarin_impute_80,method=meth,
                                  maxit=maxit, m=nimp, printFlag = FALSE,
                                  seed=123)

# add imputed WT to warfarin_pk_80

for (i in 1:80) {
  warfarin_pk_80[[paste0("WT", sprintf("%02d", i))]] <- (complete(warfarin_imputed_80, action = i))$WT
}

# Average WT01 to WT80 in warfarin_pk_80 by ID 

WT01_80 <- warfarin_pk_80 %>%
  group_by(ID) %>%
  summarise(across(WT01:WT80, mean, na.rm = TRUE))

View(WT01_80)

# Removing columns WT01 to WT80 from warfarin_pk_80

warfarin_pk_80 <- warfarin_pk_80 %>% select(-WT01:-WT80)

# Merge WT01_80 with warfarin_pk_80 by ID

warfarin_pk_80 <- merge(warfarin_pk_80, WT01_80, by = "ID", all.x = TRUE)

# Repeat the above steps from the rounding up one

# Round up WT01 to WT80 in warfarin_pk_80 to 1 significant digit

for (i in 1:80) {
  warfarin_pk_80[[paste0("WT", sprintf("%02d", i))]] <- round(warfarin_pk_80[[paste0("WT", sprintf("%02d", i))]], 1)
}

# Removing WT26 to WT80 columns in warfarin_pk_80 dataframe to make warfarin_MI_80_01 dataframe

warfarin_MI_80_01 <- warfarin_pk_80 %>% select(-WT26:-WT80)

# Similarly, warfarin_MI_80_02 is produced by removing WT01 to WT25 and WT 51 to 80 columns in warfarin_pk_80 dataframe

warfarin_MI_80_02 <- warfarin_pk_80 %>% select(-WT01:-WT25, -WT51:-WT80)

# Lastly, warfarin_MI_80_03 is produced by removing WT01 to WT50 columns in warfarin_pk_80 dataframe

warfarin_MI_80_03 <- warfarin_pk_80 %>% select(-WT01:-WT50)

# Export it with the name warfarin_MI_80_01.csv, warfarin_MI_80_02.csv and warfarin_MI_80_03.csv

write.csv(warfarin_MI_80_01, "warfarin_MI_80_01.csv",quote = F, row.names = F)

write.csv(warfarin_MI_80_02, "warfarin_MI_80_02.csv",quote = F, row.names = F)

write.csv(warfarin_MI_80_03, "warfarin_MI_80_03.csv",quote = F, row.names = F)

# Pooling parameters to check -------------------------------------------------

## My very first 60 imputations ------------------------------------------------

# Read in the dataset
setwd("C:/Users/u0164053/OneDrive - KU Leuven/PhD in KU Leuven/Projects/Multiple Imputation Tutorial in Pharmacometrics/Case study/MI")
Pooled_War<-read.csv("Pooled_War_Revised.csv",sep=",")

# Creating Var (Variance) column 
Pooled_War$Var<-(Pooled_War$RSE*Pooled_War$Estimates/100)^2


# Then, pooling WT_CL 

WT_CL<-mean(Pooled_War$Estimates[Pooled_War$Parameters == "WT_CL"]) #0.4794

W_WT_CL<-mean(Pooled_War$Var[Pooled_War$Parameters == "WT_CL"]) #0.0255248 #within imputation variance

B_WT_CL<-var(Pooled_War$Estimates[Pooled_War$Parameters == "WT_CL"]) #0.008850446 #between imputation variance

WT_CL # 0.60028
W_WT_CL
B_WT_CL

# Total variance

m<-50
T_WT_CL<-W_WT_CL+(1+1/m)*B_WT_CL #0.07616805
T_WT_CL

# Confidence interval - T distribution

m<-50 #number of imputations
k<-12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect)
n<-32

lambda_WT_CL<-(B_WT_CL+B_WT_CL/m)/T_WT_CL #proportion of variation attributable to the missing data

r_WT_CL<-(lambda_WT_CL)/(1-lambda_WT_CL) #relative increase in variance due to nonresponse

v_old_WT_CL<-(m-1)*(1+1/r_WT_CL^2) #old degree of freedom

v_com=n-k #degree of freedom of parameter estimate (WT_CL) in the hypothetically complete data

v_obs_WT_CL=((v_com+1)/(v_com+3))*v_com*(1-lambda_WT_CL) #observed data degrees of freedom that accounts for the missing information

v_WT_CL<-(v_old_WT_CL*v_obs_WT_CL)/(v_old_WT_CL+v_obs_WT_CL)

t_crit_WT_CL <- qt(0.975, v_WT_CL) 

lower_bound_WT_CL<-WT_CL-t_crit_WT_CL*sqrt(T_WT_CL) #1.425197
upper_bound_WT_CL<-WT_CL+t_crit_WT_CL*sqrt(T_WT_CL) #1.958517

# Here I approximate it with normal distribution to calculate RSE

sqrt(T_WT_CL)/WT_CL*100 #46.0%

lower_bound_WT_CL
upper_bound_WT_CL



# Then, pooling WT_V 

WT_V<-mean(Pooled_War$Estimates[Pooled_War$Parameters == "WT_V"]) #0.4794

W_WT_V<-mean(Pooled_War$Var[Pooled_War$Parameters == "WT_V"]) #0.0255248 #within imputation variance

B_WT_V<-var(Pooled_War$Estimates[Pooled_War$Parameters == "WT_V"]) #0.008850446 #between imputation variance

WT_V # 0.82766
W_WT_V
B_WT_V

# Total variance

m<-50
T_WT_V<-W_WT_V+(1+1/m)*B_WT_V #0.07616805
T_WT_V

# Confidence interval - T distribution

m<-50 #number of imputations
k<-12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect)
n<-32

lambda_WT_V<-(B_WT_V+B_WT_V/m)/T_WT_V #proportion of variation attributable to the missing data

r_WT_V<-(lambda_WT_V)/(1-lambda_WT_V) #relative increase in variance due to nonresponse

v_old_WT_V<-(m-1)*(1+1/r_WT_V^2) #old degree of freedom

v_com=n-k #degree of freedom of parameter estimate (WT_V) in the hypothetically complete data

v_obs_WT_V=((v_com+1)/(v_com+3))*v_com*(1-lambda_WT_V) #observed data degrees of freedom that accounts for the missing information

v_WT_V<-(v_old_WT_V*v_obs_WT_V)/(v_old_WT_V+v_obs_WT_V)

t_crit_WT_V <- qt(0.975, v_WT_V) 

lower_bound_WT_V<-WT_V-t_crit_WT_V*sqrt(T_WT_V) 
upper_bound_WT_V<-WT_V+t_crit_WT_V*sqrt(T_WT_V) 

# Here I approximate it with normal distribution to calculate RSE

sqrt(T_WT_V)/WT_V*100 # 24% 

lower_bound_WT_V
upper_bound_WT_V


## For the 30% missing WT scenario ---------------------------------------------

# Read in the dataset
setwd("C:/Users/u0164053/OneDrive - KU Leuven/PhD in KU Leuven/Projects/Multiple Imputation Tutorial in Pharmacometrics/Case study/MI pooling/30% WT missing")
Pooled_War_30<-read.csv("Pooled_War_30.csv",sep=",")

# Creating Var (Variance) column 
Pooled_War_30$Var<-(Pooled_War_30$RSE*Pooled_War_30$Estimates/100)^2


# Then, pooling WT_CL 

WT_CL<-mean(Pooled_War_30$Estimates[Pooled_War_30$Parameters == "WT_CL"]) #0.6358

W_WT_CL<-mean(Pooled_War_30$Var[Pooled_War_30$Parameters == "WT_CL"]) #0.08308769 #within imputation variance

B_WT_CL<-var(Pooled_War_30$Estimates[Pooled_War_30$Parameters == "WT_CL"]) #0.004717062 #between imputation variance

WT_CL # 0.60028
W_WT_CL
B_WT_CL

# Total variance

m<-30
T_WT_CL<-W_WT_CL+(1+1/m)*B_WT_CL #0.08796199
T_WT_CL

# Confidence interval - T distribution

m<-30 #number of imputations
k<-12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect)
n<-32

lambda_WT_CL<-(B_WT_CL+B_WT_CL/m)/T_WT_CL #proportion of variation attributable to the missing data

r_WT_CL<-(lambda_WT_CL)/(1-lambda_WT_CL) #relative increase in variance due to nonresponse

v_old_WT_CL<-(m-1)*(1+1/r_WT_CL^2) #old degree of freedom

v_com=n-k #degree of freedom of parameter estimate (WT_CL) in the hypothetically complete data

v_obs_WT_CL=((v_com+1)/(v_com+3))*v_com*(1-lambda_WT_CL) #observed data degrees of freedom that accounts for the missing information

v_WT_CL<-(v_old_WT_CL*v_obs_WT_CL)/(v_old_WT_CL+v_obs_WT_CL)

t_crit_WT_CL <- qt(0.975, v_WT_CL) 

lower_bound_WT_CL<-WT_CL-t_crit_WT_CL*sqrt(T_WT_CL) #1.425197
upper_bound_WT_CL<-WT_CL+t_crit_WT_CL*sqrt(T_WT_CL) #1.958517

# Here I approximate it with normal distribution to calculate RSE

sqrt(T_WT_CL)/WT_CL*100 #46.64735%

lower_bound_WT_CL
upper_bound_WT_CL
# Final result: 0.6358 (46.64735%)


# Then, pooling WT_V

WT_V<-mean(Pooled_War_30$Estimates[Pooled_War_30$Parameters == "WT_V"]) #1.085667

W_WT_V<-mean(Pooled_War_30$Var[Pooled_War_30$Parameters == "WT_V"]) #0.02165006 #within imputation variance

B_WT_V<-var(Pooled_War_30$Estimates[Pooled_War_30$Parameters == "WT_V"]) #0.001411609 #between imputation variance

WT_V # 1.085667
W_WT_V
B_WT_V

# Total variance

m<-30
T_WT_V<-W_WT_V+(1+1/m)*B_WT_V #0.02310873
T_WT_V

# Confidence interval - T distribution

m<-30 #number of imputations
k<-12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect)
n<-32

lambda_WT_V<-(B_WT_V+B_WT_V/m)/T_WT_V #proportion of variation attributable to the missing data

r_WT_V<-(lambda_WT_V)/(1-lambda_WT_V) #relative increase in variance due to nonresponse

v_old_WT_V<-(m-1)*(1+1/r_WT_V^2) #old degree of freedom

v_com=n-k #degree of freedom of parameter estimate (WT_V) in the hypothetically complete data

v_obs_WT_V=((v_com+1)/(v_com+3))*v_com*(1-lambda_WT_V) #observed data degrees of freedom that accounts for the missing information

v_WT_V<-(v_old_WT_V*v_obs_WT_V)/(v_old_WT_V+v_obs_WT_V)

t_crit_WT_V <- qt(0.975, v_WT_V) 

lower_bound_WT_V<-WT_V-t_crit_WT_V*sqrt(T_WT_V) 
upper_bound_WT_V<-WT_V+t_crit_WT_V*sqrt(T_WT_V) 

# Here I approximate it with normal distribution to calculate RSE

sqrt(T_WT_V)/WT_V*100 # 14.0% 

lower_bound_WT_V
upper_bound_WT_V
# Final result: 1.085667 (14.0%)


## For the 50% missing WT scenario ---------------------------------------------

# Read in the dataset, remove NA
setwd("C:/Users/u0164053/OneDrive - KU Leuven/PhD in KU Leuven/Projects/Multiple Imputation Tutorial in Pharmacometrics/Case study/MI pooling/50% WT missing")
Pooled_War_50<-read.csv("Pooled_War_50.csv",sep=",")
Pooled_War_50<-na.omit(Pooled_War_50) # Here I will remove NAs, which mean paramter estimates from non-converged models

# Creating Var (Variance) column 
Pooled_War_50$Var<-(Pooled_War_50$RSE*Pooled_War_50$Estimates/100)^2


# Then, pooling WT_CL 

WT_CL<-mean(Pooled_War_50$Estimates[Pooled_War_50$Parameters == "WT_CL"]) #0.612551

W_WT_CL<-mean(Pooled_War_50$Var[Pooled_War_50$Parameters == "WT_CL"]) #within imputation variance

B_WT_CL<-var(Pooled_War_50$Estimates[Pooled_War_50$Parameters == "WT_CL"]) #between imputation variance

WT_CL 
W_WT_CL
B_WT_CL

# Total variance

m<-49
T_WT_CL<-W_WT_CL+(1+1/m)*B_WT_CL 
T_WT_CL

# Confidence interval - T distribution

m<-49 #number of imputations
k<-12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect)
n<-32

lambda_WT_CL<-(B_WT_CL+B_WT_CL/m)/T_WT_CL #proportion of variation attributable to the missing data

r_WT_CL<-(lambda_WT_CL)/(1-lambda_WT_CL) #relative increase in variance due to nonresponse

v_old_WT_CL<-(m-1)*(1+1/r_WT_CL^2) #old degree of freedom

v_com=n-k #degree of freedom of parameter estimate (WT_CL) in the hypothetically complete data

v_obs_WT_CL=((v_com+1)/(v_com+3))*v_com*(1-lambda_WT_CL) #observed data degrees of freedom that accounts for the missing information

v_WT_CL<-(v_old_WT_CL*v_obs_WT_CL)/(v_old_WT_CL+v_obs_WT_CL)

t_crit_WT_CL <- qt(0.975, v_WT_CL) 

lower_bound_WT_CL<-WT_CL-t_crit_WT_CL*sqrt(T_WT_CL) 
upper_bound_WT_CL<-WT_CL+t_crit_WT_CL*sqrt(T_WT_CL) 

# Here I approximate it with normal distribution to calculate RSE

sqrt(T_WT_CL)/WT_CL*100 #45.5612%

lower_bound_WT_CL
upper_bound_WT_CL
# Final result: 0.612551 (45.5612%)


# Then, pooling WT_V

WT_V<-mean(Pooled_War_50$Estimates[Pooled_War_50$Parameters == "WT_V"]) #0.8474082

W_WT_V<-mean(Pooled_War_50$Var[Pooled_War_50$Parameters == "WT_V"]) #within imputation variance

B_WT_V<-var(Pooled_War_50$Estimates[Pooled_War_50$Parameters == "WT_V"]) #between imputation variance

WT_V 
W_WT_V
B_WT_V

# Total variance

m<-49
T_WT_V<-W_WT_V+(1+1/m)*B_WT_V 
T_WT_V

# Confidence interval - T distribution

m<-49 #number of imputations
k<-12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect)
n<-32

lambda_WT_V<-(B_WT_V+B_WT_V/m)/T_WT_V #proportion of variation attributable to the missing data

r_WT_V<-(lambda_WT_V)/(1-lambda_WT_V) #relative increase in variance due to nonresponse

v_old_WT_V<-(m-1)*(1+1/r_WT_V^2) #old degree of freedom

v_com=n-k #degree of freedom of parameter estimate (WT_V) in the hypothetically complete data

v_obs_WT_V=((v_com+1)/(v_com+3))*v_com*(1-lambda_WT_V) #observed data degrees of freedom that accounts for the missing information

v_WT_V<-(v_old_WT_V*v_obs_WT_V)/(v_old_WT_V+v_obs_WT_V)

t_crit_WT_V <- qt(0.975, v_WT_V) 

lower_bound_WT_V<-WT_V-t_crit_WT_V*sqrt(T_WT_V) 
upper_bound_WT_V<-WT_V+t_crit_WT_V*sqrt(T_WT_V) 

# Here I approximate it with normal distribution to calculate RSE

sqrt(T_WT_V)/WT_V*100 # 22.5247% 

lower_bound_WT_V
upper_bound_WT_V
# Final result: 0.8474082 (22.5247%)


## For the 80% missing WT scenario ---------------------------------------------

# Read in the dataset, remove NA
setwd("C:/Users/u0164053/OneDrive - KU Leuven/PhD in KU Leuven/Projects/Multiple Imputation Tutorial in Pharmacometrics/Case study/MI pooling/80% WT missing")
Pooled_War_80<-read.csv("Pooled_War_80.csv",sep=",")
Pooled_War_80<-na.omit(Pooled_War_80) # Here I will remove NAs, which mean paramter estimates from non-converged models

# Creating Var (Variance) column 
Pooled_War_80$Var<-(Pooled_War_80$RSE*Pooled_War_80$Estimates/100)^2


# Then, pooling WT_CL 

WT_CL<-mean(Pooled_War_80$Estimates[Pooled_War_80$Parameters == "WT_CL"]) #0.3210117

W_WT_CL<-mean(Pooled_War_80$Var[Pooled_War_80$Parameters == "WT_CL"]) #within imputation variance

B_WT_CL<-var(Pooled_War_80$Estimates[Pooled_War_80$Parameters == "WT_CL"]) #between imputation variance

WT_CL 
W_WT_CL
B_WT_CL

# Total variance

m<-73
T_WT_CL<-W_WT_CL+(1+1/m)*B_WT_CL 
T_WT_CL

# Confidence interval - T distribution

m<-73 #number of imputations
k<-12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect)
n<-32

lambda_WT_CL<-(B_WT_CL+B_WT_CL/m)/T_WT_CL #proportion of variation attributable to the missing data

r_WT_CL<-(lambda_WT_CL)/(1-lambda_WT_CL) #relative increase in variance due to nonresponse

v_old_WT_CL<-(m-1)*(1+1/r_WT_CL^2) #old degree of freedom

v_com=n-k #degree of freedom of parameter estimate (WT_CL) in the hypothetically complete data

v_obs_WT_CL=((v_com+1)/(v_com+3))*v_com*(1-lambda_WT_CL) #observed data degrees of freedom that accounts for the missing information

v_WT_CL<-(v_old_WT_CL*v_obs_WT_CL)/(v_old_WT_CL+v_obs_WT_CL)

t_crit_WT_CL <- qt(0.975, v_WT_CL) 

lower_bound_WT_CL<-WT_CL-t_crit_WT_CL*sqrt(T_WT_CL) #1.425197
upper_bound_WT_CL<-WT_CL+t_crit_WT_CL*sqrt(T_WT_CL) #1.958517

# Here I approximate it with normal distribution to calculate RSE

sqrt(T_WT_CL)/WT_CL*100 #158.6361%

lower_bound_WT_CL
upper_bound_WT_CL
# Final result: 0.3210117 (158.6361%)


# Then, pooling WT_V

WT_V<-mean(Pooled_War_80$Estimates[Pooled_War_80$Parameters == "WT_V"]) #0.9687808

W_WT_V<-mean(Pooled_War_80$Var[Pooled_War_80$Parameters == "WT_V"]) #within imputation variance

B_WT_V<-var(Pooled_War_80$Estimates[Pooled_War_80$Parameters == "WT_V"]) #between imputation variance

WT_V 
W_WT_V
B_WT_V

# Total variance

m<-73
T_WT_V<-W_WT_V+(1+1/m)*B_WT_V 
T_WT_V

# Confidence interval - T distribution

m<-73 #number of imputations
k<-12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect)
n<-32

lambda_WT_V<-(B_WT_V+B_WT_V/m)/T_WT_V #proportion of variation attributable to the missing data

r_WT_V<-(lambda_WT_V)/(1-lambda_WT_V) #relative increase in variance due to nonresponse

v_old_WT_V<-(m-1)*(1+1/r_WT_V^2) #old degree of freedom

v_com=n-k #degree of freedom of parameter estimate (WT_V) in the hypothetically complete data

v_obs_WT_V=((v_com+1)/(v_com+3))*v_com*(1-lambda_WT_V) #observed data degrees of freedom that accounts for the missing information

v_WT_V<-(v_old_WT_V*v_obs_WT_V)/(v_old_WT_V+v_obs_WT_V)

t_crit_WT_V <- qt(0.975, v_WT_V) 

lower_bound_WT_V<-WT_V-t_crit_WT_V*sqrt(T_WT_V) 
upper_bound_WT_V<-WT_V+t_crit_WT_V*sqrt(T_WT_V) 

# Here I approximate it with normal distribution to calculate RSE

sqrt(T_WT_V)/WT_V*100 # 32.35046% 

lower_bound_WT_V
upper_bound_WT_V
# Final result: 0.9687808 (32.35046%)

