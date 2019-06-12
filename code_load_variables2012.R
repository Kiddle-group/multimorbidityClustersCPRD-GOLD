# Last updated: 29 May 2019 #
# This file accompanies the paper "Identification of clusters of multimorbid patients in UK general practice: a latent class analysis."


setwd("C:\\MM_proc_data\\ffdb")

# LIBRARIES #
#**********************************************************************#
source('C:\\MMcluster\\code_misc.R') # load helper functions
#install.packages("CALIBERdatamanage", repos="http://R-Forge.R-project.org")
library(CALIBERdatamanage)
library(ggplot2)        # for generating visualizations
#install.packages("dplyr")
library(dplyr)      # for data manipulation
install.packages("AMR")
library(AMR)      # for data manipulation
source('C:\\MMcluster\\kiddle_misc.R') # load helper functions
require(data.table)
install.packages("devtools")
library(devtools)
source("http://pcwww.liv.ac.uk/~william/R/crosstab.r")
library(reshape2)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.8")
library(ComplexHeatmap)
require(gridExtra)
install.packages('survminer')
source("https://bioconductor.org/biocLite.R")
library(survminer)
library(survival)
#install.packages('ranger')
library(ranger)
#install.packages(c("ggplot2", "ggfortify"))
library(ggfortify)
library(grid)

#*********************************************************************#


#rm(list=ls())
patient2012 <- as.data.table(patient)

tmp <- patient2012[,c('crd','uts')]

patient2012$startid <- apply(tmp, 1, 'max')

patient2012 <- patient2012[patient2012$startid <= '2012-01-01',]

# identify patients who transferred out of practice before start ID, NONE, GREAT!

out_before_in <- which(patient2012$startid > as.Date(patient2012$tod,format='%d/%m/%Y')) 

rm(out_before_in) #remove tmp files

# age 
  patient2012$age <- 212 - patient2012$yob

# load and merge IMD data

  imd <- read.table('Q:\\kiddlegroup\\MM clustering project\\!Raw Linkage Data\\patient_imd2015_16_057RA.txt',sep='\t',header=T)
  patient2012 <- merge(patient2012,imd[,-2],by='patid',all.x=T)

#smoking, make an exception, use memory instead of hard drive to avoid strange error
# smokestatus final order: 1= current, 2= never, 3 = ex
  
  # find instances of medcodes
  codes <- read.dta13('Q:\\kiddlegroup\\CodelistsV1.1\\smoking_codes.dta')[,c(1,3)]
  relevant <- as.data.table(find_medcodes(codes,lm_time=2012))
  
  # recode to make consistent with additional data
  codes[codes[,"smokstatus"] == 1,"smokstatus"] <- 3
  codes[codes[,"smokstatus"] == 2,"smokstatus"] <- 1
  
  # add data to relevant
  relevant <- merge(relevant,codes,by='medcode')
  relevant$post_lm <- as.numeric(relevant$eventdate - as.Date('2012-01-01'))
  
  # find instances use enttype and linked additional data
  
  clinical_relevant2 <- as.data.table.ffdf(clinical2012_before[clinical2012_before$enttype == 4,c('patid','adid','eventdate')])
  additional_relevant <- as.data.table.ffdf(additional[additional$enttype==4,c('patid','adid','data1')])
  colnames(additional_relevant)[3] <- 'smokstatus'
  
  additional_relevant$smokstatus <- as.numeric(as.character(additional_relevant$smokstatus))
  
  # use clinical enttype to identify presence of additional data, then track it down
  clinical_relevant2 <- merge(clinical_relevant2,additional_relevant,by=c('patid','adid'))
  
  clinical_relevant2$post_lm <- as.numeric(clinical_relevant2$eventdate - as.Date('2012-01-01'))
  
  # merge data derived using medcodes and enttype
  smoke <- rbind(relevant[,c('patid','smokstatus','post_lm')],clinical_relevant2[,c('patid','smokstatus','post_lm')])
  
  rm(relevant,clinical_relevant2)
  
  # neat trick to extract most recent data for each patient
  setorderv(smoke,c('patid','post_lm'),order=c(1,-1))
  smoke$order <- ave(rep(1,nrow(smoke)),smoke[,1],FUN=seq_along )
  smoke <- smoke[smoke$order==1,]
  
  # add most recent smoking data to patient dataset
  patient2012 <- merge(patient2012,smoke[,c('patid','smokstatus')],by='patid',all.x=T)
  
  rm(smoke, codes)

# BMI
  
  # find medcode instances
  clinical_relevant <- as.data.table.ffdf(clinical2012_before[clinical2012_before$enttype == 13,c('patid','eventdate','adid','enttype')])
  additional_relevant <- as.data.table.ffdf(additional[additional$enttype==13,c('patid','adid','data3')])
  bmi_both <- merge(clinical_relevant,additional_relevant,by=c('patid','adid'))
  bmi_both$data3 <- as.numeric(as.character(bmi_both$data3))
  bmi_both <- bmi_both[!is.na(bmi_both$data3)]
  bmi_both <- bmi_both[bmi_both$data3<=70,]# BMI > 70 possible but very unlikely (max val 325000!!)
  
  bmi_both <- find_latest(bmi_both,ffdf=F)
  colnames(bmi_both)[5] <- 'bmi'
  
  patient2012 <- merge(patient2012,bmi_both[,c('patid','bmi')],by='patid',all.x=T)
  
  rm(clinical_relevant,additional_relevant,bmi_both)

  
# need to identify co-morbidities 
#
# using codes and algorithms from https://www.phpc.cam.ac.uk/pcu/cprd_cam/codelists/v11/
#

print('...identify comorbid')

# CKD
  print('..extract eGFR for CKD')
  
  # find egfr data
  egfr <- test2012_before[test2012_before$enttype==466,c('patid','eventdate','data1','data2')] 
  
  #load('../output/proc_data/egfr.RData')
  
  #data1 - Operator (OPR) 0, 1 < , 2 <=,3 =, 4 >, 5 >=, 6 ~
  #data2 - eGFR Value
  #data3 - Unit of measure (SUM) - 90 is mL/min
  #data4 - Qualifier	(TQU) NA
  #data5 - Normal range from NA	
  #data6 - Normal range to NA
  #data7 - Normal range basis (POP) NA
  
  egfr <- egfr[!duplicated(egfr),] # remove duplicates
  
  egfr <- egfr[!is.na(egfr$data2),] # remove missing values
  
  # remove rare confusing cases
  ind <- !(egfr$data1==1 & egfr$data2==60) 
  egfr <- egfr[ind,] 
  
  # remove unclassifiable cases
  ind <- !(egfr$data1==1 & egfr$data2==90) 
  egfr <- egfr[ind,] 
  ind <- !(egfr$data1==2 & (egfr$data2==90 | egfr$data2==60 )) 
  egfr <- egfr[ind,] 
  
  # remove missing data
  ind <- !(egfr$data2==0 ) 
  egfr <- egfr[ind,] 
  
  # remove unrealistic range data (i.e.error in records)
  ind <- !(egfr$data2>200 ) 
  egfr <- egfr[ind,] 
  
  egfr <- as.data.table.ffdf(egfr)
  
  setorder(egfr,patid,-eventdate)
  egfr[,'test'] <- ave(rep(1,nrow(egfr)),egfr[,1],FUN=seq_along )
  egfr[,'num_test'] <- ave(rep(1,nrow(egfr)),egfr[,1],FUN=sum )
  
  egfr <- egfr[egfr$num_test > 1,]
  egfr <- egfr[egfr$test < 3,]
  
  # neat trick to get max of all values for individual
  egfr[,'max'] <- ave(egfr$data2,egfr[,1],FUN=max )
  egfr <- egfr[egfr$test == 1,]
  
  # extract egfr, and whether above or below relevant cutpoints
  egfr$egfr <- egfr$max
  egfr$c_ckd_copd <- egfr$max < 60 # used in CMS
  egfr$cci_ckd_copd <- egfr$max < 30 # used in CCI
  
  # add egfr/ckd data to patient dataset
  patient2012<- merge(patient2012,egfr[,c('patid','egfr','c_ckd_copd','cci_ckd_copd')],by='patid',all.x=T)
  patient2012$c_ckd_copd[is.na(patient2012$c_ckd_copd)] <- F
  patient2012$cci_ckd_copd[is.na(patient2012$cci_ckd_copd)] <- F
  
  rm(egfr)
  
  # in case of crashes, save progress
  save.image(file='C:\\MM_proc_data\\beforesimpledef.RData')

# for each simple defintion, find first instances and add data to patient dataset
print('..simple definitions')

# alcohol problems: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_ALC138_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_ap_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# anorexia bulimia: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_ANO139_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_ab_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# atrial fibrillation: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_ATR143_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_atr_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# blindness and low vision:read code ever recoded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_BLI144_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_bli_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# bronchiectasis: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_BRO145_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_bro_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# cancer 1.1 is moved to later!

# chronic liver disease: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_CLD148_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_cld_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# chronic sinusitis: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_SIN149_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_sin_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data
  
# COPD: read code ever recorded
  
  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_COP151_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_copd_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# coronary heart disease: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_CHD126_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_chd_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# dementia: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_DEM131_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_dem_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# diabetes: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_DIB128_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_dia_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# diverticular disease of intestine: read coded ever recorded
  
  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_DIV154_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_div_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# hearing loss: read code every recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_HEL157_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_hel_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Heart failure: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_HEF158_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_hf_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Hypertension: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_HYP159_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_hyp_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Inflammatory bowel disease (IBD): read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_IBD160_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_ibd_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Learning disability: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_LEA163_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_lea_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Multiple sclerosis: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_MSC165_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_ms_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Parkinson's disease: read code ever recorded
  
  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_PRK169_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_prk_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data
  
# Peptic UIcer Disease: read code ever recorded
  
  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_PEP135_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_pud_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data  
  
  
# Peripheral vascular disease: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_PVD168_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_pvd_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Prostate disorders : read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_PRO170_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_pro_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Psychoactive substance misuse (NOT ALCOHOL): read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_PSM173_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_ops_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Rheumatoid arthritis: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_RHE174_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_rhe_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Stroke & transient ischaemic attack: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_STR130_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_str_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# Thyroid disorders: read code ever recorded

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_THY179_MC_V1-1_Oct2018.csv') # load codes
  first <- ever_recorded(codes,lm_time=2012) # find first instance per patient
  colnames(first)[2] <- "c_thy_date" # date of first instance per patient
  patient2012 <- merge(patient2012,first[,c(1,2)],by='patid',all.x = TRUE)
  rm(codes,first) #remove tmp data

# in case of crashes, save progress
  save.image(file='C:\\MM_proc_data\\2012_readcodeeverrecordeddone.RData')
    

###################################################################
#print('..simple, by COPD diagnosis')
# legacy code

# for ever recorded co-morbidities, turn date into above
#tmp <- vapply(patient2012[,33:56],dateNumeric,numeric(dim(patient2012)[1]))

# dates after first COPD date?
#tmp2 <- tmp < rep(dateNumeric(patient2012[,'first_copd']),24)

# set missing to FALSE as well, means never mention of readcode
#tmp2[which(is.na(tmp2),arr.ind=T)] <- FALSE

# extract co-morbidity labels from column names of tmp2
#tmp_mat <- matrix(unlist(strsplit(colnames(tmp2),'_')),3,24)

# new labels to show indicator of co-morbidity by first COPD date
#colnames(tmp2) <- paste('c',tmp_mat[2,],sep='_')
###################################################################

# make column names for morbidities  
#patient2012 <- subset(patient2012, select = -c(c_can)) # drop c_can column, will add later 

colnum <- 32

# identify comorbidities with a missing date, i.e. that haven't been diagnosed
tmp <- vapply(patient2012[,(colnum+1):(colnum+27)],is.na,numeric(dim(patient2012)[1]))

# use comorbidity date column names to make comorbidity present column names
tmp_mat <- matrix(unlist(strsplit(colnames(tmp),'_')),3,27)
colnames(tmp) <- paste('c',tmp_mat[2,],sep='_')

# add to dataset if not missing, i.e. if a diagnosis has been made
patient2012 <- cbind(patient2012,!tmp)

patient2012 <- subset(patient2012, select = -c(c_ckd)) # drop c_ckd column (wrong entry)
colnames(patient2012)[colnames(patient2012) =="c_ckd_copd"] <- "c_ckd"

rm(tmp,tmp_mat)

# cancer v1.1, any read code in last 5 years

  # identify instances, add info
  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_CAN146_MC_V1-1_Oct2018.csv')
  relevant <- as.data.table.ffdf(find_medcodes(codes,lm_time=2012))
  # relevant <- merge(relevant,all_copd,by='patid') #all_copd does not exist - no merge executed
  
  # set to TRUE if recorded in last 5 years
  relevant$before_diagnosis <- as.numeric(as.Date('2012-01-01') - relevant$eventdate)
  relevant <- relevant[(relevant$before_diagnosis<=5*365),]
  relevant <- relevant[!duplicated(relevant[,'patid']),]
  relevant$c_can <- T
  
  patient2012 <- merge(patient2012,relevant[,c('patid','c_can')],by='patid',all.x=T)
  patient2012$c_can[is.na(patient2012$c_can)] <- F

# save in case of crash
save.image(file='C:\\MM_proc_data\\readcodeeverrecordeddone.RData')

print('..complex definitions')
# the rest of the comorbidities are not 'ever recorded', and so require more work to identify
# for each comorbidity we give the definition at the beggining

# Anxiety (read code in last 12 months OR >= 4 anxiolytic/hypnotic presciption in last 12 months)
  
  # load med/prod codes
  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_ANX140_MC_V1-1_Oct2018.csv') 
  prods <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_ANX141_PC_V1-1_Oct2018.csv')
  
  # find instances of medcodes
  relevant <- as.data.table.ffdf(find_medcodes(codes,lm_time=2012))
  relevant <- relevant[relevant$eventdate >= '2011-01-01',]
  #relevant <- find_latest(relevant,ffdf=F)
  
  relevant$c_anx <- T 
  
  # find instances of prodcodes

  therapy_relevant<-bigTherapyExtract(prods,chunkSize=10^7)
  therapy_relevant$before_diagnosis <- as.ff(as.numeric(as.Date('2012-01-01') - as.ram(therapy_relevant$eventdate)))
  # require to be within 12 months
  therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])
  
  # how many per patient? need >=4 
  therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
  therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
  therapy_relevant$c_anx <- T
  
  # combine as can be either, but leave unique patids and then merge
  rels <- rbind(relevant[,c('patid','c_anx')],therapy_relevant[,c('patid','c_anx')])
  rels <- rels[!duplicated(rels[,'patid']),]
  
  patient2012 <- merge(patient2012,rels,by='patid',all.x=T)
  patient2012$c_anx[is.na(patient2012$c_anx)] <- F


# Asthma (read code ever recorded AND presciption in last 12 months)
# Quint codes and rules, basically COPD and asthma meds overlap, and 0-2 yr likely misdiagnosed 
  #codes <- read.dta13('../data/other_files/codes/asthma.dta') # this is for COPD project, do not use for MMcluster
  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_AST142_MC_V1-1_Oct2018.csv') # CPRD @ cambridge codes, use this
  prods <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_AST127_PC_V1-1_Oct2018.csv')
  
  # keep only rows with right codes
  relevant <- as.data.table.ffdf(find_medcodes(codes,lm_time=2012))
  relevant <- relevant[!duplicated(relevant[,'patid']),]
  
  # num days between landmark time and read code recording
  #relevant$before_diagnosis <- as.numeric(as.Date('2012-01-01') - relevant$eventdate)
  #relevant <- relevant[(relevant$before_diagnosis>=0)&(relevant$before_diagnosis<=1*365),]

  # treatment in last year
  therapy_relevant<-bigTherapyExtract(prods,chunkSize=10^7)
  therapy_relevant$before_diagnosis <- as.ff(as.numeric(as.Date('2012-01-01') - as.ram(therapy_relevant$eventdate)))
  therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])
  therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  
  # combine treatment and diagnosis info and add to dataset
  rels <- merge(relevant,therapy_relevant,by='patid')
  rels$c_ast <- T
  
  patient2012 <- merge(patient2012,rels[,c('patid','c_ast')],by='patid',all.x=T)
  patient2012$c_ast[is.na(patient2012$c_ast)] <- F

# Constipation - 4 or more laxative prescriptions in last year
  prods <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_CON150_PC_V1-1_Oct2018.csv')
  
  # identify instances within last year, count and require at least 4
  therapy_relevant<-bigTherapyExtract(prods,chunkSize=10^7)
  therapy_relevant$before_diagnosis <- as.ff(as.numeric(as.Date('2012-01-01') - as.ram(therapy_relevant$eventdate)))
  therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])
  
  therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
  therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
  therapy_relevant$c_con <- T
  
  patient2012 <- merge(patient2012,therapy_relevant[,c('patid','c_con')],by='patid',all.x=T)
  patient2012$c_con[is.na(patient2012$c_con)] <- F


# Depression (read code in last 12 months OR >= 4 presciption in last 12 months)

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_DEP152_MC_V1-1_Oct2018.csv') # load codes
  prods <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_DEP153_PC_V1-1_Oct2018.csv')
  
  # keep only rows with right codes
  relevant <- as.data.table.ffdf(find_medcodes(codes,lm_time=2012))
  #relevant <- merge(relevant,by='patid') #no all_copd - no merge
  
  # require to be within year of landmark time
  relevant$before_diagnosis <- as.numeric(as.Date('2012-01-01') - relevant$eventdate)
  relevant <- relevant[relevant$before_diagnosis<=365,]
  relevant$c_dep <- T 
  
  # identify treatment instances and require to be within a year
  therapy_relevant<-bigTherapyExtract(prods,chunkSize=10^7)
  therapy_relevant$before_diagnosis <- as.ff(as.numeric(as.Date('2012-01-01') - as.ram(therapy_relevant$eventdate)))
  therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])
  
  # count and require to be at least 4
  therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
  therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
  therapy_relevant$c_dep <- T
  
  # combine as can be either, but leave unique patids
  rels <- rbind(relevant[,c('patid','c_dep')],therapy_relevant[,c('patid','c_dep')])
  rels <- rels[!duplicated(rels[,'patid']),]
  
  patient2012 <- merge(patient2012,rels,by='patid',all.x=T)
  patient2012$c_dep[is.na(patient2012$c_dep)] <- F


# Epilepsy (read code ever recorded AND presciption in last 12 months)

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_EPI155_MC_V1-1_Oct2018.csv') # load codes
  prods <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_EPI156_PC_V1-1_Oct2018.csv')
  
  # keep only rows with right codes
  relevant <- as.data.table.ffdf(find_medcodes(codes,lm_time=2012))
  relevant <- relevant[!duplicated(relevant[,'patid']),]
  
  # identify treatment instances in year before landmark time
  therapy_relevant<-bigTherapyExtract(prods,chunkSize=10^7)
  therapy_relevant$before_diagnosis <- as.ff(as.numeric(as.Date('2012-01-01') - as.ram(therapy_relevant$eventdate)))
  therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])
  
  # any?
  therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
  therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  
  # merge, i.e. if not in both will have no entry
  both <- merge(relevant,therapy_relevant,by='patid')
  both$c_epi <- T
  patient2012 <- merge(patient2012,both[,c('patid','c_epi')],by='patid',all.x=T)
  patient2012$c_epi[is.na(patient2012$c_epi)] <- F


# painful condition (4 or more prescriptions in last 12 months OR (4 or more anti-epleptics in absence of elilepsy ever recorded) 
  
  prods <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_PNC166_PC_V1-1_Oct2018.csv')
  
  # identify treatment in the preceding year of the landmark time using new list
  therapy_relevant2 <- bigTherapyExtract(prods,chunkSize=10^7)
  therapy_relevant2$before_diagnosis <- as.ff(as.numeric(as.Date('2012-01-01') - as.ram(therapy_relevant2$eventdate)))
  therapy_relevant2 <- as.data.table.ffdf(therapy_relevant2[therapy_relevant2$before_diagnosis <= 365,])
  
  # any? either in this list or previous: 4 or more pain pain prods / (4 or more anti-elipsy & no elipesy medcode)
  therapy_relevant2$num <- ave(rep(1,nrow(therapy_relevant2)),therapy_relevant2[,1],FUN=sum )
  therapy_relevant2 <- therapy_relevant2[!duplicated(therapy_relevant2[,'patid']),]
  therapy_relevant2 <- therapy_relevant2[therapy_relevant2$num>=4,]
  therapy_relevant2$c_pnc <- T
  
  therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
  therapy_relevant <- therapy_relevant[!(therapy_relevant$patid %in% relevant$patid)]
  therapy_relevant$c_pnc <- T
  
  # combine as can be either, but leave unique patids
  rels <- rbind(therapy_relevant[,c('patid','c_pnc')],therapy_relevant2[,c('patid','c_pnc')])
  rels <- rels[!duplicated(rels[,'patid']),]
  patient2012 <- merge(patient2012,rels[,c('patid','c_pnc')],by='patid',all.x=T)
  patient2012$c_pnc[is.na(patient2012$c_pnc)] <- F


# IBS: irritable bowel syndrome (read code ever recorded OR 4 or more presciption in last 12 months)
  
  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_IBS161_MC_V1-1_Oct2018.csv') # load codes
  prods <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_IBS162_PC_V1-1_Oct2018.csv')
  
  # keep only rows with right codes
  relevant <- as.data.table.ffdf(find_medcodes(codes,lm_time=2012))
  relevant$c_ibs <- T 
  
  # identify treatment in preceding year
  therapy_relevant<-bigTherapyExtract(prods,chunkSize=10^7)
  therapy_relevant$before_diagnosis <- as.ff(as.numeric(as.Date('2012-01-01') - as.ram(therapy_relevant$eventdate)))
  therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])
  
  # four or more?
  therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
  therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
  therapy_relevant$c_ibs <- T
  
  # combine as can be either, but leave unique patids
  rels <- rbind(relevant[,c('patid','c_ibs')],therapy_relevant[,c('patid','c_ibs')])
  rels <- rels[!duplicated(rels[,'patid']),]
  
  patient2012 <- merge(patient2012,rels,by='patid',all.x=T)
  patient2012$c_ibs[is.na(patient2012$c_ibs)] <- F


# Migrane - 4 or more migrane prescriptions in last year

  prods <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_MIG164_PC_V1-1_Oct2018.csv')
  
  # identify tratment instances
  therapy_relevant<-bigTherapyExtract(prods,chunkSize=10^7)
  therapy_relevant$before_diagnosis <- as.ff(as.numeric(as.Date('2012-01-01') - as.ram(therapy_relevant$eventdate)))
  therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])
  
  # four or more?
  therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
  therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
  therapy_relevant$c_mig <- T
  
  patient2012 <- merge(patient2012,therapy_relevant[,c('patid','c_mig')],by='patid',all.x=T)
  patient2012$c_mig[is.na(patient2012$c_mig)] <- F


# psoriasis or eczema (read code ever recorded AND 4 or more prescriptions in last 12 months) 
  
  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_PSO171_MC_V1-1_Oct2018.csv') # load codes
  prods <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_PSO172_PC_V1-1_Oct2018.csv')
  
  # keep only rows with right codes
  relevant <- as.data.table(find_medcodes(codes,lm_time=2012))
  relevant <- relevant[!duplicated(relevant[,'patid']),]
  
  # identify treatment instances in preceding year
  therapy_relevant<-bigTherapyExtract(prods,chunkSize=10^7)
  therapy_relevant$before_diagnosis <- as.ff(as.numeric(as.Date('2012-01-01') - as.ram(therapy_relevant$eventdate)))
  therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])
  
  # >= 4 prescriptions
  therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
  therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
  
  # merge diagnoses and treatments, so only patients with both remain
  both <- merge(relevant,therapy_relevant,by='patid')
  both$c_pso <- T
  patient2012 <- merge(patient2012,both[,c('patid','c_pso')],by='patid',all.x=T)
  patient2012$c_pso[is.na(patient2012$c_pso)] <- F
  

#AECOPD only happen after COPD diagnosis, so don't worry about for now
#aecopd <- read.dta13('../../data/other_files/AECOPD/all_aecopd2.dta')
#aepatient2012 <- aecopd[aecopd$patid %in% patient2012$patid,]
#aecopd_12m <- aepatient2012[which((aepatient2012$eventdate <= '2013-01-01') & (aepatient2012$eventdate >= '2012-01-01')),]
#aecopd_12m_hes <- aecopd_12m[aecopd_12m$hes_aecopd>0,]

#MRC, dyspnea occurs rarely before diagnosis
#
#codes <- read.xlsx('../../data/other_files/codes/Excel codelists/MRC.xlsx',1) # load codes (COPD only)
# relevant <- find_medcodes(codes,lm_time=2012) # find first instance per patient
# relevant <- find_latest(relevant)
# rownames(codes) <- codes$medcode
# mrc_data <- codes[as.character(relevant$medcode[-1]),3:7]
# mrc_data$patid <- relevant$patid[-1]
# patient2012 <- merge(patient2012,mrc_data,by='patid',all.x=T)


#FEV1 in last 12 months, remove missing FEV1 and deal with units
  # fev1 <- read.dta13('../../data/other_files/FEV/all_fev1.dta') # data from Quint group
  # fev1_nona <- fev1[!is.na(fev1$fev1),]
  # fev1_nona$fev1[fev1_nona$data3 == 89] <- fev1_nona$fev1[fev1_nona$data3 == 89]/1000
  # fev1_nona <- fev1_nona[fev1_nona$fev1 < 5,] # remove outlying values
  # 
  # # extract just values from eligible patients within year before COPD diagnosis
  # fev1_eligible <- fev1_nona[fev1_nona$patid %in% patient2012$patid,]
  # fev1_eligible$post_lm<- dateNumeric(fev1_eligible$eventdate) - dateNumeric('2012-01-01')
  # fev1_before <- fev1_eligible[(fev1_eligible$post_lm <= 0) & (fev1_eligible$post_lm >= -365),]
  # 
  # # reorder by patient, days before diagnosis and highest reading for the day
  # fev1_before <- fev1_before[order(fev1_before$patid,fev1_before$post_diagnosis,fev1_before$fev1,decreasing = c(F,T,T)),]
  # 
  # # nice trick to give per patient order
  # fev1_before$order <- ave(rep(1,nrow(fev1_before)),fev1_before[,1],FUN=seq_along )
  # 
  # # extract latest observation
  # fev1_before <- fev1_before[fev1_before$order==1,]
  # 
  # # extract FEV1 percent predicted (i.e. airflow obstruction)
  # fev1_before$airflow <- fev1_before$percent_fev
  # 
  # # merge with patient data
  # patient2012 <- merge(patient2012,fev1_before[,c('patid','airflow')],by='patid',all.x=T)
  # 
  # # set up not recorded indicator
  # missing <- is.na(patient2012$airflow)
  # patient2012$airflow[missing] <- 0
  # patient2012$airflowMiss <- missing
  # 
  # rm(fev1,fev1_eligible,fev1_nona,fev1_before,missing) # remove tmp files

# scizophrenia (read code ever recorded OR lithium (a prodcode) ever recorded) 

  codes <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_SCZ175_MC_V1-1_Oct2018.csv') # load codes
  prods <- read.csv('Q:\\kiddlegroup\\CodelistsV1.1\\CPRDCAM_SCZ176_PC_V1-1_Oct2018.csv')
  
  # medcodes
  relevant <- as.data.table(find_medcodes(codes,lm_time=2012))
  relevant <- relevant[!duplicated(relevant[,'patid']),]
  relevant$c_scz <- T
  
  # prodcodes
  therapy_relevant<-bigTherapyExtract(prods,chunkSize=10^7)
  therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  therapy_relevant$c_scz <- T
  
  # either
  rels <- rbind(relevant[,c('patid','c_scz')],therapy_relevant[,c('patid','c_scz')])
  rels <- rels[!duplicated(rels[,'patid']),]
  
  patient2012 <- merge(patient2012,rels,by='patid',all.x=T)
  patient2012$c_scz[is.na(patient2012$c_scz)] <- F

rm(both,codes,imd,ind,linked,missing,prods,relevant,rels,therapy_relevant,therapy_relevant2)

save(patient2012,file='C:\\MM_proc_data\\variables.RData',replace)
save.image(file='C:\\MM_proc_data\\allmorbiditiesready.RData')

# *****************************************************************************************#

##DATA CHECKING## 

# Prevalence: percentage of people with each condition

setwd("C:\\MM_proc_data")

load("C:\\MM_proc_data\\variables.RData")

head(patient2012)
dim(patient2012)

#note gender: 1=male, 2=female

# 38 multimorbidiies in the dataset (using CPRD Oct 2018 data)

# *****************************************************************************************#
# Find potential predictors and outcomes #

# Covariates
summary (patient2012$age)
prop.table(table(patient2012$imd2015_5))
prop.table(table(patient2012$smokstatus))
prop.table(table(patient2012$gender))
table(patient2012$gender)

# Outcome: Treatment burden: total number of prodcodes prescribed at least 3 OR 4 times in 1 yr
#         Note adjustments (version2): count unique drug ingredients using BNF, ignoring packsizes (redo therapy_2012 before and after files)
# first load therapy2012_after data into R

head(therapy2012_after)
  # patient2012 <- subset(patient2012, select = -c(ctprod12mons.x, ctprod12mons.y, ctprod12mons))
  # therapy_relevant<-subset(therapy_relevant,select=-c(ctprod12mons ctprod6mons))
  
  # 12 mons after landmark time:  total nubmer of unique prescriptions
  
  # #**************************************************************#
  # subset test
  therapy_relevanttoy <- therapy_relevant[1:500,] 
  therapy_relevanttoy <-therapy_relevanttoy %>%                   
                         group_by(patid,BNF) %>%          
                         mutate(ct12mons_BNF = n()) 
  therapy_relevanttoy <-therapy_relevanttoy %>%                   
                         group_by(patid) %>%          
                         mutate(ct12mons4plus = uniqueN(BNF[ct12mons_BNF>=4]))
  table(therapy_relevanttoy$ct12mons4plus)
  
  # #**************************************************************#
# VERSION1: use unique prodcodes

  therapy_relevant <- as.data.table.ffdf(therapy2012_after[therapy2012_after$eventdate >= '2012-01-01'&therapy2012_after$eventdate <= '2013-01-01',])
  therapy_relevant<-therapy_relevant %>%                   
                         group_by(patid,prodcode) %>%          
                         mutate(ctprod12mons = n()) 
  therapy_relevant <-therapy_relevant %>%                   
                         group_by(patid) %>%          
                         mutate(ctprod12mons4plus = uniqueN(prodcode[ctprod12mons>=4]))  
   
  therapy_relevant1 <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  patient2012<- merge(patient2012,therapy_relevant1[,c('patid','ctprod12mons4plus')],by='patid',all.x=T)
  patient2012$ctprod12mons4plus[is.na(patient2012$ctprod12mons4plus)]<-0
  
  table(patient2012$ctprod12mons4plus)
  hist(patient2012[patient2012$ctprod12mons4plus!=0]$ctprod12mons4plus)

#Version 2: use BNF codes for unique identification
  therapy_relevant<-therapy_relevant %>%                   
    group_by(patid,BNF) %>%          
    mutate(ctprod12mons_BNF = n()) 
  therapy_relevant <-therapy_relevant %>%                   
    group_by(patid) %>%          
    mutate(ctprod12mons4plus_BNF = uniqueN(BNF[ctprod12mons_BNF>=4]))  
  
  therapy_relevant2 <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
  patient2012<- merge(patient2012,therapy_relevant2[,c('patid','ctprod12mons4plus_BNF')],by='patid',all.x=T)
  patient2012$ctprod12mons4plus_BNF[is.na(patient2012$ctprod12mons4plus_BNF)]<-0
  
  table(patient2012$ctprod12mons4plus_BNF)
  hist(patient2012$ctprod12mons4plus_BNF)
  
  
# OUtcome: Total number of GP consultations? (consultation file not avaiable - installed now)
  # 12 mons after landmark time:  total nubmer of GP consultations (count number of unique days)
  consultation_relevant <- as.data.table.ffdf(consultation2012_after[consultation2012_after$eventdate >= '2012-01-01'&consultation2012_after$eventdate <= '2013-01-01',])
  consultation_relevant[,'ctconsult12mons'] <- as.numeric(ave(as.character(consultation_relevant$eventdate), as.factor(consultation_relevant$patid), FUN=function(x) length(unique(x))))
  consultation_relevant <- consultation_relevant[!duplicated(consultation_relevant[,'patid']),]
  patient2012<- merge(patient2012,consultation_relevant[,c('patid','ctconsult12mons')],by='patid',all.x=T)
  patient2012$ctconsult12mons[is.na(patient2012$ctconsult12mons)]<-0
  
  table(patient2012$ctconsult12mons)
  hist(patient2012$ctconsult12mons)
  
  
  # 6 mons after landmark time:  total nubmer of GP consultations (count number of unique days)
  consultation_relevant <- as.data.table.ffdf(consultation2012_after[consultation2012_after$eventdate >= '2012-01-01'&consultation2012_after$eventdate <= '2012-06-30',])
  consultation_relevant[,'ctconsult6mons'] <- as.numeric(ave(as.character(consultation_relevant$eventdate), as.factor(consultation_relevant$patid), FUN=function(x) length(unique(x))))
  consultation_relevant <- consultation_relevant[!duplicated(consultation_relevant[,'patid']),]
  patient2012<- merge(patient2012,consultation_relevant[,c('patid','ctconsult6mons')],by='patid',all.x=T)
  patient2012$ctconsult6mons[is.na(patient2012$ctconsult6mons)]<-0
  
  table(patient2012$ctconsult6mons)
  hist(patient2012$ctconsult6mons)
  

# Outcome: Hospital admission2: Total number of hospital admission spells in 1 yr / 6 mons (treat discharge as admission date)
  hosp<- read.table('Q:\\kiddlegroup\\MM clustering project\\!Raw Linkage Data\\HES_diagnosis_hosp_integrated_16_057RA.txt',sep='\t',header=T)
  hosp$discharged <- as.Date(hosp$discharged,format='%d/%m/%Y')
  hosp_relevant <- hosp[hosp$discharged >='2012-01-01' & hosp$discharged<='2012-12-31',]
  
  hosp_relevant<-hosp_relevant %>%                   
    group_by(patid) %>%          
    mutate(cthospspl_12mons = uniqueN(spno)) 
  
  tmp <- hosp_relevant[!duplicated(hosp_relevant[,'patid']),]
  
  patient2012 <- merge(patient2012,tmp[,c("patid","cthospspl_12mons")],by='patid',all.x=T)
  patient2012$cthospspl_12mons[is.na(patient2012$cthospspl_12mons)]<-0 # patients with no information on HES data DO NOT have hospitalisation episodoes

  rm(tmp)
  
# Outcomes: Median number of months between consecutive consultations in 6 months and 1 yr 
#  (this may show the "intensity of health service needs"
# patients with only 0 OR 1 consultation in the window has interconsult=0   
  consultation_relevant <- as.data.table.ffdf(consultation2012_after[consultation2012_after$eventdate >= '2012-01-01'&consultation2012_after$eventdate <= '2012-12-31',])
  consultation_relevant<-arrange(consultation_relevant,patid,eventdate)
  
  consultation_relevant<-consultation_relevant %>%
        group_by(patid) %>%
        mutate(interconsult12mons = as.numeric(eventdate - lag(eventdate)))
  summary(consultation_relevant$interconsult12mons)
  
  consultation_relevant<-consultation_relevant %>%
    group_by(patid) %>%
    mutate(count12mons = uniqueN(eventdate)) # used unique to avoid double couting for mulitple consultations in one day
  
  consultation_relevant$interconsult12mons <- ifelse(consultation_relevant$count12mons<=1 & is.na(consultation_relevant$interconsult12mons), 0, consultation_relevant$interconsult12mons)   
  # after the previous line, NA is only present for the first entry of each patient, which won't hurt dist.
  
  consultation_relevant<-consultation_relevant %>%
    group_by(patid) %>%
    mutate(interconsult12monsmed = median(interconsult12mons,na.rm = TRUE))
  summary(consultation_relevant$interconsult12monsmed)
  
  consultation_relevant <- consultation_relevant[!duplicated(consultation_relevant[,'patid']),]
  patient2012<- merge(patient2012,consultation_relevant[,c('patid','interconsult12monsmed')],by='patid',all.x=T)
  patient2012$interconsult12monsmed[is.na(patient2012$interconsult12monsmed)]<-0
  patient2012$interconsult12monsmed<-patient2012$interconsult12monsmed/30 # transform to months
  hist(patient2012$interconsult12monsmed)
  
  ###Now work on 6 month###
  consultation_relevant <- as.data.table.ffdf(consultation2012_after[consultation2012_after$eventdate >= '2012-01-01'&consultation2012_after$eventdate <= '2012-6-30',])
  consultation_relevant<-arrange(consultation_relevant,patid,eventdate)
  
  consultation_relevant<-consultation_relevant %>%
    group_by(patid) %>%
    mutate(interconsult6mons = as.numeric(eventdate - lag(eventdate)))
  summary(consultation_relevant$interconsult6mons)
  
  consultation_relevant<-consultation_relevant %>%
    group_by(patid) %>%
    mutate(count6mons = uniqueN(eventdate)) # used unique to avoid double couting for mulitple consultations in one day
  
  consultation_relevant$interconsult6mons <- ifelse(consultation_relevant$count6mons<=1 & is.na(consultation_relevant$interconsult6mons), 0, consultation_relevant$interconsult6mons)   
  # after the previous line, NA is only present for the first entry of each patient, which won't hurt dist.
  
  consultation_relevant<-consultation_relevant %>%
    group_by(patid) %>%
    mutate(interconsult6monsmed = median(interconsult6mons,na.rm = TRUE))
  summary(consultation_relevant$interconsult6monsmed)
  
  consultation_relevant <- consultation_relevant[!duplicated(consultation_relevant[,'patid']),]
  patient2012<- merge(patient2012,consultation_relevant[,c('patid','interconsult6monsmed')],by='patid',all.x=T)
  patient2012$interconsult6monsmed[is.na(patient2012$interconsult6monsmed)]<-0
  patient2012$interconsult6monsmed<-patient2012$interconsult6monsmed/30 # transform to months
  
  hist(patient2012$interconsult6monsmed)
  
  # drop individuals with gender =3 (neithermale nor female)
  patient2012<-patient2012[!(patient2012$gender=="3"),]
  
  # exclude 271 individuals with IMD ineligible
  patient2012<-patient2012[!is.na(patient2012$imd2015_5),]
  
# Outcome: death
  death<- read.table('Q:\\kiddlegroup\\MM clustering project\\!Raw Linkage Data\\set 16\\GOLD_linked\\death_patient_16_057RA2.txt',sep='\t',header=T)
  death$dod <- as.Date(death$dod,format='%d/%m/%Y')
  baseline<- as.Date("2012-01-01")
  death$dur<- round((death$dod - baseline)/30) # duration in months
  death<-death[death$dur>0,] # keep only individuals who died after 2012
  death$d2012<-(death$dod<="2012-12-31") #binary indicator for death in 1 yr
  death$d2013<-(death$dod<="2013-12-31") #binary indicator for death in 2 yrs
  death$d2016<-(death$dod<="2016-12-31") #binary indicator for death in 5 yrs
  death$d2018<-(death$dod<="2018-12-31") #binary indicator for death in 7 yrs
  patient2012 <- merge(patient2012,death[,c("patid","dod","dur","d2012","d2013","d2016","d2018")],by='patid',all.x=T)
  patient2012$d2012[is.na(patient2012$d2012)]<-0
  patient2012$d2013[is.na(patient2012$d2013)]<-0
  patient2012$d2016[is.na(patient2012$d2016)]<-0
  patient2012$d2018[is.na(patient2012$d2018)]<-0
  patient2012$dur[is.na(patient2012$dur)]<-76
  
# Note dur= NA means they are still alive-> we set them the largest observed duration+1!
  
  
  save(patient2012,file='C:\\MM_proc_data\\variables.RData',replace)
  save.image(file='C:\\MM_proc_data\\final.RData')

#*********************************************************************************#
# Plots and descriptive stats # 
  
# Descriptive tables
  
  colmorbid<-c("c_ckd","c_ap","c_ab","c_atr", "c_bli", "c_bro", "c_cld",
               "c_sin","c_copd","c_chd", "c_dem", "c_dia", "c_div", "c_hel", 
               "c_hf", "c_hyp", "c_ibd", "c_lea","c_ms", "c_prk", "c_pud","c_pvd", 
               "c_pro", "c_ops", "c_rhe", "c_str", "c_thy", "c_can", "c_anx", "c_ast",
               "c_con", "c_dep", "c_epi", "c_pnc", "c_ibs", "c_mig", "c_pso", "c_scz")
  patient2012$num_multimorb <- rowSums(patient2012[,..colmorbid],na.rm=TRUE) #.. is for data table to recognised columnlist index#
  hist(patient2012$num_multimorb)
  patient2012$binary_multimorb <- (patient2012$num_multimorb>=2)
  
  patient2012$agegroup <- cut(patient2012$age, breaks=c(-Inf, 24, 34, 44,54,64,74,84,Inf), labels=c("18-24","25-34","35-44","45-54", "55-64","65-74","75-84","85+"))
  
  # tabulation of sample size and percentages
  table(patient2012$gender)
  100*prop.table(table(patient2012$gender))
  
  table(patient2012$agegroup)
  100*prop.table(table(patient2012$agegroup))
  
  table(patient2012$imd2015_5)
  100*prop.table(table(patient2012$imd2015_5))
  
  patient2012 %>% group_by(gender) %>% summarise(sum(binary_multimorb)/n() * 100)
  patient2012 %>% group_by(agegroup) %>% summarise(sum(binary_multimorb)/n() * 100)
  patient2012 %>% group_by(imd2015_5) %>% summarise(sum(binary_multimorb)/n() * 100)
  
  # median number of morbidities by group
  summary(patient2012$num_multimorb,digits = 2) 
  
  group_by(patient2012, gender) %>%
    summarise(min=signif(min(num_multimorb),digits = 2),
              max=signif(max(num_multimorb),digits = 2),
              median = signif(median(num_multimorb),digits = 2)
              )
  
  group_by(patient2012, agegroup) %>%
    summarise(min=signif(min(num_multimorb),digits = 2),
              max=signif(max(num_multimorb),digits = 2),
              median = signif(median(num_multimorb),digits = 2)
              )
  group_by(patient2012, smokstatus) %>%
    summarise(min=signif(min(num_multimorb),digits = 2),
              max=signif(max(num_multimorb),digits = 2),
              median = signif(median(num_multimorb),digits = 2)
    )
  
  group_by(patient2012, imd2015_5) %>%
    summarise(min=signif(min(num_multimorb),digits = 2),
              max=signif(max(num_multimorb),digits = 2),
              median = signif(median(num_multimorb),digits = 2)
    )
  
  #Simple independence tests for relationship between multimorbidity (binary or number) and each demographic characteristics
  print(chisq.test(table(patient2012$gender, patient2012$binary_multimorb)))
  print(chisq.test(table(patient2012$agegroup, patient2012$binary_multimorb)))
  print(chisq.test(table(patient2012$imd2015_5, patient2012$binary_multimorb)))
  
  # prevalence of morbidities + index morbidities' associated co-morbidities
  tb_mltimorb<-subset(patient2012,select=c("c_ckd","c_ap","c_ab","c_atr", "c_bli", "c_bro", "c_cld",
                                           "c_sin","c_copd","c_chd", "c_dem", "c_dia", "c_div", "c_hel", 
                                           "c_hf", "c_hyp", "c_ibd", "c_lea","c_ms", "c_prk", "c_pud","c_pvd", 
                                           "c_pro", "c_ops", "c_rhe", "c_str", "c_thy", "c_can", "c_anx", "c_ast",
                                           "c_con", "c_dep", "c_epi", "c_pnc", "c_ibs", "c_mig", "c_pso", "c_scz"))
  colmorbidlist<-colnames(tb_mltimorb)
  
  colcount <- data.frame(colSums(tb_mltimorb == 1))
  colcount[2]<-colcount[1]/nrow(patient2012)*100
  names(colcount)<-c("Freq.","Prop.")
  prev<-colcount[order(-colcount$Prop.),]
  prev<- tibble::rownames_to_column(prev, "Morbidity")
  
  colmorbid<-prev[1:10,]$Morbidity # ordered rownames and work with only top 10
  for (i in colmorbid){
    print (i)
    tmp<-(tb_mltimorb[[i]]==1)*tb_mltimorb
    colcount <- data.frame(colSums(tmp == 1))
    colcount[2]<-colcount[1]/nrow(tb_mltimorb[tb_mltimorb[[i]]==1,])*100
    names(colcount)<-c("Freq.","Prop.")
    tmp$ttlcomorb<- rowSums(tmp==1,na.rm=TRUE)-1
    tmp$ttlcomorb<-ifelse(tmp$ttlcomorb<0,0,tmp$ttlcomorb)
    print(summary(tmp[tmp[[i]]==1]$ttlcomorb))
    comorb_ix_i<-colcount[order(-colcount$Prop.),]
    comorb_ix_i<- tibble::rownames_to_column(comorb_ix_i, "Morbidity")
    print(comorb_ix_i)
    rm(i,tmp,colcount)
  }
  
  
# Prevalence plots
  prev$Morbidityfull=c("Hypertension", "Painful condition", "Hearing loss","Depression","Irritable bowel syndrome",
                       "Asthma","Diabetes","Coronary heart disease","Thyroid disorders","Anxiety",
                       "Diverticular disease of intestine","Prostate disorder","Chronic kidney disease",
                       "Chronic sinusitis","Stroke & transient ischaemic attack","COPD","Atrial fibrillation",
                       "Rheumatoid arthritis","Constipation","Cancer","Peptic ulcer disease","Alcohol problems",
                       "Psychoactive substance misuse","Blindness and low vision","Psoriasis or eczema",
                       "Heart failure","Schizophrenia","Inflammatory bowel disease","Peripheral vascular disease",
                       "Dementia","Anorexia or bulimia","Chronic liver disease","Epilepsy","Migraine","Learning disability",
                       "Bronchiectasis","Multiple sclerosis","Parkinson's diseases")
  setwd("C:\\MM_proc_data\\Results")
  pdf(file="PrevalentDistr.pdf")
    p<-ggplot(prev, aes(x=reorder(Morbidityfull,Prop.), y=Prop.))+
      geom_bar(stat="identity",fill="dark grey")+
      coord_flip()+
      labs(title="Prevalence of morbidities (N= 391,669)",y="Proportions (%)",x="Morbidities")+
      geom_text(aes(label=sprintf("%0.2f", round(Prop., digits = 2))), vjust=0.2,size=2.5)+
      theme_minimal()
    p
  dev.off()
    
  
# Histograms of outcomes
  # Polypharmacy (prodcode) in 1 year
pdf(file="PolypharmaDistr1.pdf")
  p<-ggplot(patient2012[!ctprod12mons4plus==0], aes(x=ctprod12mons4plus)) + 
    geom_histogram(fill="dark grey",bins = 50)+
    labs(title="Distribution of non-zero treatment burden in 1 year (N= 138, 295)",
         subtitle="253,374/391,669 (64.7%) patients are with zero counts and included in the study",
         y="Frequency",x="Number of unique drugs prescribed at least 4 times in 1 year")
    theme_minimal()
    p
  dev.off()
  
  # Polypharmacy (BNF) in 1 year MORE ACCURATE!
  pdf(file="PolypharmaDistr2.pdf")
  p<-ggplot(patient2012[!ctprod12mons4plus_BNF==0], aes(x=ctprod12mons4plus_BNF)) + 
    geom_histogram(fill="dark grey",bins = 50)+
    labs(title="Distribution of non-zero treatment burden in 1 year (N= 142,721)",
         subtitle="248,948/391,669 (63.6%) patients are with zero counts and included in the study",
         y="Frequency",x="Number of drugs (in BNF codes) prescribed at least 4 times in 1 year")
  theme_minimal()
  p
  dev.off()
  
  # Number of GP conslutations in 1 yr
  pdf(file="GPconsultDistr.pdf")
  p<-ggplot(patient2012[!ctconsult12mons==0], aes(x=ctconsult12mons)) + 
    geom_histogram(fill="dark grey",bins = 50)+
    labs(title="Distribution of non-zero GP consultations in 1 year (N= 301,969)",
         subtitle="89,700/391,669 (22.9%) patients are with zero counts and included in the study",
         y="Frequency",x="Number of GP conslutations in 1 year")
  theme_minimal()
  p
  dev.off()
  
  # Median of inter-visit durations
  pdf(file="IntervisitDistr.pdf")
  p<-ggplot(patient2012[!interconsult12monsmed==0], aes(x=interconsult12monsmed)) + 
    geom_histogram(fill="dark grey",bins = 50)+
    labs(title="Distribution of non-zero median durations of inter-GP visits in 1 year (N= 273,845)",
         subtitle="117,824/391,669 (30.1%) patients are with zero counts and included in the study",
         y="Frequency",x="Median duration (in years) of inter-GP visits in 1 year")
  theme_minimal()
  p
  dev.off()
  
  # total number of hospitalisation spells
  pdf(file="HospspellDistr.pdf")
  p<-ggplot(patient2012[!cthospspl_12mons==0], aes(x=cthospspl_12mons)) + 
    geom_histogram(fill="dark grey",bins = 50)+
    labs(title="Distribution of non-zero hospitalisation spells in 1 year (N= 62,304)",
         subtitle="329,365/391,669 (84.1%) patients are with zero counts and included in the study",
         y="Frequency",x="Number of hospitalisation spells in 1 year")
    theme_minimal()
  p
  dev.off()
  
  # Distribution of time to death since 2012; cencoring time =2013end (2yr) and 2016end (5yr) (kaplan meier)
  pdf(file="TimetodeathDistr5yr.pdf")
  km_AG_fit5yr <- survfit(Surv(dur, d2016) ~ agegroup, data=patient2012)
  autoplot(km_AG_fit5yr, 
           main="Plot of Kaplan Meier survival curve for 5-year mortality by age group",
           xlab = "Time to death since 2012 (months)", ylab = "Proportion of individuals who have survived",
           censor.shape = '*',
           xlim = c(0,60)
          )
  dev.off()
  
  pdf(file="TimetodeathDistr2yr.pdf")
  km_AG_fit2yr <- survfit(Surv(dur, d2013) ~ agegroup, data=patient2012)
  autoplot(km_AG_fit2yr, 
           main="Plot of Kaplan Meier survival curve for 2-year mortality by age group",
           xlab = "Time to death since 2012 (months)", ylab = "Proportion of individuals who have survived",
           censor.shape = '*',
           xlim = c(0,24)
            )
  dev.off()

  # A set of survival plots for 5-year mortality by "multimorbidity or not" 
  pdf(file="TimetodeathDistr5yrmm.pdf")
  km3_fit5yrmm<- survfit(Surv(dur, d2016) ~ binary_multimorb, data=patient2012[patient2012$agegroup=="35-44",])
  km5_fit5yrmm<- survfit(Surv(dur, d2016) ~ binary_multimorb, data=patient2012[patient2012$agegroup=="55-64",])
  km7_fit5yrmm<- survfit(Surv(dur, d2016) ~ binary_multimorb, data=patient2012[patient2012$agegroup=="75-84",])
  
  summary(km3_fit5yrmm)
  a<-autoplot(km3_fit5yrmm, 
           main="Age group: 35-44",
           legend.title="Morbidities",
           xlab = "Time to death since 2012 (months)", ylab = "Proportion of individuals who have survived",
           censor.shape = '*',
           xlim = c(0,60)
  )
  
  b<-autoplot(km5_fit5yrmm, 
           main="Age group: 55-64",
           xlab = "Time to death since 2012 (months)", ylab = "Proportion of individuals who have survived",
           censor.shape = '*',
           xlim = c(0,60)
           
  )
  c<-autoplot(km7_fit5yrmm, 
              main="Age group: 75-84",
              xlab = "Time to death since 2012 (months)", ylab = "Proportion of individuals who have survived",
              censor.shape = '*',
              xlim = c(0,60)
  )
  grid.arrange(a, b,c, ncol=3,
               top = textGrob("Plot of Kaplan Meier survival curve for 5-year mortality by morbidities"))
  dev.off()
  
  
  #A set of survival plots for 2-year mortality by "multimorbidity or not" 
  pdf(file="TimetodeathDistr2yrmm.pdf")
  km3_fit2yrmm<- survfit(Surv(dur, d2013) ~ binary_multimorb, data=patient2012[patient2012$agegroup=="35-44",])
  km5_fit2yrmm<- survfit(Surv(dur, d2013) ~ binary_multimorb, data=patient2012[patient2012$agegroup=="55-64",])
  km7_fit2yrmm<- survfit(Surv(dur, d2013) ~ binary_multimorb, data=patient2012[patient2012$agegroup=="75-84",])
  
  
  a<-autoplot(km3_fit2yrmm, 
              main="Age group: 35-44",
              xlab = "Months", ylab = "Proportion of individuals who have survived",
              censor.shape = '*',
              xlim = c(0,24)
  )
  
  b<-autoplot(km5_fit2yrmm, 
              main="Age group: 55-64",
              xlab = "Months",ylab="",
              censor.shape = '*',
              xlim = c(0,24)
  )
  c<-autoplot(km7_fit2yrmm, 
              main="Age group: 75-84",
              xlab = "Months",ylab="",
              censor.shape = '*',
              xlim = c(0,24)
  )
  grid.arrange(a, b, c,ncol=3,
               top = textGrob("Plot of Kaplan Meier survival curve for 2-year mortality by morbidities"))
  dev.off()
  
  #rm(km_5564_fit5yrmm)
  save(patient2012,file='C:\\MM_proc_data\\variables.RData',replace)
  save.image(file='C:\\MM_proc_data\\final.RData')
  
#*********************************************************************************#
#  Heat maps for association between morbidities among individuals who have at least two multimorbidities
  cormat1 <- round(cor(patient2012[binary_multimorb==1& agegroup=="35-44",colmorbidlist,with=FALSE]),2)
  cormat2 <- round(cor(patient2012[binary_multimorb==1& agegroup=="55-64",colmorbidlist,with=FALSE]),2)
  cormat3 <- round(cor(patient2012[binary_multimorb==1& agegroup=="75-84",colmorbidlist,with=FALSE]),2)
  
  # Reorder the correlation matrix
  cormat1 <- reorder_cormat(cormat1)
  cormat2 <- reorder_cormat(cormat2)
  cormat3 <- reorder_cormat(cormat3)
  
  # Melt the correlation matrix
  melted_cormat1 <- melt(cormat1, na.rm = TRUE)
  melted_cormat2 <- melt(cormat2, na.rm = TRUE)
  melted_cormat3 <- melt(cormat3, na.rm = TRUE)
  
  
  table(melted_cormat1$value)
  table(melted_cormat2$value)
  table(melted_cormat3$value)
  
  
  # Create a ggheatmap
  plot1<-ggplot(melted_cormat1, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(min(melted_cormat1$value),max(melted_cormat1[!melted_cormat1$value==1,]$value)), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                     size = 8, hjust = 1))+
    labs(title="Correlation between morbidities (age group: 35-44)",
         y="Morbidities",x="Morbidities")+
    coord_fixed()
  plot2<-ggplot(melted_cormat2, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(min(melted_cormat2$value),max(melted_cormat2[!melted_cormat2$value==1,]$value)), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                     size = 8, hjust = 1))+
    labs(title="Correlation between morbidities (age group: 55-64)",
         y="Morbidities",x="Morbidities")+
    coord_fixed()
  plot3<-ggplot(melted_cormat3, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(min(melted_cormat3$value),max(melted_cormat3[!melted_cormat3$value==1,]$value)), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                     size = 8, hjust = 1))+
    labs(title="Correlation between morbidities (age group: 75-84)",
         y="Morbidities",x="Morbidities")+
    coord_fixed()
  dev.off()
  #plot jointly
 pdf(file="Heatmaps=35-44.pdf")
 plot1
 dev.off()
 pdf(file="Heatmaps=55-64.pdf")
 plot2
 dev.off()
 pdf(file="Heatmaps=75-84.pdf")
 plot3
 dev.off()
  
  rm(mat1,mat2, mat3, mat4, mat5, mat6, mat7, mat8, ht1,ht2)
  save(patient2012,file='C:\\MM_proc_data\\variables.RData',replace)
  save.image(file='C:\\MM_proc_data\\final.RData')
