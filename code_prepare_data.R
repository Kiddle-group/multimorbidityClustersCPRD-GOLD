# Last updated: 29 May 2019 #
# Script to load and merge relevant CPRD data
# This file accompanies the paper "Identification of clusters of multimorbid patients in UK general practice: a latent class analysis."
# This file contains codes to load massive CPRD patient data bits by bits into RAM using ffdf.

# Load required packages and helper functions

# install.packages('readstata13') # allows reading in of stata dta files
# install.packages('rms') # scoring functions and more
# install.packages('survival') # survival 
# install.packages('ROCR') # ROC curves and prediction metrics
# install.packages('xlsx') # read from excel files
# install.packages('ff')
# install.packages('ffbase')
# install.packages('data.table')
# install.packages('bit64')


library(readstata13) # allows reading in of stata dta files
library(rms) # scoring functions and more
library(survival) # survival 
library(ROCR) # ROC curves and prediction metrics
library(xlsx) # read from excel files
library(ff)
library(ffbase)
library(data.table)
library(bit64)
source('C:\\MMcluster\\code_misc.R') # load helper functions


# if ffdf aren't made yet
#source('code_make_ffdf.R')

setwd("C:\\MM_proc_data\\ffdb")

# if ffdf are made but need to be loaded
load.ffdf('C:\\MM_proc_data\\ffdb\\clinical')
load.ffdf('C:\\MM_proc_data\\ffdb\\therapy')
load.ffdf('C:\\MM_proc_data\\ffdb\\test')
load.ffdf('C:\\MM_proc_data\\ffdb\\patient')
load.ffdf('C:\\MM_proc_data\\ffdb\\immunisation')
load.ffdf('C:\\MM_proc_data\\ffdb\\referral')
load.ffdf('C:\\MM_proc_data\\ffdb\\additional')
load.ffdf('C:\\MM_proc_data\\ffdb\\consultation')


load.ffdf('C:\\MM_proc_data\\ffdb\\clinical2012_before')
load.ffdf('C:\\MM_proc_data\\ffdb\\therapy2012_before')
load.ffdf('C:\\MM_proc_data\\ffdb\\test2012_before')
load.ffdf('C:\\MM_proc_data\\ffdb\\immunisation2012_before')
load.ffdf('C:\\MM_proc_data\\ffdb\\referral2012_before')
load.ffdf('C:\\MM_proc_data\\ffdb\\consultation2012_after')


load.ffdf('C:\\MM_proc_data\\ffdb\\therapy2012_after')


# NOTE: if the previous load.ffdfs are loaded, the following codes shall not be repeated
clinical2012_before <- clinical[clinical$eventdate <= '2012-01-01' ,]
save.ffdf(clinical2012_before, dir = 'C:\\MM_proc_data\\ffdb\\clinical2012_before')
head(clinical2012_before)

test2012_before <- test[test$eventdate <= '2012-01-01' ,]
save.ffdf(test2012_before, dir = 'C:\\MM_proc_data\\ffdb\\test2012_before')
head(test2012_before)

referral2012_before <- referral[referral$eventdate <= '2012-01-01' ,]
save.ffdf(referral2012_before, dir = 'C:\\MM_proc_data\\ffdb\\referral2012_before')
head(referral2012_before)

immunisation2012_before <- immunisation[immunisation$eventdate <= '2012-01-01' ,]
save.ffdf(immunisation2012_before, dir = 'C:\\MM_proc_data\\ffdb\\immunisation2012_before')
head(immunisation2012_before)

therapy2012_before <- therapy[therapy$eventdate <= '2012-01-01' ,]
save.ffdf(therapy2012_before, dir = 'C:\\MM_proc_data\\ffdb\\therapy2012_before')
head(therapy2012_before)

#***************************************************************************************#

#For outcomes

therapy2012_after <- therapy[therapy$eventdate >='2012-01-01' ,]
save.ffdf(therapy2012_after, dir = 'C:\\MM_proc_data\\ffdb\\therapy2012_after')
head(therapy2012_after)

#first load consultation file!
consultation2012_after <- consultation[consultation$eventdate >='2012-01-01' ,]
save.ffdf(consultation2012_after, dir = 'C:\\MM_proc_data\\ffdb\\consultation2012_after')
head(consultation2012_after)


