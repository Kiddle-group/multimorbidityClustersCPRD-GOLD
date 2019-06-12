# Last updated: 29 May 2019 #
# This file accompanies the paper "Identification of clusters of multimorbid patients in UK general practice: a latent class analysis."
# This file contains codes to extract LTCs, patient demographics and outcomes from the CPRD data #

setwd('Q:\\kiddlegroup\\MM clustering project\\unzipped raw')

patient <- read.table.ffdf(file = 'comorb_Extract_Patient_001.txt',sep='\t',header=T)
patient$crd <- as.Date(patient$crd,format='%d/%m/%Y')
practice <- read.table.ffdf(file = 'comorb_Extract_Practice_001.txt',sep='\t',header=T)
practice$uts <- as.Date(practice$uts,format='%d/%m/%Y')

str <- as.character(as.ram(patient$patid))

patient$pracid <- as.ff(as.numeric(substr(str,nchar(str)-2,nchar(str))))
patient <- merge(patient,practice,all.x=T,by='pracid')
save.ffdf(patient,dir = 'C:\\MM_proc_data\\ffdb\\patient')

additional <- read.table.ffdf(file = 'comorb_Extract_Additional_001.txt',sep='\t',header=T)
save.ffdf(additional,dir = 'C:\\MM_proc_data\\ffdb\\additional')

referral <- read.table.ffdf(file = 'comorb_Extract_Referral_001.txt',sep='\t',header=T)
referral$eventdate <- as.Date(referral$eventdate,format='%d/%m/%Y')

save.ffdf(referral,dir = 'C:\\MM_proc_data\\ffdb\\referral')

immunisation <- read.table.ffdf(file = 'comorb_Extract_Immunisation_001.txt',sep='\t',header=T)
immunisation$eventdate <- as.Date(immunisation$eventdate,format='%d/%m/%Y')
save.ffdf(immunisation,dir = 'C:\\MM_proc_data\\ffdb\\immunisation')

f <- list.files()

files <- grep(pattern = 'Clinical',f,value=T)

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
tmp <- fread(files[1],header=T,sep='\t',stringsAsFactors = T)
tmp$eventdate <- as.Date(tmp$eventdate,format='%d/%m/%Y')
tmp$sysdate <- as.Date(tmp$sysdate,format='%d/%m/%Y')

clinical <- as.ffdf(tmp)
print(dim(clinical))

for (i in 2:length(files)){
  
  print(i)
  
  tmp <- fread(files[i],header=T,sep='\t')
  tmp$eventdate <- as.Date(tmp$eventdate,format='%d/%m/%Y',stringsAsFactors=T)
  tmp$sysdate <- as.Date(tmp$sysdate,format='%d/%m/%Y',stringsAsFactors=T)
  
  clinical <- ffdfappend(clinical,tmp)
  
  print(dim(clinical))
}

save.ffdf(clinical,dir = 'C:\\MM_proc_data\\ffdb\\clinical')


files <- grep(pattern = 'Test',f,value=T)

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
tmp <- fread(files[1],header=T,sep='\t')
tmp$eventdate <- as.Date(tmp$eventdate,format='%d/%m/%Y')
tmp$sysdate <- as.Date(tmp$sysdate,format='%d/%m/%Y')
tmp <- subset(tmp, select = -c(textid,staffid,data8))
test <- as.ffdf(tmp)

for (i in 2:length(files)){
  
  print(i)
  
  tmp <- fread(files[i],header=T,sep='\t')
  tmp$eventdate <- as.Date(tmp$eventdate,format='%d/%m/%Y',stringsAsFactors=T)
  tmp$sysdate <- as.Date(tmp$sysdate,format='%d/%m/%Y',stringsAsFactors=T)
  tmp <- subset(tmp, select = -c(textid,staffid,data8))
  test <- ffdfappend(test,tmp)
  
}

save.ffdf(test,dir = 'C:\\MM_proc_data\\ffdb\\test')

files <- grep(pattern = 'Therapy',f,value=T)

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
tmp <- fread(files[1],header=T,sep='\t')
tmp$eventdate <- as.Date(tmp$eventdate,format='%d/%m/%Y')
tmp$sysdate <- as.Date(tmp$sysdate,format='%d/%m/%Y')
tmp <- subset(tmp, select = -c(textid,staffid))

therapy <- as.ffdf(tmp)

for (i in 2:length(files)){
  
  print(i)
  
  tmp <- fread(files[i],header=T,sep='\t')
  tmp$eventdate <- as.Date(tmp$eventdate,format='%d/%m/%Y',stringsAsFactors=T)
  tmp$sysdate <- as.Date(tmp$sysdate,format='%d/%m/%Y',stringsAsFactors=T)
  tmp <- subset(tmp, select = -c(textid,staffid))
  therapy <- ffdfappend(therapy,tmp)
  
}

# add columns for productcodes
prod_browse <- fread('Q:\\kiddlegroup\\code browser data\\product.txt',skip=2,header = F,stringsAsFactors = T)
colnames(prod_browse)[1] <- 'prodcode'
colnames(prod_browse)[5] <- 'drugname'
prod_browse <- as.ffdf(prod_browse)
therapy <- merge(therapy,prod_browse,by='prodcode',all.x=T)
colnames(therapy)[21] <- 'BNF'
save.ffdf(therapy,dir = 'C:\\MM_proc_data\\ffdb\\therapy',overwrite=TRUE)

setwd('Q:\\kiddlegroup\\MM clustering project\\unzipped raw')

files <- grep(pattern = 'Consultation',f,value=T)

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
tmp <- fread(files[1],header=T,sep='\t')
tmp$eventdate <- as.Date(tmp$eventdate,format='%d/%m/%Y')
tmp$sysdate <- as.Date(tmp$sysdate,format='%d/%m/%Y')
tmp <- subset(tmp, select = -c(staffid))

consultation <- as.ffdf(tmp)

for (i in 2:length(files)){
  
  print(i)
  
  tmp <- fread(files[i],header=T,sep='\t')
  tmp$eventdate <- as.Date(tmp$eventdate,format='%d/%m/%Y',stringsAsFactors=T)
  tmp$sysdate <- as.Date(tmp$sysdate,format='%d/%m/%Y',stringsAsFactors=T)
  tmp <- subset(tmp, select = -c(staffid))
  consultation<- ffdfappend(consultation,tmp)
  
}

save.ffdf(consultation,dir = 'C:\\MM_proc_data\\ffdb\\consultation',overwrite=TRUE)