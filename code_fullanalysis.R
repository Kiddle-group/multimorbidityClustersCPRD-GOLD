## set working dir ##
setwd("C:\\MM_proc_data")

# Last updated: 11 June 2019 #
# Note: 
# Final dataset = file='C:\\MM_proc_data\\final_patid_strata_final_filterdur.RData'
# This code accompanies the paper "Identification of clusters of multimorbid patients in UK general practice: a latent class analysis."
# This files contains codes for the entire analysis (dependencies: load_variables.r, make_ffdf.r, misc.r, prepare_data.r, mplus.r)

## LIBRARIES ##
###################################################################################

source('C:\\MMcluster\\code_misc.R') # load helper functions
#install.packages("CALIBERdatamanage", repos="http://R-Forge.R-project.org")
library(CALIBERdatamanage)
library(ggplot2)        # for generating visualizations
#install.packages("plyr")
library(plyr)      # for data manipulation
library(dplyr)      # for data manipulation
#install.packages("AMR")
library(AMR)      # for data manipulation
library(data.table)
#install.packages("gridGraphics")
library(gridGraphics)
source("http://pcwww.liv.ac.uk/~william/R/crosstab.r")
library(reshape2)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap", version = "3.8")
library(ComplexHeatmap)
library(gridExtra)
#install.packages('survminer')
source("https://bioconductor.org/biocLite.R")
library(survminer)
library(survival)
#install.packages('ranger')
library(ranger)
#install.packages(c("ggplot2", "ggfortify"))
library(ggfortify)
library(grid)
#install.packages("lcmm")
library(lcmm)
#install.packages("randomForest")
library(randomForest)
#install.packages("poLCA")
library("poLCA")
#install.packages(c("ggparallel","igraph","tidyr","knitr"))
library("ggparallel")
library("igraph")
library("tidyr")
library("knitr")
#install.packages("ztable")
library("ztable")
#install.packages("forcats")
library("forcats")
#install.packages("gapminder")
#install.packages("ps")
#library(ps,lib.loc="C:/Users/yajing/Documents/R/win-library/3.5")
library("gapminder")
#install.packages("tidyverse")
library("tidyverse")
#install.packages(c('devtools','processx'))
library(devtools)
#install.packages("jtools")
library(jtools)
#install_github("MattisvdBergh/LCT")
library(LCTpackage)
#install.packages(c("cluster.datasets"), dependencies = TRUE)
library(cluster.datasets)
#install.packages("CDM")
library(CDM)
#install.packages("ggstance")
library(ggstance)
#install_github("wbonat/mcglm")
#install.packages("mcglm")
#install.packages("fastDummies")
library(fastDummies)
#devtools::install_github("ewenharrison/finalfit")
#library(finalfit)
#install.packages("memisc")
library(memisc)
library(Hmisc)
#biocLite("survcomp")
library(survcomp)
#install.packages("DHARMa")
library(DHARMa)
#install.packages("countreg", repos="http://R-Forge.R-project.org")
library(countreg)
#install.packages("depmixS4")
library(depmixS4)
#install.packages("rio")
library(rio)
#install.packages("MplusAutomation")
library(MplusAutomation)
#install.packages("scales")
#install.packages("cowplot")
library(scales)
library(cowplot)
#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
library("rhdf5")
#install.packages("devEMF")
library("devEMF")
#install.packages("pheatmap")
library(pheatmap)
#install.packages("sjPlot")
#install.packages("sjlabelled")
#install.packages("sjmisc")

library(sjPlot)
library(sjlabelled)
library(sjmisc)
#install.packages("devEMF")
#install.packages("pglm")
library(pglm)
library(devEMF)
source('C:\\MMcluster\\mplus.R') # load helper functions
library(caret)
#install.packages("qwraps2")
library(qwraps2)
#install.packages("htmlTable")
library("htmlTable")

#install.packages(c("RCurl","XML"))
library(RCurl)
library(XML)
#install.packages("pdist")
library("pdist")

#########################################################################

# LOAD data #
load(file="final_patid_nonstratified.RData")
rm(list=setdiff(ls(), c("patient2012","patient2012_copy","patient2012_copy_mm","prev","colmorbidlist","tb_mltimorb","morbidptn")))

# LOAD saved data #
save.image(file='final_patid_strata.RData')                                                                                                                                                                    


# Descriptive stats for the whole population #
a<-patient2012
summary(a$num_multimorb)
summary(a[a$agegroup=="85+",]$num_multimorb)
summary(a[a$imd2015_5=="1",]$num_multimorb)

## The following analysis works on training sets only: split at patient level ##
## Note: work only on indivdiuals with 2+ morbidities #
## Stratified clustering ##

#############################################################################################
# Latent class models (finite mixture models) #

# Note:  binary 0/1 needs to be 1/2. NA values need to be in a separate category  
# Setting nrep>1 automates the search for the global-rather than just a local-maximum of the log-likelihood function. 
#     poLCA returns the parameter estimates corresponding to the model with the greatest log-likelihood.
# poLCA synatx is given below that can replace Mplus syntax. However, we note that poLCA requires a lot more iterations
#	to obtain the best likelihood as Mplus did (2000 replications and more).
# We have run both poLCA and Mplus modules to confirm final results.  
##############################################################################################

# list of extracted LTCs (have been extracted using load_variables.R) #
f <-cbind(c_ckd,c_ap,c_ab,c_atr,c_bli,c_bro,c_cld,c_sin,c_copd,c_chd,c_dem,c_dia,c_div,c_hel,c_hf,c_hyp,c_ibd,c_lea,c_ms,
          c_prk,c_pud,c_pvd,c_pro,c_ops,c_rhe,c_str,c_thy,c_can,c_anx,c_ast,c_con,c_dep,c_epi,c_pnc,c_ibs,c_mig,c_pso,c_scz)~1

# generate prevalence table by age group: 18-44, 45-64, 65-84, 85+ #
patient2012$strataid<- cut(patient2012$age, breaks=c(-Inf, 44,64,84,Inf), labels=c("18-44","45-64","65-84","85+"))

# generate stratified data #
tb_mltimorb_1844<-subset(patient2012,strataid=="18-44",select=c("c_ckd","c_ap","c_ab","c_atr", "c_bli", "c_bro", "c_cld",
                                         "c_sin","c_copd","c_chd", "c_dem", "c_dia", "c_div", "c_hel", 
                                         "c_hf", "c_hyp", "c_ibd", "c_lea","c_ms", "c_prk", "c_pud","c_pvd", 
                                         "c_pro", "c_ops", "c_rhe", "c_str", "c_thy", "c_can", "c_anx", "c_ast",
                                         "c_con", "c_dep", "c_epi", "c_pnc", "c_ibs", "c_mig", "c_pso", "c_scz"))

tb_mltimorb_4564<-subset(patient2012,strataid=="45-64",select=c("c_ckd","c_ap","c_ab","c_atr", "c_bli", "c_bro", "c_cld",
                                         "c_sin","c_copd","c_chd", "c_dem", "c_dia", "c_div", "c_hel", 
                                         "c_hf", "c_hyp", "c_ibd", "c_lea","c_ms", "c_prk", "c_pud","c_pvd", 
                                         "c_pro", "c_ops", "c_rhe", "c_str", "c_thy", "c_can", "c_anx", "c_ast",
                                         "c_con", "c_dep", "c_epi", "c_pnc", "c_ibs", "c_mig", "c_pso", "c_scz"))

tb_mltimorb_6584<-subset(patient2012,strataid=="65-84",select=c("c_ckd","c_ap","c_ab","c_atr", "c_bli", "c_bro", "c_cld",
                                         "c_sin","c_copd","c_chd", "c_dem", "c_dia", "c_div", "c_hel", 
                                         "c_hf", "c_hyp", "c_ibd", "c_lea","c_ms", "c_prk", "c_pud","c_pvd", 
                                         "c_pro", "c_ops", "c_rhe", "c_str", "c_thy", "c_can", "c_anx", "c_ast",
                                         "c_con", "c_dep", "c_epi", "c_pnc", "c_ibs", "c_mig", "c_pso", "c_scz"))

tb_mltimorb_85plus<-subset(patient2012,strataid=="85+",select=c("c_ckd","c_ap","c_ab","c_atr", "c_bli", "c_bro", "c_cld",
                                         "c_sin","c_copd","c_chd", "c_dem", "c_dia", "c_div", "c_hel", 
                                         "c_hf", "c_hyp", "c_ibd", "c_lea","c_ms", "c_prk", "c_pud","c_pvd", 
                                         "c_pro", "c_ops", "c_rhe", "c_str", "c_thy", "c_can", "c_anx", "c_ast",
                                         "c_con", "c_dep", "c_epi", "c_pnc", "c_ibs", "c_mig", "c_pso", "c_scz"))
# initial checks for prevalence (use 85plus as an example, repeated for all age strata) #
colmorbidlist<-colnames(tb_mltimorb_85plus)

colcount <- data.frame(colSums(tb_mltimorb_85plus == 1))
colcount[2]<-colcount[1]/nrow(tb_mltimorb_85plus)*100
names(colcount)<-c("Freq.","Prop.")
prev_85plus<-colcount[order(-colcount$Prop.),]
prev_85plus<- tibble::rownames_to_column(prev_85plus, "Morbidity")

# Focus on the multimorbid patients only (2+ LTCs) #
# Split at patient level for EACH age strata #

# generate strata indicator in the multimorbid population #
patient2012_copy_mm$dur<-as.numeric(patient2012_copy_mm$dur)
patient2012_copy_mm$strataid<- cut(patient2012_copy_mm$age, breaks=c(-Inf, 44,64,84,Inf), labels=c("18-44","45-64","65-84","85+"))

split_pat<-0.8 # 80% training data

# stratified random sampling: random sampling within strata #
set.seed(102)
train.index <- createDataPartition(patient2012_copy_mm$strataid, p = split_pat, list = FALSE)
train_pat_mm <- patient2012_copy_mm[train.index,]
test_pat_mm  <- patient2012_copy_mm[-train.index,]

# Check raw distribution of varibles in whole multimorbid datasets across 4 age strata #
a<-patient2012_copy_mm
options(qwraps2_markup="latex",digits=1)
summary_statistics <-
  list(
     "Patients" = 
       list(
        "Patients" = ~length(patid)
           ),
    "Gender" = 
      list(
        "Female" = ~ qwraps2::n_perc0(gender == 2,na_rm = TRUE),
        "Missing" = ~sum(is.na(gender))
          ),
    "# morbidities" =
      list(
        "Median (Q1, Q3)" = ~qwraps2::median_iqr(num_multimorb, na_rm = TRUE,show_n="never"),
        "Missing" = ~sum(is.na(num_multimorb))
      ),
    "Age" =
      list(
        "Median (Q1, Q3)" = ~qwraps2::median_iqr(age, na_rm = TRUE,show_n="never"),
        "Missing" = ~sum(is.na(age))
      ),
    "BMI" =
      list(
        "Median (Q1, Q3)" = ~qwraps2::median_iqr(bmi, na_rm = TRUE,show_n="never"),
        "Missing" = ~sum(is.na(bmi))
      ),
    "Smoking status" =
      list(
        "Current smoker" = ~ qwraps2::n_perc0(smokstatus == 1,na_rm = TRUE),
        "Never smoker" = ~ qwraps2::n_perc0(smokstatus == 2,na_rm = TRUE),
        "Ex smoker" = ~ qwraps2::n_perc0(smokstatus == 3,na_rm = TRUE),
        "Missing" = ~sum(is.na(smokstatus))
      ),
    "Index of multiple deprivation in quintiles" =
      list(
        "1 (least depreviation)" = ~ qwraps2::n_perc0(imd2015_5 == 1,na_rm = TRUE),
        "2 " = ~ qwraps2::n_perc0(imd2015_5 == 2,na_rm = TRUE),
        "3 " = ~ qwraps2::n_perc0(imd2015_5 == 3,na_rm = TRUE),
        "4 " = ~ qwraps2::n_perc0(imd2015_5 == 4,na_rm = TRUE),
        "5 (greatest depreviation)" = ~ qwraps2::n_perc0(imd2015_5 == 5,na_rm = TRUE),
        "Missing" = ~sum(is.na(imd2015_5))
      )
  )

a_summary <-qwraps2::summary_table(dplyr::group_by(a, strataid), summary_statistics)
a_summary <-qwraps2::summary_table(dplyr::group_by(a), summary_statistics)


# added demographic summaries for the whole population
# a<-patient2012 %>%
#   group_by(smokstatus) %>%
#   mutate(median= quantile(num_multimorb,0.5),
#          Q1=quantile(num_multimorb,0.25),
#          Q3=quantile(num_multimorb,0.75))
# a <- a %>%
#   group_by(smokstatus) %>%
#   summarise(median = mean(median), Q1= mean(Q1), Q3 = mean(Q3),
#             num=n(),
#             num_prop=n()/nrow(.)*100,
#             multimorb=sum(binary_multimorb)/nrow(.)*100)


# output summary stats to tables

capture.output(a_summary, file = "C:\\MM_proc_data\\Results\\summarystats_mm_strata.txt")

# check distribution #
prop.table(table(patient2012_copy_mm$strataid)) #0.1351989 0.3188471 0.4371837 0.1087703 
prop.table(table(train_pat_mm$strataid)) # 0.1351978 0.3188438 0.4371819 0.1087765 
prop.table(table(test_pat_mm$strataid)) # 0.1352032 0.3188604 0.4371908 0.1087456

# check total number of patterns of morbidity that involve each disease #
mapply(train_pat_mm[,c(colmorbidlist)], FUN=sum, 1)
mapply(test_pat_mm[,c(colmorbidlist)], FUN=sum, 1)

# Check all morbidity patterns in the dataset
fpattern <- paste("~", paste(colmorbidlist, collapse="+"))

morbidptn_train_mm<-ddply(train_pat_mm[,c(colmorbidlist)],as.formula(fpattern),summarise, freq=length(c_ckd))
morbidptn_train_mm$ttlptn<-nrow(morbidptn_train_mm)

morbidptn_test_mm<-ddply(test_pat_mm[,c(colmorbidlist)],as.formula(fpattern),summarise, freq=length(c_ckd))
morbidptn_test_mm$ttlptn<-nrow(morbidptn_test_mm)

# Make word-clear emf plots #
setwd("C:\\MM_proc_data\\Results")
emf(file="PrevalentDistr.emf",emfPlus = F)
p<-ggplot(prev, aes(x=reorder(Morbidityfull,Prop.), y=Prop.))+
  geom_bar(stat="identity",fill="dark grey")+
  coord_flip()+
  labs(title="Prevalence of morbidities (N= 391,669)",y="Proportions (%)",x="Morbidities")+
  geom_text(aes(label=sprintf("%0.2f", round(Prop., digits = 2))), vjust=0.2,size=2.5)+
  theme_minimal()
p
dev.off()


######################################################################################## 

# Prepare data for poLCA / Mplus#
# generate sub dataframes for analysis (ease loading times) #

train_pat_mm_1844<-train_pat_mm[train_pat_mm$strataid=="18-44",]
train_pat_mm_4564<-train_pat_mm[train_pat_mm$strataid=="45-64",]  
train_pat_mm_6584<-train_pat_mm[train_pat_mm$strataid=="65-84",]
train_pat_mm_85plus<-train_pat_mm[train_pat_mm$strataid=="85+",]    

test_pat_mm_1844<-test_pat_mm[test_pat_mm$strataid=="18-44",]
test_pat_mm_4564<-test_pat_mm[test_pat_mm$strataid=="45-64",]  
test_pat_mm_6584<-test_pat_mm[test_pat_mm$strataid=="65-84",]
test_pat_mm_85plus<-test_pat_mm[test_pat_mm$strataid=="85+",]  


dflist<-list(train_pat_mm_1844,train_pat_mm_4564,train_pat_mm_6584,train_pat_mm_85plus,   
             test_pat_mm_1844,test_pat_mm_4564,test_pat_mm_6584,test_pat_mm_85plus)
name.dflist<-c("train_pat_mm_1844","train_pat_mm_4564","train_pat_mm_6584","train_pat_mm_85plus",   
               "test_pat_mm_1844","test_pat_mm_4564","test_pat_mm_6584","test_pat_mm_85plus")

# check number of patterns in the training for each age strata #
for (i in 1:length(dflist)){
  print(i)
  fpattern <- paste("~", paste(colmorbidlist, collapse="+"))
  myname<-ddply(dflist[[i]][,c(colmorbidlist)],as.formula(fpattern),summarise, freq=length(c_ckd))
  myname$ttlptn<-nrow(myname)
  assign(paste0("morbidptn_",name.dflist[i]), myname)
}

for (i in 1:length(dflist)){
  myfile <- file.path("H:\\Mplus_MM\\unsup_patidstrata", paste0(name.dflist[i],".dat"))
  prepareMplusData(as.data.frame(dflist[[i]]), myfile, 
                 keepCols=c("c_ckd", "c_ap", "c_ab", "c_atr", "c_bli", "c_bro",
                            "c_cld", "c_sin", "c_copd", "c_chd", "c_dem", "c_dia",
                            "c_div","c_hel", "c_hf", "c_hyp", "c_ibd", "c_lea",
                            "c_ms", "c_prk","c_pud", "c_pvd", "c_pro", "c_ops", 
                            "c_rhe", "c_str", "c_thy", "c_can", "c_anx", "c_ast",
                            "c_con", "c_dep", "c_epi", "c_pnc", "c_ibs", "c_mig",
                            "c_pso", "c_scz","ctconsult12mons", "cthospspl_12mons", 
                            "ctprod12mons4plus_BNF", "agegroup", "gender",
                            "imd2015_5","dur","d2012","d2013","d2016","d2018","patid"))
}

# Run Mplus (alternaive poLCA) in high performance computing with 20 cores #
# We include in the main scripts final Mplus codes due to its better likelihood solutions. For poLCA users:
#################################################################################################
# Step 1: run LCA on each of the age strata
#	ncls<-20 # number of models to run
#
#	set.seed(54321)
#	for(i in 1:ncls){ # note class=1 -> log-linear independence model
#  	print(i)
#  	temp<- poLCA(f, train_prac_mm, nclass=i, maxiter=5000, 
#               tol=1e-5, na.rm=FALSE,  
#               nrep=20000, verbose=TRUE, calc.se=TRUE) 
#  	assign(paste0("lc_mm", i), temp)
#	}  
# Step 2: compute a set of model-selection statistics and choose the best model
#	results_mm <- data.frame(Modell=1,
#                         log_likelihood=lc_mm1$llik,
#                         df = lc_mm1$resid.df,
#                         BIC=lc_mm1$bic,
#                         SABIC=  (-2*lc_mm1$llik) + (log((lc_mm1$N + 2)/24)) * (lc_mm1$bic-(-2*lc_mm1$llik))/log(lc_mm1$N),
#                         likelihood_ratio=lc_mm1$Gsq)
# 	results_mm$R2_entropy
#	results_mm$Modell<-as.integer(results_mm$Modell)
#	results_mm[1,7]<-c("-")
#
# 	# Build lists of latent class objects
#	my.lcalist <- lapply(paste('lc_mm', seq(1,cls,1), sep=''), get)

# 	# For all results, i in 1:20; for restrictions on class size, i in 1:10
#	for(i in 2:10){
#  	results_mm[i,1]<-i
# 	results_mm[i,2]<-my.lcalist[[i]]$llik
#  	results_mm[i,3]<-my.lcalist[[i]]$resid.df
#  	results_mm[i,4]<-my.lcalist[[i]]$bic
#  	results_mm[i,5]<-(-2*my.lcalist[[i]]$llik) + (log((my.lcalist[[i]]$N + 2)/24)) * (my.lcalist[[i]]$bic-(-2*my.lcalist[[i]]$llik))/log(my.lcalist[[i]]$N) #SAbic
#  	results_mm[i,6]<-my.lcalist[[i]]$Gsq  
#  	error_prior<-entropy(my.lcalist[[i]]$P) # class proportions model i
#  	error_post<-mean(apply(my.lcalist[[i]]$posterior,1, entropy),na.rm = TRUE)
#  	results_mm[i,7]<-round(((error_prior-error_post) / error_prior),3)
#	}
#
# 	# Combining results to a dataframe
#	colnames(results_mm)<-c("Model","log-likelihood","resid. df","BIC","SABIC","likelihood-ratio","Entropy")
#	lca_results_mm<-results_mm
#	ztable::ztable(lca_results_mm)
# 	bestmodel_mm<-lc_mm4 # for example
#
# Step 3: Rerun the best mode with the corresponding starting values & and to avoid label switching by order the class sizes using these values
#	probs.start.new <- poLCA.reorder(bestmodel_mm$probs.start,order(bestmodel_mm$P,decreasing=TRUE))
#	bestmodel_mm<- poLCA(f,train_prac_mm,nclass=5, maxiter=5000,  # for example, best model is class=5
#                     probs.start=probs.start.new,nrep=1,
#                     tol=1e-5,na.rm=FALSE,verbose=TRUE, calc.se=TRUE)
#
#	lcmodel_mm <- reshape2::melt(bestmodel_mm$probs, level=2)
#################################################################################

createModels("H:\\Mplus_MM\\unsup_patidstrata\\1844_template.txt")
createModels("H:\\Mplus_MM\\unsup_patidstrata\\4564_template.txt")
createModels("H:\\Mplus_MM\\unsup_patidstrata\\6584_template.txt")
createModels("H:\\Mplus_MM\\unsup_patidstrata\\85plus_template.txt")

# first check if any morbidity is NOT present in the data (e.g. no young ppl has dementia!)
colSums(train_pat_mm_1844[,colmorbidlist]==2) # 2=TRUE # 
colSums(train_pat_mm_4564[,colmorbidlist]==2) # 2=TRUE # 
colSums(train_pat_mm_6584[,colmorbidlist]==2) # 2=TRUE # 
colSums(train_pat_mm_85plus[,colmorbidlist]==2) # 2=TRUE # 

colSums(test_pat_mm_1844[,colmorbidlist]==2) # 2=TRUE # 
colSums(test_pat_mm_4564[,colmorbidlist]==2) # 2=TRUE # 
colSums(test_pat_mm_6584[,colmorbidlist]==2) # 2=TRUE # 
colSums(test_pat_mm_85plus[,colmorbidlist]==2) # 2=TRUE # note, all c_mig=0 in the test data

runModels("H:\\Mplus_MM\\unsup_patidstrata",recursive=TRUE) 

## NOW pull OUTPUT FILES from high performance computer to local PC ##

# To ensure convergence
# we use a larger number of random sets (8000) of starting values. 
# Restrict analysis to 1-10 classes as initial results show higher class number lead to <1% class size
# After training and test sessions are all completed, copy all files in H:\ to C:\\MM_proc_data\\clusterMM

# MODELS# 

#####################
# age strata 18-44 # 
#####################

#extract model-fit and class separation statistics#
a<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\1844",
              filefilter="lca-unsup1844.*",what="summaries")
a_unsup<- do.call("rbind.fill",sapply(a,"[", "summaries"))

# combining results to a dataframe
sum_unsup<-a_unsup[,c("Parameters","LL","BIC","aBIC", "Entropy")]
sum_unsup$Model<-as.integer(a_unsup$Parameters/38)

sum_unsup<-sum_unsup[order(sum_unsup$Model),]
lca_results2_unsup1844<-sum_unsup
ztable::ztable(lca_results2_unsup1844) # this was poLCA users' last line for model derivation. Within each strata, the following lines are common to both Mplus and poLCA users.


#convert to long format
lca_results2_unsup1844

lca_results2_unsup1844<-lca_results2_unsup1844 %>%
  mutate(likelihood_ratio = (-2)*(dplyr::lag(LL)-LL))

lca_results2_unsup1844<-lca_results2_unsup1844[,c("Model","Parameters","LL","BIC","aBIC","likelihood_ratio","Entropy")]

sum2_unsup<-tidyr::gather(lca_results2_unsup1844,fitstatsname,fitstatsvalue,4:7)
sum2_unsup[,c("Model","Parameters","LL","fitstatsvalue")] <- 
  lapply(sum2_unsup[,c("Model","Parameters","LL","fitstatsvalue")], as.numeric)

#plot

setwd("C:\\MM_proc_data\\Results")
pdf(file="C:\\MM_proc_data\\Results\\LCAselectionstats_unsup1844_patid_1-10.pdf")
lcafit<-ggplot(sum2_unsup) + 
  geom_point(aes(x=Model,y=fitstatsvalue),size=3) +
  geom_line(aes(Model, fitstatsvalue, group = 1)) +
  scale_x_discrete("# of classes", limits=as.character(sum2_unsup$Model))+
  theme_bw()+
  labs(y="Statistics", 
       title = "Model selection statistics",
       subtitle = "Fit statistics and classification quality") + 
  facet_grid(fitstatsname ~. ,scales = "free") +
  theme_bw(base_size = 10, base_family = "") + 
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        legend.title = element_text(size = 10, face = 'bold'),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text=  element_text(size=10),
        axis.line = element_line(colour = "black")) # thicker line
lcafit
dev.off()

# 2. Plot the best model's profile (conditional probabilities) #
# Need to reorder class to fix the label #
# Best model =5-class
modelParam<- readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\1844\\lca-unsup1844_5.out",recursive=T)
a<-modelParam$parameters$probability.scale
a$param<-tolower(a$param)

# extract estimaed prob. from Mplus#
lcmodel_mm_train_unsup1844<-cbind(colmorbidlist,as.data.frame(mplus.get.estimated_probabilities("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\1844\\lca-unsup1844_5.gh5","process1",2,2)))

# merge full morbidity names into lca results
colnames(lcmodel_mm_train_unsup1844) <- c("Morbidity","class 1:","class 2:","class 3:","class 4:","class 5:")
lcmodel_mm_train_unsup1844 <- merge(lcmodel_mm_train_unsup1844,prev[,c('Morbidity','Morbidityfull')],by='Morbidity',all.x=TRUE)
lcmodel_mm_train_unsup1844<-melt(lcmodel_mm_train_unsup1844)
lcmodel_mm_train_unsup1844$clslabel<-as.numeric(lcmodel_mm_train_unsup1844$variable)

classsize_unsup1844<-modelParam$class_counts$mostLikely[,c("class","proportion")]
colnames(classsize_unsup1844)<-c("clslabel","clsprop")

lcmodel_mm_train_unsup1844<- merge(lcmodel_mm_train_unsup1844,classsize_unsup1844[,c("clsprop","clslabel")],by='clslabel')

lcmodel_mm_train_unsup1844$Var1<-paste("class:",lcmodel_mm_train_unsup1844$clslabel,
                                       round(100*lcmodel_mm_train_unsup1844$clsprop,2), "%")
classsize_unsup1844$clsprop<-round(100*classsize_unsup1844$clsprop,2)

pdf(file="C:\\MM_proc_data\\Results\\condplot1844_mm_morbid_unsup_patid.pdf")
condplot <- ggplot(lcmodel_mm_train_unsup1844,aes(x = Morbidityfull, y = value)) +
  geom_bar(stat = "identity")+
  facet_grid(clslabel ~ .)+
  labs(x = "Morbidities",y="Probabilities given class membership",
       title="Plots of condtional probabilities of morbidity diagnosis")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0, hjust=1,size=6),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),                    
        panel.grid.major.y=element_blank())
guides(fill = guide_legend(reverse=TRUE))
condplot
dev.off()

# Interpret the classes: list top 3 most prevalent morbidities for each class (top3morb) #
# Note, later codes added another notation by using the top3 most DISTINCTIVE morbidites (top3morb_3) # 

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_train_unsup1844 #  
temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value)) %>% # sort pr(true) for each disease by descending order
  top_n(3, value) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>% 
  group_by(clslabel) %>%
  mutate(top3morb = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value)
classsize_unsup1844<- merge(classsize_unsup1844,temp1[,c('clslabel','value','top3morb')], by='clslabel')
colnames(classsize_unsup1844)[3]<-"pr(True)"
classsize_unsup1844<-classsize_unsup1844[order(-classsize_unsup1844$clsprop),]
classsize_unsup1844

## try new cluster labels ##

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_train_unsup1844

# merge with prev data for overall prevalence of each disease #
temp1 <- merge(temp1, prev_1844[,c("Morbidity","Prop.")],by="Morbidity")
colnames(temp1)[5]<-"rawvalue"
temp1$value2<-(temp1$rawvalue*100>temp1$Prop.)*temp1$rawvalue #only posterior prob. bigger than overall prevalence is counted
# alternative 
temp1$value3<-(temp1$rawvalue*100>temp1$Prop.)*(temp1$rawvalue*100-temp1$Prop.)/100 #difference in relative prevalence

temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value3)) %>% # sort pr(true) for each disease by descending order, use differece in prevalence!
  top_n(3, value3) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>%
  group_by(clslabel) %>%
  mutate(top3morb_3 = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value3)

a<- merge(classsize_unsup1844,temp1[,c('clslabel','value2',"value3",'top3morb_3')], by='clslabel')
colnames(a)[5]<-"pr(True)_2"
colnames(a)[6]<-"pr(True)_3"
a<-a[order(-a$clsprop),]
classsize_unsup1844<-a
##########

# Mplus results feed into R (poLCA users do not need this line. Results have been inclued in the poLCA object) #
results_train_unsup1844<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\1844\\lca-unsup1844_5.out", 
                                 what="savedata")$savedata

#####################
# age strata 45-64 #
#####################

#extract model-fit and class separation statistics#
a<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\4564",
              filefilter="lca-unsup4564.*",what="summaries")
a_unsup<- do.call("rbind.fill",sapply(a,"[", "summaries"))

# combining results to a dataframe
sum_unsup<-a_unsup[,c("Parameters","LL","BIC","aBIC", "Entropy")]
sum_unsup$Model<-as.integer(a_unsup$Parameters/38)

sum_unsup<-sum_unsup[order(sum_unsup$Model),]
lca_results2_unsup4564<-sum_unsup
ztable::ztable(lca_results2_unsup4564)


#convert to long format
lca_results2_unsup4564

lca_results2_unsup4564<-lca_results2_unsup4564 %>%
  mutate(likelihood_ratio = (-2)*(dplyr::lag(LL)-LL))

lca_results2_unsup4564<-lca_results2_unsup4564[,c("Model","Parameters","LL","BIC","aBIC","likelihood_ratio","Entropy")]

sum2_unsup<-tidyr::gather(lca_results2_unsup4564,fitstatsname,fitstatsvalue,4:7)
sum2_unsup[,c("Model","Parameters","LL","fitstatsvalue")] <- 
  lapply(sum2_unsup[,c("Model","Parameters","LL","fitstatsvalue")], as.numeric)

#plot

setwd("C:\\MM_proc_data\\Results")
pdf(file="C:\\MM_proc_data\\Results\\LCAselectionstats_unsup4564_patid_1-10.pdf")
lcafit<-ggplot(sum2_unsup) + 
  geom_point(aes(x=Model,y=fitstatsvalue),size=3) +
  geom_line(aes(Model, fitstatsvalue, group = 1)) +
  scale_x_discrete("# of classes", limits=as.character(sum2_unsup$Model))+
  theme_bw()+
  labs(y="Statistics", 
       title = "Model selection statistics",
       subtitle = "Fit statistics and classification quality") + 
  facet_grid(fitstatsname ~. ,scales = "free") +
  theme_bw(base_size = 10, base_family = "") + 
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        legend.title = element_text(size = 10, face = 'bold'),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text=  element_text(size=10),
        axis.line = element_line(colour = "black")) # thicker line
lcafit
dev.off()

# 2. Plot the best model's profile (conditional probabilities) #
# Need to reorder class to fix the label #
# Best model =5-class
modelParam<- readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\4564\\lca-unsup4564_5.out",recursive=T)
a<-modelParam$parameters$probability.scale
a$param<-tolower(a$param)

# extract estimaed prob. from Mplus#
lcmodel_mm_train_unsup4564<-cbind(colmorbidlist,as.data.frame(mplus.get.estimated_probabilities("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\4564\\lca-unsup4564_5.gh5","process1",2,2)))

# merge full morbidity names into lca results
colnames(lcmodel_mm_train_unsup4564) <- c("Morbidity","class 1:","class 2:","class 3:","class 4:","class 5:")
lcmodel_mm_train_unsup4564 <- merge(lcmodel_mm_train_unsup4564,prev[,c('Morbidity','Morbidityfull')],by='Morbidity',all.x=TRUE)
lcmodel_mm_train_unsup4564<-melt(lcmodel_mm_train_unsup4564)
lcmodel_mm_train_unsup4564$clslabel<-as.numeric(lcmodel_mm_train_unsup4564$variable)

classsize_unsup4564<-modelParam$class_counts$mostLikely[,c("class","proportion")]
colnames(classsize_unsup4564)<-c("clslabel","clsprop")

lcmodel_mm_train_unsup4564<- merge(lcmodel_mm_train_unsup4564,classsize_unsup4564[,c("clsprop","clslabel")],by='clslabel')

lcmodel_mm_train_unsup4564$Var1<-paste("class:",lcmodel_mm_train_unsup4564$clslabel,
                                       round(100*lcmodel_mm_train_unsup4564$clsprop,2), "%")
classsize_unsup4564$clsprop<-round(100*classsize_unsup4564$clsprop,2)

pdf(file="C:\\MM_proc_data\\Results\\condplot4564_mm_morbid_unsup_patid.pdf")
condplot <- ggplot(lcmodel_mm_train_unsup4564,aes(x = Morbidityfull, y = value)) +
  geom_bar(stat = "identity")+
  facet_grid(clslabel ~ .)+
  labs(x = "Morbidities",y="Probabilities given class membership",
       title="Plots of condtional probabilities of morbidity diagnosis")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0, hjust=1,size=6),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),                    
        panel.grid.major.y=element_blank())
guides(fill = guide_legend(reverse=TRUE))
condplot
dev.off()

# Interpret the classes: list top 3 morbidities for each class #

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_train_unsup4564 #  
temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value)) %>% # sort pr(true) for each disease by descending order
  top_n(3, value) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>% 
  group_by(clslabel) %>%
  mutate(top3morb = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value)
classsize_unsup4564<- merge(classsize_unsup4564,temp1[,c('clslabel','value','top3morb')], by='clslabel')
colnames(classsize_unsup4564)[3]<-"pr(True)"
classsize_unsup4564<-classsize_unsup4564[order(-classsize_unsup4564$clsprop),]
classsize_unsup4564

## try new cluster labels ##

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_train_unsup4564

# merge with prev data for overall prevalence of each disease #
temp1 <- merge(temp1, prev_4564[,c("Morbidity","Prop.")],by="Morbidity")
colnames(temp1)[5]<-"rawvalue"
temp1$value2<-(temp1$rawvalue*100>temp1$Prop.)*temp1$rawvalue #only posterior prob. bigger than overall prevalence is counted
# alternative 
temp1$value3<-(temp1$rawvalue*100>temp1$Prop.)*(temp1$rawvalue*100-temp1$Prop.)/100 #difference in relative prevalence

temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value3)) %>% # sort pr(true) for each disease by descending order, use differece in prevalence!
  top_n(3, value3) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>%
  group_by(clslabel) %>%
  mutate(top3morb_3 = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value3)

a<- merge(classsize_unsup4564,temp1[,c('clslabel','value2',"value3",'top3morb_3')], by='clslabel')
colnames(a)[5]<-"pr(True)_2"
colnames(a)[6]<-"pr(True)_3"
a<-a[order(-a$clsprop),]
classsize_unsup4564<-a
##########


# Mplus results feed into R #
results_train_unsup4564<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\4564\\lca-unsup4564_5.out", 
                                    what="savedata")$savedata

#####################
# age strata 65-84 #
#####################

#extract model-fit and class separation statistics#
a<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\6584",
              filefilter="lca-unsup6584.*",what="summaries")
a_unsup<- do.call("rbind.fill",sapply(a,"[", "summaries"))

# combining results to a dataframe
sum_unsup<-a_unsup[,c("Parameters","LL","BIC","aBIC", "Entropy")]
sum_unsup$Model<-as.integer(a_unsup$Parameters/38)

sum_unsup<-sum_unsup[order(sum_unsup$Model),]
lca_results2_unsup6584<-sum_unsup
ztable::ztable(lca_results2_unsup6584)


#convert to long format
lca_results2_unsup6584

lca_results2_unsup6584<-lca_results2_unsup6584 %>%
  mutate(likelihood_ratio = (-2)*(dplyr::lag(LL)-LL))

lca_results2_unsup6584<-lca_results2_unsup6584[,c("Model","Parameters","LL","BIC","aBIC","likelihood_ratio","Entropy")]

sum2_unsup<-tidyr::gather(lca_results2_unsup6584,fitstatsname,fitstatsvalue,4:7)
sum2_unsup[,c("Model","Parameters","LL","fitstatsvalue")] <- 
  lapply(sum2_unsup[,c("Model","Parameters","LL","fitstatsvalue")], as.numeric)

#plot

setwd("C:\\MM_proc_data\\Results")
pdf(file="C:\\MM_proc_data\\Results\\LCAselectionstats_unsup6584_patid_1-10.pdf")
lcafit<-ggplot(sum2_unsup) + 
  geom_point(aes(x=Model,y=fitstatsvalue),size=3) +
  geom_line(aes(Model, fitstatsvalue, group = 1)) +
  scale_x_discrete("# of classes", limits=as.character(sum2_unsup$Model))+
  theme_bw()+
  labs(y="Statistics", 
       title = "Model selection statistics",
       subtitle = "Fit statistics and classification quality") + 
  facet_grid(fitstatsname ~. ,scales = "free") +
  theme_bw(base_size = 10, base_family = "") + 
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        legend.title = element_text(size = 10, face = 'bold'),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text=  element_text(size=10),
        axis.line = element_line(colour = "black")) # thicker line
lcafit
dev.off()

# 2. Plot the best model's profile (conditional probabilities) #
# Need to reorder class to fix the label #
# Best model =6-class
modelParam<- readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\6584\\lca-unsup6584_6.out",recursive=T)
a<-modelParam$parameters$probability.scale
a$param<-tolower(a$param)

# extract estimaed prob. from Mplus#
lcmodel_mm_train_unsup6584<-cbind(colmorbidlist,as.data.frame(mplus.get.estimated_probabilities("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\6584\\lca-unsup6584_6.gh5","process1",2,2)))

# merge full morbidity names into lca results
colnames(lcmodel_mm_train_unsup6584) <- c("Morbidity","class 1:","class 2:","class 3:","class 4:","class 5:","class 6:")
lcmodel_mm_train_unsup6584 <- merge(lcmodel_mm_train_unsup6584,prev[,c('Morbidity','Morbidityfull')],by='Morbidity',all.x=TRUE)
lcmodel_mm_train_unsup6584<-melt(lcmodel_mm_train_unsup6584)
lcmodel_mm_train_unsup6584$clslabel<-as.numeric(lcmodel_mm_train_unsup6584$variable)

classsize_unsup6584<-modelParam$class_counts$mostLikely[,c("class","proportion")]
colnames(classsize_unsup6584)<-c("clslabel","clsprop")

lcmodel_mm_train_unsup6584<- merge(lcmodel_mm_train_unsup6584,classsize_unsup6584[,c("clsprop","clslabel")],by='clslabel')

lcmodel_mm_train_unsup6584$Var1<-paste("class:",lcmodel_mm_train_unsup6584$clslabel,
                                       round(100*lcmodel_mm_train_unsup6584$clsprop,2), "%")
classsize_unsup6584$clsprop<-round(100*classsize_unsup6584$clsprop,2)

pdf(file="C:\\MM_proc_data\\Results\\condplot6584_mm_morbid_unsup_patid.pdf")
condplot <- ggplot(lcmodel_mm_train_unsup6584,aes(x = Morbidityfull, y = value)) +
  geom_bar(stat = "identity")+
  facet_grid(clslabel ~ .)+
  labs(x = "Morbidities",y="Probabilities given class membership",
       title="Plots of condtional probabilities of morbidity diagnosis")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0, hjust=1,size=6),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),                    
        panel.grid.major.y=element_blank())
guides(fill = guide_legend(reverse=TRUE))
condplot
dev.off()

# Interpret the classes: list top 3 morbidities for each class #

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_train_unsup6584 #  
temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value)) %>% # sort pr(true) for each disease by descending order
  top_n(3, value) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>% 
  group_by(clslabel) %>%
  mutate(top3morb = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value)
classsize_unsup6584<- merge(classsize_unsup6584,temp1[,c('clslabel','value','top3morb')], by='clslabel')
colnames(classsize_unsup6584)[3]<-"pr(True)"
classsize_unsup6584<-classsize_unsup6584[order(-classsize_unsup6584$clsprop),]
classsize_unsup6584

## try new cluster labels ##

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_train_unsup6584

# merge with prev data for overall prevalence of each disease #
temp1 <- merge(temp1, prev_6584[,c("Morbidity","Prop.")],by="Morbidity")
colnames(temp1)[5]<-"rawvalue"
temp1$value2<-(temp1$rawvalue*100>temp1$Prop.)*temp1$rawvalue #only posterior prob. bigger than overall prevalence is counted
# alternative 
temp1$value3<-(temp1$rawvalue*100>temp1$Prop.)*(temp1$rawvalue*100-temp1$Prop.)/100 #difference in relative prevalence

temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value3)) %>% # sort pr(true) for each disease by descending order, use differece in prevalence!
  top_n(3, value3) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>%
  group_by(clslabel) %>%
  mutate(top3morb_3 = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value3)

a<- merge(classsize_unsup6584,temp1[,c('clslabel','value2',"value3",'top3morb_3')], by='clslabel')
colnames(a)[5]<-"pr(True)_2"
colnames(a)[6]<-"pr(True)_3"
a<-a[order(-a$clsprop),]
classsize_unsup6584<-a
##########


# Mplus results feed into R #
results_train_unsup6584<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\6584\\lca-unsup6584_6.out", 
                                    what="savedata")$savedata

#####################
# age strata 85+ #
#####################

#extract model-fit and class separation statistics#
a<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\85plus",
              filefilter="lca-unsup85plus.*",what="summaries")
a_unsup<- do.call("rbind.fill",sapply(a,"[", "summaries"))

# combining results to a dataframe
sum_unsup<-a_unsup[,c("Parameters","LL","BIC","aBIC", "Entropy")]
sum_unsup$Model<-as.integer(a_unsup$Parameters/38)

sum_unsup<-sum_unsup[order(sum_unsup$Model),]
lca_results2_unsup85plus<-sum_unsup
ztable::ztable(lca_results2_unsup85plus)


#convert to long format
lca_results2_unsup85plus

lca_results2_unsup85plus<-lca_results2_unsup85plus %>%
  mutate(likelihood_ratio = (-2)*(dplyr::lag(LL)-LL))

lca_results2_unsup85plus<-lca_results2_unsup85plus[,c("Model","Parameters","LL","BIC","aBIC","likelihood_ratio","Entropy")]

sum2_unsup<-tidyr::gather(lca_results2_unsup85plus,fitstatsname,fitstatsvalue,4:7)
sum2_unsup[,c("Model","Parameters","LL","fitstatsvalue")] <- 
  lapply(sum2_unsup[,c("Model","Parameters","LL","fitstatsvalue")], as.numeric)

#plot

setwd("C:\\MM_proc_data\\Results")
pdf(file="C:\\MM_proc_data\\Results\\LCAselectionstats_unsup85plus_patid_1-10.pdf")
lcafit<-ggplot(sum2_unsup) + 
  geom_point(aes(x=Model,y=fitstatsvalue),size=3) +
  geom_line(aes(Model, fitstatsvalue, group = 1)) +
  scale_x_discrete("# of classes", limits=as.character(sum2_unsup$Model))+
  theme_bw()+
  labs(y="Statistics", 
       title = "Model selection statistics",
       subtitle = "Fit statistics and classification quality") + 
  facet_grid(fitstatsname ~. ,scales = "free") +
  theme_bw(base_size = 10, base_family = "") + 
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        legend.title = element_text(size = 10, face = 'bold'),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text=  element_text(size=10),
        axis.line = element_line(colour = "black")) # thicker line
lcafit
dev.off()

# 2. Plot the best model's profile (conditional probabilities) #
# Need to reorder class to fix the label #
# Best model =4-class
modelParam<- readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\85plus\\lca-unsup85plus_4.out",recursive=T)
a<-modelParam$parameters$probability.scale
a$param<-tolower(a$param)

# extract estimaed prob. from Mplus#
lcmodel_mm_train_unsup85plus<-cbind(colmorbidlist,as.data.frame(mplus.get.estimated_probabilities("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\85plus\\lca-unsup85plus_4.gh5","process1",2,2)))

# merge full morbidity names into lca results
colnames(lcmodel_mm_train_unsup85plus) <- c("Morbidity","class 1:","class 2:","class 3:","class 4:")
lcmodel_mm_train_unsup85plus <- merge(lcmodel_mm_train_unsup85plus,prev[,c('Morbidity','Morbidityfull')],by='Morbidity',all.x=TRUE)
lcmodel_mm_train_unsup85plus<-melt(lcmodel_mm_train_unsup85plus)
lcmodel_mm_train_unsup85plus$clslabel<-as.numeric(lcmodel_mm_train_unsup85plus$variable)

classsize_unsup85plus<-modelParam$class_counts$mostLikely[,c("class","proportion")]
colnames(classsize_unsup85plus)<-c("clslabel","clsprop")

lcmodel_mm_train_unsup85plus<- merge(lcmodel_mm_train_unsup85plus,classsize_unsup85plus[,c("clsprop","clslabel")],by='clslabel')

lcmodel_mm_train_unsup85plus$Var1<-paste("class:",lcmodel_mm_train_unsup85plus$clslabel,
                                       round(100*lcmodel_mm_train_unsup85plus$clsprop,2), "%")
classsize_unsup85plus$clsprop<-round(100*classsize_unsup85plus$clsprop,2)

pdf(file="C:\\MM_proc_data\\Results\\condplot85plus_mm_morbid_unsup_patid.pdf")
condplot <- ggplot(lcmodel_mm_train_unsup85plus,aes(x = Morbidityfull, y = value)) +
  geom_bar(stat = "identity")+
  facet_grid(clslabel ~ .)+
  labs(x = "Morbidities",y="Probabilities given class membership",
       title="Plots of condtional probabilities of morbidity diagnosis")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0, hjust=1,size=6),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),                    
        panel.grid.major.y=element_blank())
guides(fill = guide_legend(reverse=TRUE))
condplot
dev.off()

# Interpret the classes: list top 3 morbidities for each class #

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_train_unsup85plus #  
temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value)) %>% # sort pr(true) for each disease by descending order
  top_n(3, value) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>% 
  group_by(clslabel) %>%
  mutate(top3morb = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value)
classsize_unsup85plus<- merge(classsize_unsup85plus,temp1[,c('clslabel','value','top3morb')], by='clslabel')
colnames(classsize_unsup85plus)[3]<-"pr(True)"
classsize_unsup85plus<-classsize_unsup85plus[order(-classsize_unsup85plus$clsprop),]
classsize_unsup85plus

## try new cluster labels ##

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_train_unsup85plus

# merge with prev data for overall prevalence of each disease #
temp1 <- merge(temp1, prev_85plus[,c("Morbidity","Prop.")],by="Morbidity")
colnames(temp1)[5]<-"rawvalue"
temp1$value2<-(temp1$rawvalue*100>temp1$Prop.)*temp1$rawvalue #only posterior prob. bigger than overall prevalence is counted
# alternative 
temp1$value3<-(temp1$rawvalue*100>temp1$Prop.)*(temp1$rawvalue*100-temp1$Prop.)/100 #difference in relative prevalence

temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value3)) %>% # sort pr(true) for each disease by descending order, use differece in prevalence!
  top_n(3, value3) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>%
  group_by(clslabel) %>%
  mutate(top3morb_3 = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value3)

a<- merge(classsize_unsup85plus,temp1[,c('clslabel','value2',"value3",'top3morb_3')], by='clslabel')
colnames(a)[5]<-"pr(True)_2"
colnames(a)[6]<-"pr(True)_3"
a<-a[order(-a$clsprop),]
classsize_unsup85plus<-a
##########


# Mplus results feed into R #
results_train_unsup85plus<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\85plus\\lca-unsup85plus_4.out", 
                                    what="savedata")$savedata

##########################################################################################################################

# Tests clustering quality on test data #
# first check if any morbidity is NOT present in the data (e.g. no young ppl has dementia?)
colSums(test_pat_mm_1844[,colmorbidlist]==2) # 2=TRUE # 
colSums(test_pat_mm_4564[,colmorbidlist]==2) # 2=TRUE # 
colSums(test_pat_mm_6584[,colmorbidlist]==2) # 2=TRUE # 
colSums(test_pat_mm_85plus[,colmorbidlist]==2) # 2=TRUE # c_mig all 0


# runModels in high performane computing:  outcome_test folder in H:\ #

#####################
# age strata 18-44, 5-class model#
#####################

# Check class profiles in the test data using model with best-class#
modelParam<- readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup1844_5.out",recursive=T)
a<-modelParam$parameters$probability.scale
a$param<-tolower(a$param)

# Open Mplus to obtain data #
lcmodel_mm_test_unsup1844<-cbind(colmorbidlist,as.data.frame(mplus.get.estimated_probabilities("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup1844_5.gh5","process1",2,2)))

# merge full morbidity names into lca results
colnames(lcmodel_mm_test_unsup1844) <- c("Morbidity","class 1:","class 2:","class 3:","class 4:","class 5:")
lcmodel_mm_test_unsup1844 <- merge(lcmodel_mm_test_unsup1844,prev[,c('Morbidity','Morbidityfull')],by='Morbidity',all.x=TRUE)
lcmodel_mm_test_unsup1844<-melt(lcmodel_mm_test_unsup1844)
lcmodel_mm_test_unsup1844$clslabel<-as.numeric(lcmodel_mm_test_unsup1844$variable)

classsize_test_unsup1844<-modelParam$class_counts$mostLikely[,c("class","proportion")]
colnames(classsize_test_unsup1844)<-c("clslabel","clsprop")

lcmodel_mm_test_unsup1844<- merge(lcmodel_mm_test_unsup1844,classsize_test_unsup1844[,c("clsprop","clslabel")],by='clslabel')

lcmodel_mm_test_unsup1844$Var1<-paste("class:",lcmodel_mm_test_unsup1844$clslabel,round(100*lcmodel_mm_test_unsup1844$clsprop,2), "%")
classsize_test_unsup1844$clsprop<-round(100*classsize_test_unsup1844$clsprop,2)

pdf(file="C:\\MM_proc_data\\Results\\condplot1844_mm_morbid_unsuptest_patid.pdf")
condplot <- ggplot(lcmodel_mm_test_unsup1844,aes(x = Morbidityfull, y = value)) +
  geom_bar(stat = "identity")+
  facet_grid(clslabel ~ .)+
  labs(x = "Morbidities",y="Probabilities given class membership",
       title="Plots of condtional probabilities of morbidity diagnosis")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0, hjust=1,size=6),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),                    
        panel.grid.major.y=element_blank())
guides(fill = guide_legend(reverse=TRUE))
condplot
dev.off()

# Interpret the classes: list top 3 morbidities for each class #

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_test_unsup1844 
temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value)) %>% # sort pr(true) for each disease by descending order
  top_n(3, value) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>% 
  group_by(clslabel) %>%
  mutate(top3morb = paste(Morbidityfull, collapse =",")) %>%
  top_n(1,value)
classsize_test_unsup1844<- merge(classsize_test_unsup1844,temp1[,c('clslabel','value','top3morb')], by='clslabel')
colnames(classsize_test_unsup1844)[3]<-"pr(True)"
classsize_test_unsup1844<-classsize_test_unsup1844[order(-classsize_test_unsup1844$clsprop),]
classsize_test_unsup1844

## try new cluster labels ##

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_test_unsup1844

# merge with prev data for overall prevalence of each disease #
temp1 <- merge(temp1, prev_1844[,c("Morbidity","Prop.")],by="Morbidity")
colnames(temp1)[5]<-"rawvalue"
temp1$value2<-(temp1$rawvalue*100>temp1$Prop.)*temp1$rawvalue #only posterior prob. bigger than overall prevalence is counted
# alternative 
temp1$value3<-(temp1$rawvalue*100>temp1$Prop.)*(temp1$rawvalue*100-temp1$Prop.)/100 #difference in relative prevalence

temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value3)) %>% # sort pr(true) for each disease by descending order, use differece in prevalence!
  top_n(3, value3) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>%
  group_by(clslabel) %>%
  mutate(top3morb_3 = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value3)

a<- merge(classsize_test_unsup1844,temp1[,c('clslabel','value2',"value3",'top3morb_3')], by='clslabel')
colnames(a)[5]<-"pr(True)_2"
colnames(a)[6]<-"pr(True)_3"
a<-a[order(-a$clsprop),]
classsize_test_unsup1844<-a
##########


# Mplus results feed into R #
results_test_unsup1844<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup1844_5.out", 
                                  what="savedata")$savedata

#####################
# age strata 45-64#
#####################

# Check class profiles in the test data using model with best-class#
modelParam<- readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup4564_5.out",recursive=T)
a<-modelParam$parameters$probability.scale
a$param<-tolower(a$param)

# Open Mplus to obtain data #
lcmodel_mm_test_unsup4564<-cbind(colmorbidlist,as.data.frame(mplus.get.estimated_probabilities("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup4564_5.gh5","process1",2,2)))

# merge full morbidity names into lca results
colnames(lcmodel_mm_test_unsup4564) <- c("Morbidity","class 1:","class 2:","class 3:","class 4:","class 5:")
lcmodel_mm_test_unsup4564 <- merge(lcmodel_mm_test_unsup4564,prev[,c('Morbidity','Morbidityfull')],by='Morbidity',all.x=TRUE)
lcmodel_mm_test_unsup4564<-melt(lcmodel_mm_test_unsup4564)
lcmodel_mm_test_unsup4564$clslabel<-as.numeric(lcmodel_mm_test_unsup4564$variable)

classsize_test_unsup4564<-modelParam$class_counts$mostLikely[,c("class","proportion")]
colnames(classsize_test_unsup4564)<-c("clslabel","clsprop")

lcmodel_mm_test_unsup4564<- merge(lcmodel_mm_test_unsup4564,classsize_test_unsup4564[,c("clsprop","clslabel")],by='clslabel')

lcmodel_mm_test_unsup4564$Var1<-paste("class:",lcmodel_mm_test_unsup4564$clslabel,round(100*lcmodel_mm_test_unsup4564$clsprop,2), "%")
classsize_test_unsup4564$clsprop<-round(100*classsize_test_unsup4564$clsprop,2)

pdf(file="C:\\MM_proc_data\\Results\\condplot4564_mm_morbid_unsuptest_patid.pdf")
condplot <- ggplot(lcmodel_mm_test_unsup4564,aes(x = Morbidityfull, y = value)) +
  geom_bar(stat = "identity")+
  facet_grid(clslabel ~ .)+
  labs(x = "Morbidities",y="Probabilities given class membership",
       title="Plots of condtional probabilities of morbidity diagnosis")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0, hjust=1,size=6),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),                    
        panel.grid.major.y=element_blank())
guides(fill = guide_legend(reverse=TRUE))
condplot
dev.off()

# Interpret the classes: list top 3 morbidities for each class #

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_test_unsup4564 
temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value)) %>% # sort pr(true) for each disease by descending order
  top_n(3, value) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>% 
  group_by(clslabel) %>%
  mutate(top3morb = paste(Morbidityfull, collapse =",")) %>%
  top_n(1,value)
classsize_test_unsup4564<- merge(classsize_test_unsup4564,temp1[,c('clslabel','value','top3morb')], by='clslabel')
colnames(classsize_test_unsup4564)[3]<-"pr(True)"
classsize_test_unsup4564<-classsize_test_unsup4564[order(-classsize_test_unsup4564$clsprop),]
classsize_test_unsup4564

## try new cluster labels ##

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_test_unsup4564

# merge with prev data for overall prevalence of each disease #
temp1 <- merge(temp1, prev_4564[,c("Morbidity","Prop.")],by="Morbidity")
colnames(temp1)[5]<-"rawvalue"
temp1$value2<-(temp1$rawvalue*100>temp1$Prop.)*temp1$rawvalue #only posterior prob. bigger than overall prevalence is counted
# alternative 
temp1$value3<-(temp1$rawvalue*100>temp1$Prop.)*(temp1$rawvalue*100-temp1$Prop.)/100 #difference in relative prevalence

temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value3)) %>% # sort pr(true) for each disease by descending order, use differece in prevalence!
  top_n(3, value3) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>%
  group_by(clslabel) %>%
  mutate(top3morb_3 = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value3)

a<- merge(classsize_test_unsup4564,temp1[,c('clslabel','value2',"value3",'top3morb_3')], by='clslabel')
colnames(a)[5]<-"pr(True)_2"
colnames(a)[6]<-"pr(True)_3"
a<-a[order(-a$clsprop),]
classsize_test_unsup4564<-a
##########


# Mplus results feed into R #
results_test_unsup4564<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup4564_5.out", 
                                   what="savedata")$savedata

#####################
# age strata 65-84#
#####################

# Check class profiles in the test data using model with best-class#
modelParam<- readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup6584_6.out",recursive=T)
a<-modelParam$parameters$probability.scale
a$param<-tolower(a$param)

# Open Mplus to obtain data #
lcmodel_mm_test_unsup6584<-cbind(colmorbidlist,as.data.frame(mplus.get.estimated_probabilities("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup6584_6.gh5","process1",2,2)))

# merge full morbidity names into lca results
colnames(lcmodel_mm_test_unsup6584) <- c("Morbidity","class 1:","class 2:","class 3:","class 4:","class 5:","class 6:")
lcmodel_mm_test_unsup6584 <- merge(lcmodel_mm_test_unsup6584,prev[,c('Morbidity','Morbidityfull')],by='Morbidity',all.x=TRUE)
lcmodel_mm_test_unsup6584<-melt(lcmodel_mm_test_unsup6584)
lcmodel_mm_test_unsup6584$clslabel<-as.numeric(lcmodel_mm_test_unsup6584$variable)

classsize_test_unsup6584<-modelParam$class_counts$mostLikely[,c("class","proportion")]
colnames(classsize_test_unsup6584)<-c("clslabel","clsprop")

lcmodel_mm_test_unsup6584<- merge(lcmodel_mm_test_unsup6584,classsize_test_unsup6584[,c("clsprop","clslabel")],by='clslabel')

lcmodel_mm_test_unsup6584$Var1<-paste("class:",lcmodel_mm_test_unsup6584$clslabel,round(100*lcmodel_mm_test_unsup6584$clsprop,2), "%")
classsize_test_unsup6584$clsprop<-round(100*classsize_test_unsup6584$clsprop,2)

pdf(file="C:\\MM_proc_data\\Results\\condplot6584_mm_morbid_unsuptest_patid.pdf")
condplot <- ggplot(lcmodel_mm_test_unsup6584,aes(x = Morbidityfull, y = value)) +
  geom_bar(stat = "identity")+
  facet_grid(clslabel ~ .)+
  labs(x = "Morbidities",y="Probabilities given class membership",
       title="Plots of condtional probabilities of morbidity diagnosis")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0, hjust=1,size=6),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),                    
        panel.grid.major.y=element_blank())
guides(fill = guide_legend(reverse=TRUE))
condplot
dev.off()

# Interpret the classes: list top 3 morbidities for each class #

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_test_unsup6584 
temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value)) %>% # sort pr(true) for each disease by descending order
  top_n(3, value) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>% 
  group_by(clslabel) %>%
  mutate(top3morb = paste(Morbidityfull, collapse =",")) %>%
  top_n(1,value)
classsize_test_unsup6584<- merge(classsize_test_unsup6584,temp1[,c('clslabel','value','top3morb')], by='clslabel')
colnames(classsize_test_unsup6584)[3]<-"pr(True)"
classsize_test_unsup6584<-classsize_test_unsup6584[order(-classsize_test_unsup6584$clsprop),]
classsize_test_unsup6584

## try new cluster labels ##

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_test_unsup6584

# merge with prev data for overall prevalence of each disease #
temp1 <- merge(temp1, prev_6584[,c("Morbidity","Prop.")],by="Morbidity")
colnames(temp1)[5]<-"rawvalue"
temp1$value2<-(temp1$rawvalue*100>temp1$Prop.)*temp1$rawvalue #only posterior prob. bigger than overall prevalence is counted
# alternative 
temp1$value3<-(temp1$rawvalue*100>temp1$Prop.)*(temp1$rawvalue*100-temp1$Prop.)/100 #difference in relative prevalence

temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value3)) %>% # sort pr(true) for each disease by descending order, use differece in prevalence!
  top_n(3, value3) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>%
  group_by(clslabel) %>%
  mutate(top3morb_3 = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value3)

a<- merge(classsize_test_unsup6584,temp1[,c('clslabel','value2',"value3",'top3morb_3')], by='clslabel')
colnames(a)[5]<-"pr(True)_2"
colnames(a)[6]<-"pr(True)_3"
a<-a[order(-a$clsprop),]
classsize_test_unsup6584<-a
##########


# Mplus results feed into R #
results_test_unsup6584<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup6584_6.out", 
                                   what="savedata")$savedata

#####################
# age strata 85+#
#####################

# Check class profiles in the test data using model with best-class#
modelParam<- readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup85plus_4.out",recursive=T)
a<-modelParam$parameters$probability.scale
a$param<-tolower(a$param)

# Open Mplus to obtain data # -36 for c_mig =0
lcmodel_mm_test_unsup85plus<-cbind(colmorbidlist[-36],as.data.frame(mplus.get.estimated_probabilities("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup85plus_4.gh5","process1",2,2)))

# merge full morbidity names into lca results
colnames(lcmodel_mm_test_unsup85plus) <- c("Morbidity","class 1:","class 2:","class 3:","class 4:")
lcmodel_mm_test_unsup85plus <- merge(lcmodel_mm_test_unsup85plus,prev[,c('Morbidity','Morbidityfull')],by='Morbidity',all.x=TRUE)
lcmodel_mm_test_unsup85plus<-melt(lcmodel_mm_test_unsup85plus)
lcmodel_mm_test_unsup85plus$clslabel<-as.numeric(lcmodel_mm_test_unsup85plus$variable)

classsize_test_unsup85plus<-modelParam$class_counts$mostLikely[,c("class","proportion")]
colnames(classsize_test_unsup85plus)<-c("clslabel","clsprop")

lcmodel_mm_test_unsup85plus<- merge(lcmodel_mm_test_unsup85plus,classsize_test_unsup85plus[,c("clsprop","clslabel")],by='clslabel')

lcmodel_mm_test_unsup85plus$Var1<-paste("class:",lcmodel_mm_test_unsup85plus$clslabel,round(100*lcmodel_mm_test_unsup85plus$clsprop,2), "%")
classsize_test_unsup85plus$clsprop<-round(100*classsize_test_unsup85plus$clsprop,2)

pdf(file="C:\\MM_proc_data\\Results\\condplot85plus_mm_morbid_unsuptest_patid.pdf")
condplot <- ggplot(lcmodel_mm_test_unsup85plus,aes(x = Morbidityfull, y = value)) +
  geom_bar(stat = "identity")+
  facet_grid(clslabel ~ .)+
  labs(x = "Morbidities",y="Probabilities given class membership",
       title="Plots of condtional probabilities of morbidity diagnosis")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0, hjust=1,size=6),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),                    
        panel.grid.major.y=element_blank())
guides(fill = guide_legend(reverse=TRUE))
condplot
dev.off()

# Interpret the classes: list top 3 morbidities for each class #

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_test_unsup85plus 
temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value)) %>% # sort pr(true) for each disease by descending order
  top_n(3, value) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>% 
  group_by(clslabel) %>%
  mutate(top3morb = paste(Morbidityfull, collapse =",")) %>%
  top_n(1,value)
classsize_test_unsup85plus<- merge(classsize_test_unsup85plus,temp1[,c('clslabel','value','top3morb')], by='clslabel')
colnames(classsize_test_unsup85plus)[3]<-"pr(True)"
classsize_test_unsup85plus<-classsize_test_unsup85plus[order(-classsize_test_unsup85plus$clsprop),]
classsize_test_unsup85plus

## try new cluster labels ##

# order pr(true) for morbidities within each class
temp1<- lcmodel_mm_test_unsup85plus

# merge with prev data for overall prevalence of each disease #
temp1 <- merge(temp1, prev_85plus[,c("Morbidity","Prop.")],by="Morbidity")
colnames(temp1)[5]<-"rawvalue"
temp1$value2<-(temp1$rawvalue*100>temp1$Prop.)*temp1$rawvalue #only posterior prob. bigger than overall prevalence is counted
# alternative 
temp1$value3<-(temp1$rawvalue*100>temp1$Prop.)*(temp1$rawvalue*100-temp1$Prop.)/100 #difference in relative prevalence

temp1<-temp1 %>%
  group_by(clslabel) %>%
  arrange(clslabel,desc(value3)) %>% # sort pr(true) for each disease by descending order, use differece in prevalence!
  top_n(3, value3) # select top 3 rows

# Report top 3 morbidities for each class
temp1<-temp1 %>%
  group_by(clslabel) %>%
  mutate(top3morb_3 = paste0(Morbidityfull, collapse =",")) %>%
  top_n(1,value3)

a<- merge(classsize_test_unsup85plus,temp1[,c('clslabel','value2',"value3",'top3morb_3')], by='clslabel')
colnames(a)[5]<-"pr(True)_2"
colnames(a)[6]<-"pr(True)_3"
a<-a[order(-a$clsprop),]
classsize_test_unsup85plus<-a
##########


# Mplus results feed into R #
results_test_unsup85plus<-readModels("C:\\MM_proc_data\\clusterMM\\unsup_patidstrata\\test\\lca-unsup85plus_4.out", 
                                   what="savedata")$savedata

save.image(file='C:\\MM_proc_data\\final_patid_strata_final.RData') 

########################################################################################################
 ## Some added code to filter REAL GP consultations (blue codes + code 61) ##
load(file="C:\\MM_proc_data\\final_patid_strata_final.RData") 

load.ffdf('C:\\MM_proc_data\\ffdb\\consultation')
consultation_relevant <- as.data.table.ffdf(consultation[consultation$eventdate >= '2012-01-01'&consultation$eventdate <= '2012-12-31',])
filterlist<-c(1,2,3,4,6,7,8,9,10,11,18,21,22,27,28,30,31,32,33,34,35,36,37,38,50,55)
consultation_relevant<-consultation_relevant %>%
  group_by(patid) %>%          
  mutate(ctconsult12mons_filterdur=uniqueN(as.character(eventdate)[constype %in% filterlist]),
         ctconsult12mons_nofilterdur=uniqueN(as.character(eventdate)))
consultation_relevant<- consultation_relevant[!duplicated(consultation_relevant[,'patid']),]
########################################################################################################

## Test sets ##
 
# a_test<- merge(patient2012_copy,a[,c('patid','ctconsult12mons_filterdur')],by='patid',all.x=T)
# a_test$ctconsult12mons_filterdur[is.na(a_test$ctconsult12mons_filterdur)]<-0
# 
# tapply(a_test$ctconsult12mons_filterdur,a_test$binary_multimorb,summary)
# tapply(a_test$ctprod12mons4plus_BNF,a_test$binary_multimorb,summary)
# tapply(a_test$cthospspl_12mons,a_test$binary_multimorb,summary)


# merge this file with ALL current results and data files # 
 dflist<-list(train_pat_mm_1844,train_pat_mm_4564,train_pat_mm_6584,train_pat_mm_85plus,
              test_pat_mm_1844,test_pat_mm_4564,test_pat_mm_6584,test_pat_mm_85plus,
              results_train_unsup1844,results_train_unsup4564,results_train_unsup6584,
              results_train_unsup85plus,
              results_test_unsup1844,results_test_unsup4564,results_test_unsup6584,
              results_test_unsup85plus,patient2012,patient2012_copy,patient2012_copy_mm,
              lcapatid_train_unsup1844,lcapatid_train_unsup4564,lcapatid_train_unsup6584,
              lcapatid_train_unsup85plus,
              lcapatid_test_unsup1844,lcapatid_test_unsup4564,lcapatid_test_unsup6584,
              lcapatid_test_unsup85plus)
 name.dflist<-c("train_pat_mm_1844","train_pat_mm_4564","train_pat_mm_6584","train_pat_mm_85plus",
                "test_pat_mm_1844","test_pat_mm_4564","test_pat_mm_6584","test_pat_mm_85plus",
                "results_train_unsup1844","results_train_unsup4564","results_train_unsup6584",
                "results_train_unsup85plus",
                "results_test_unsup1844","results_test_unsup4564","results_test_unsup6584",
                "results_test_unsup85plus","patient2012","patient2012_copy","patient2012_copy_mm",
                "lcapatid_train_unsup1844","lcapatid_train_unsup4564","lcapatid_train_unsup6584",
                "lcapatid_train_unsup85plus",
                "lcapatid_test_unsup1844","lcapatid_test_unsup4564","lcapatid_test_unsup6584",
                "lcapatid_test_unsup85plus")

 for (i in 1:length(dflist)){
   print(i)
   dflist[[i]]<- merge(dflist[[i]],consultation_relevant[,c('patid','ctconsult12mons_filterdur')],by='patid',all.x=T)
   dflist[[i]]$ctconsult12mons_filterdur[is.na(dflist[[i]]$ctconsult12mons_filterdur)]<-0
   assign(name.dflist[i], dflist[[i]])
 }

############################################################################################
# Regression models to explore the relationship between demographics, 
# clusters and outcomes. 

# merge each results_train_unsup2 dataset with demographics 
demogrlist<-c("gender","imd2015_5","smokstatus","agegroup")

# merge results into the original dataset #

# 18-44 TRAIN #
colnames(results_train_unsup1844)[44:45]<-c("clslabel","patid")
b<-subset(results_train_unsup1844,!duplicated(results_train_unsup1844$patid))
a<-merge(train_pat_mm_1844,b[,c("CPROB1", "CPROB2", "CPROB3", "CPROB4",
                           "CPROB5","clslabel","patid")],by="patid")
a<-merge(a,classsize_unsup1844,by="clslabel")

colnames(a)[1]<-"modalcls"
lcapatid_train_unsup1844<-a  # Note: train, trainpat are no longer the correct indicators for this stratified split
rm(a,b)

# 45-64 TRAIN #
colnames(results_train_unsup4564)[44:45]<-c("clslabel","patid")
b<-subset(results_train_unsup4564,!duplicated(results_train_unsup4564$patid))
a<-merge(train_pat_mm_4564,b[,c("CPROB1", "CPROB2", "CPROB3", "CPROB4",
                                "CPROB5","clslabel","patid")],by="patid")
a<-merge(a,classsize_unsup4564,by="clslabel")

colnames(a)[1]<-"modalcls"
lcapatid_train_unsup4564<-a  # Note: train, trainpat are no longer the correct indicators for this stratified split
rm(a,b)

# 65-84 TRAIN #
colnames(results_train_unsup6584)[45:46]<-c("clslabel","patid")
b<-subset(results_train_unsup6584,!duplicated(results_train_unsup6584$patid))
a<-merge(train_pat_mm_6584,b[,c("CPROB1", "CPROB2", "CPROB3", "CPROB4",
                                "CPROB5","CPROB6","clslabel","patid")],by="patid")
a<-merge(a,classsize_unsup6584,by="clslabel")

colnames(a)[1]<-"modalcls"
lcapatid_train_unsup6584<-a  # Note: train, trainpat are no longer the correct indicators for this stratified split
rm(a,b)

# 85plus TRAIN #
colnames(results_train_unsup85plus)[43:44]<-c("clslabel","patid")
b<-subset(results_train_unsup85plus,!duplicated(results_train_unsup85plus$patid))
a<-merge(train_pat_mm_85plus,b[,c("CPROB1", "CPROB2", "CPROB3", "CPROB4",
                                "clslabel","patid")],by="patid")
a<-merge(a,classsize_unsup85plus,by="clslabel")

colnames(a)[1]<-"modalcls"
lcapatid_train_unsup85plus<-a  # Note: train, trainpat are no longer the correct indicators for this stratified split
rm(a,b)

# repeate merging results with main data process for test set #

# 18-44 TEST #
colnames(results_test_unsup1844)[44:45]<-c("clslabel","patid")
b<-subset(results_test_unsup1844,!duplicated(results_test_unsup1844$patid))
a<-merge(test_pat_mm_1844,b[,c("CPROB1", "CPROB2", "CPROB3", "CPROB4",
                                "CPROB5","clslabel","patid")],by="patid")
a<-merge(a,classsize_test_unsup1844,by="clslabel")

colnames(a)[1]<-"modalcls"
lcapatid_test_unsup1844<-a  # Note: test, testpat are no longer the correct indicators for this stratified split
rm(a,b)

# 45-64 TEST #
colnames(results_test_unsup4564)[44:45]<-c("clslabel","patid")
b<-subset(results_test_unsup4564,!duplicated(results_test_unsup4564$patid))
a<-merge(test_pat_mm_4564,b[,c("CPROB1", "CPROB2", "CPROB3", "CPROB4",
                                "CPROB5","clslabel","patid")],by="patid")
a<-merge(a,classsize_test_unsup4564,by="clslabel")

colnames(a)[1]<-"modalcls"
lcapatid_test_unsup4564<-a  # Note: test, testpat are no longer the correct indicators for this stratified split
rm(a,b)

# 65-84 TEST #
colnames(results_test_unsup6584)[45:46]<-c("clslabel","patid")
b<-subset(results_test_unsup6584,!duplicated(results_test_unsup6584$patid))
a<-merge(test_pat_mm_6584,b[,c("CPROB1", "CPROB2", "CPROB3", "CPROB4",
                                "CPROB5","CPROB6","clslabel","patid")],by="patid")
a<-merge(a,classsize_test_unsup6584,by="clslabel")

colnames(a)[1]<-"modalcls"
lcapatid_test_unsup6584<-a  # Note: test, testpat are no longer the correct indicators for this stratified split
rm(a,b)

# 85plus TEST #
colnames(results_test_unsup85plus)[42:43]<-c("clslabel","patid") # careful: here index change to 4243 as c_mig is all zero for test data
b<-subset(results_test_unsup85plus,!duplicated(results_test_unsup85plus$patid))
a<-merge(test_pat_mm_85plus,b[,c("CPROB1", "CPROB2", "CPROB3", "CPROB4",
                                  "clslabel","patid")],by="patid")
a<-merge(a,classsize_test_unsup85plus,by="clslabel")

colnames(a)[1]<-"modalcls"
lcapatid_test_unsup85plus<-a  # Note: test, testpat are no longer the correct indicators for this stratified split
rm(a,b)

# set up categorical predictors 
lcapatid_train_unsup1844<-fastDummies::dummy_cols(lcapatid_train_unsup1844,select_columns=c("gender","imd2015_5","modalcls", "smokstatus"))
lcapatid_train_unsup4564<-fastDummies::dummy_cols(lcapatid_train_unsup4564,select_columns=c("gender","imd2015_5","modalcls", "smokstatus"))
lcapatid_train_unsup6584<-fastDummies::dummy_cols(lcapatid_train_unsup6584,select_columns=c("gender","imd2015_5","modalcls", "smokstatus"))
lcapatid_train_unsup85plus<-fastDummies::dummy_cols(lcapatid_train_unsup85plus,select_columns=c("gender","imd2015_5","modalcls", "smokstatus"))

lcapatid_test_unsup1844<-fastDummies::dummy_cols(lcapatid_test_unsup1844,select_columns=c("gender","imd2015_5","modalcls", "smokstatus"))
lcapatid_test_unsup4564<-fastDummies::dummy_cols(lcapatid_test_unsup4564,select_columns=c("gender","imd2015_5","modalcls", "smokstatus"))
lcapatid_test_unsup6584<-fastDummies::dummy_cols(lcapatid_test_unsup6584,select_columns=c("gender","imd2015_5","modalcls", "smokstatus"))
lcapatid_test_unsup85plus<-fastDummies::dummy_cols(lcapatid_test_unsup85plus,select_columns=c("gender","imd2015_5","modalcls", "smokstatus"))

#############################################################################################################################################

# Plot GLM results for unsupervised data : training#

# add a new co-morbid column for all lcapatid_df #
# merge this file with ALL current results and data files # 
dflist1<-list(lcapatid_train_unsup1844,lcapatid_train_unsup4564,lcapatid_train_unsup6584,
             lcapatid_train_unsup85plus,
             lcapatid_test_unsup1844,lcapatid_test_unsup4564,lcapatid_test_unsup6584,
             lcapatid_test_unsup85plus)
name.dflist1<-c("lcapatid_train_unsup1844","lcapatid_train_unsup4564","lcapatid_train_unsup6584",
               "lcapatid_train_unsup85plus",
               "lcapatid_test_unsup1844","lcapatid_test_unsup4564","lcapatid_test_unsup6584",
               "lcapatid_test_unsup85plus")

dflist2<-list(classsize_unsup1844,classsize_unsup4564,classsize_unsup6584,
              classsize_unsup85plus,
              classsize_test_unsup1844,classsize_test_unsup4564,classsize_test_unsup6584,
              classsize_test_unsup85plus)
name.dflist2<-c("classsize_unsup1844","classsize_unsup4564","classsize_unsup6584",
                "classsize_unsup85plus",
                "classsize_test_unsup1844","classsize_test_unsup4564","classsize_test_unsup6584",
                "classsize_test_unsup85plus")


for (i in 1:length(dflist1)){
  print(i)
  dflist1[[i]]<- merge(dflist1[[i]],dflist2[[i]][,c('clslabel','pr(True)_3',"top3morb_3")],
                       by.x='modalcls',by.y="clslabel",all.x=T)
  assign(name.dflist1[i], dflist1[[i]])
}

#######

## 18-44 ##

# relabel predictors for meaningful plots ref. = clslabel3 for training (IBS), 2 (never smoke)=ref for smoking#
df<-lcapatid_train_unsup1844[,c("ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
                     "dur","d2013","d2016",
                     "modalcls_1","modalcls_2","modalcls_3","modalcls_4",
                     "modalcls_5",
                     "gender_2",
                     "imd2015_5_2","imd2015_5_3","imd2015_5_4","imd2015_5_5",
                     "smokstatus_1","smokstatus_3","bmi","age")]

# know the class labels #
lcapatid_train_unsup1844 %>% group_by(modalcls) %>% summarise (top3morb_3=paste(top3morb_3,collapse=""),clsprop=mean(clsprop))

# Careful, label is not ordered, need to match with class label
colnames(df)[7:20]<-c("AST,IBS,DEP","PSM,ALC,DEP", "PNC,HL,HYP","DEP,ANX,PNC",
                      "IBS,DEP,HL",
                      "Female",
                      "IMD_2","IMD_3","IMD_4","IMD_5",
                      "Current_smoker", "Ex_smoker","BMI","age"
                      )

# set variable label for response
set_label(df$ctconsult12mons_filterdur) <- "Consultations"
set_label(df$cthospspl_12mons) <- "Hospitalisation"
set_label(df$ctprod12mons4plus_BNF) <- "Repeat prescriptions"

# Note, need to leave out 1 ref. class IBS,DEP,HL.
# Here we orderthe covariates from largest to smallest in cluster size #

m1_1844<-glm.nb(ctconsult12mons_filterdur~., data = df[,c("ctconsult12mons_filterdur","DEP,ANX,PNC","PNC,HL,HYP",
                                                "AST,IBS,DEP","PSM,ALC,DEP","Female",
                                                 "IMD_2","IMD_3","IMD_4","IMD_5",
                                                "Current_smoker", "Ex_smoker","BMI","age")])
m2_1844<-glm.nb(cthospspl_12mons~., data = df[,c("cthospspl_12mons","DEP,ANX,PNC","PNC,HL,HYP",
                                            "AST,IBS,DEP","PSM,ALC,DEP","Female",
                                            "IMD_2","IMD_3","IMD_4","IMD_5",
                                            "Current_smoker", "Ex_smoker","BMI","age")])
m3_1844<-glm.nb(ctprod12mons4plus_BNF~., data = df[,c("ctprod12mons4plus_BNF","DEP,ANX,PNC","PNC,HL,HYP",
                                                 "AST,IBS,DEP","PSM,ALC,DEP","Female",
                                                 "IMD_2","IMD_3","IMD_4","IMD_5",
                                                 "Current_smoker", "Ex_smoker","BMI","age")])
m4_1844<-glm(d2013~., data = df[,c("d2013","DEP,ANX,PNC","PNC,HL,HYP",
                                           "AST,IBS,DEP","PSM,ALC,DEP","Female",
                                           "IMD_2","IMD_3","IMD_4","IMD_5",
                                           "Current_smoker", "Ex_smoker","BMI","age")])
m5_1844<-glm(d2016~., data = df[,c("d2016","DEP,ANX,PNC","PNC,HL,HYP",
                                           "AST,IBS,DEP","PSM,ALC,DEP","Female",
                                           "IMD_2","IMD_3","IMD_4","IMD_5",
                                           "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

#summarise GLM exp results #
tab_model(m1_1844,m2_1844,m3_1844,m4_1844,m5_1844,transform = "exp",collapse.ci = TRUE,
                 p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          dv.labels = c("GP consultation", "Hospitalisation","Repeat prescription",
                        "2-year mortality","5-year mortality"),
          file = "C:\\MM_proc_data\\Results\\expresults_train_patid_1844.html")

# summarise log-odds results #
results_train_patid_1844<-mtable("No. GP visits"=m1_1844,
                      "No. hospital spells"=m2_1844,
                      "No. repeat prescriptions"=m3_1844,
                      "Mortality 2-yr"=m4_1844,
                      "Mortality 5-yr"=m5_1844,
                      summary.stats = TRUE)
capture.output(results_train_patid_1844, file = "C:\\MM_proc_data\\Results\\results_train_patid_1844.txt")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_unsup_service_patid_1844.pdf")
a<-plot_coefs(m1_1844,m2_1844,m3_1844,scale = TRUE,
              model.names = c("# Consultation","# Hospitalisation spells",
                              "# repeat prescriptions"),
              legend.title = "Service use",ci_level = 0.95)
a
dev.off()

pdf(file="C:\\MM_proc_data\\Results\\glm_unsup_mort_patid_1844.pdf")
b<-plot_summs(m4_1844,m5_1844,scale = TRUE,
              model.names = c("2-year mortality",
                              "5-year mortality"),
              legend.title = "Mortality")
b
dev.off()

# relationship between clusters and demographics (cluster size big to small) #
c1_1844_profile<-glm(`DEP,ANX,PNC`~., data = df[,c("DEP,ANX,PNC","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                          "Current_smoker", "Ex_smoker","BMI","age")])
c2_1844_profile<-glm(`PNC,HL,HYP`~., data = df[,c("PNC,HL,HYP","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                     "Current_smoker", "Ex_smoker","BMI","age")])
c3_1844_profile<-glm(`AST,IBS,DEP`~., data = df[,c("AST,IBS,DEP","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c4_1844_profile<-glm(`IBS,DEP,HL`~., data = df[,c("IBS,DEP,HL","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
c5_1844_profile<-glm(`PSM,ALC,DEP`~., data = df[,c("PSM,ALC,DEP","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

tab_model(c1_1844_profile,c2_1844_profile,c3_1844_profile,c4_1844_profile,c5_1844_profile,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          file = "C:\\MM_proc_data\\Results\\results_train_patid_1844_profile.html")



# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_patid_1844_profile.pdf")
a<-plot_coefs(c1_1844_profile,c2_1844_profile,c3_1844_profile,c4_1844_profile,c5_1844_profile,scale = TRUE,
              model.names = c("DEP,ANX,PNC","PNC,HL,HYP", "AST,IBS,DEP","IBS,DEP,HL","PSM,ALC,DEP"),
              legend.title = "MM clusters",ci_level = 0.95)
a
dev.off()


## 45-64 ##
# relabel predictors for meaningful plots ref. = clslabel3 for training (IBS)#
df<-lcapatid_train_unsup4564[,c("ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
                                "dur","d2013","d2016",
                                "modalcls_1","modalcls_2","modalcls_3","modalcls_4",
                                "modalcls_5",
                                "gender_2",
                                "imd2015_5_2","imd2015_5_3","imd2015_5_4","imd2015_5_5",
                                "smokstatus_1", "smokstatus_3","bmi","age")]
# know the class labels #
lcapatid_train_unsup4564 %>% group_by(modalcls) %>% summarise (top3morb_3=paste(top3morb_3,collapse=""),clsprop=mean(clsprop))


# Careful, label is not ordered, need to match with class label
colnames(df)[7:20]<-c("DEP,PNC,ANX","AST,COPD,PNC", "HYP,DIA,PNC","IBS,HL,PNC",
                      "ALC,PSM,PNC",
                      "Female",
                      "IMD_2","IMD_3","IMD_4","IMD_5",
                      "Current_smoker", "Ex_smoker","BMI","age")

# set variable label for response
set_label(df$ctconsult12mons_filterdur) <- "Consultations"
set_label(df$cthospspl_12mons) <- "Hospitalisation"
set_label(df$ctprod12mons4plus_BNF) <- "Repeat prescriptions"


# Note, need to leave out 1 ref. class IBS,HL,PNC
# Here we orderthe covariates from largest to smallest in cluster size #

m1_4564<-glm.nb(ctconsult12mons_filterdur~., data = df[,c("ctconsult12mons_filterdur","HYP,DIA,PNC",
                                                "DEP,PNC,ANX","AST,COPD,PNC","ALC,PSM,PNC","Female",
                                                "IMD_2","IMD_3","IMD_4","IMD_5",
                                                "Current_smoker", "Ex_smoker","BMI","age")])
m2_4564<-glm.nb(cthospspl_12mons~., data = df[,c("cthospspl_12mons","HYP,DIA,PNC",
                                                 "DEP,PNC,ANX","AST,COPD,PNC","ALC,PSM,PNC","Female",
                                                 "IMD_2","IMD_3","IMD_4","IMD_5",
                                                 "Current_smoker", "Ex_smoker","BMI","age")])
m3_4564<-glm.nb(ctprod12mons4plus_BNF~., data = df[,c("ctprod12mons4plus_BNF","HYP,DIA,PNC",
                                                      "DEP,PNC,ANX","AST,COPD,PNC","ALC,PSM,PNC","Female",
                                                      "IMD_2","IMD_3","IMD_4","IMD_5",
                                                      "Current_smoker", "Ex_smoker","BMI","age")])
m4_4564<-glm(d2013~., data = df[,c("d2013","HYP,DIA,PNC",
                                   "DEP,PNC,ANX","AST,COPD,PNC","ALC,PSM,PNC","Female",
                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                   "Current_smoker", "Ex_smoker","BMI","age")])
m5_4564<-glm(d2016~., data = df[,c("d2016","HYP,DIA,PNC",
                                   "DEP,PNC,ANX","AST,COPD,PNC","ALC,PSM,PNC","Female",
                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                   "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))


#summarise GLM exp results #
tab_model(m1_4564,m2_4564,m3_4564,m4_4564,m5_4564,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          dv.labels = c("GP consultation", "Hospitalisation","Repeat prescription",
                        "2-year mortality","5-year mortality"),
          file = "C:\\MM_proc_data\\Results\\expresults_train_patid_4564.html")

# summarise log-odds results #

results_train_patid_4564<-mtable("No. GP visits"=m1_4564,
                                 "No. hospital spells"=m2_4564,
                                 "No. repeat prescriptions"=m3_4564,
                                 "Mortality 2-yr"=m4_4564,
                                 "Mortality 5-yr"=m5_4564,
                                 summary.stats = TRUE)
capture.output(results_train_patid_4564, file = "C:\\MM_proc_data\\Results\\results_train_patid_4564.txt")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_unsup_service_patid_4564.pdf")
a<-plot_coefs(m1_4564,m2_4564,m3_4564,scale = TRUE,
              model.names = c("# Consultation","# Hospitalisation spells",
                              "# repeat prescriptions"),
              legend.title = "Service use",ci_level = 0.95)
a
dev.off()

pdf(file="C:\\MM_proc_data\\Results\\glm_unsup_mort_patid_4564.pdf")
b<-plot_summs(m4_4564,m5_4564,scale = TRUE,
              model.names = c("2-year mortality",
                              "5-year mortality"),
              legend.title = "Mortality")
b
dev.off()

# relationship between clusters and demographics (order in terms of class size) #
c1_4564_profile<-glm(`HYP,DIA,PNC`~., data = df[,c("HYP,DIA,PNC",
                                                   "Female",
                                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c2_4564_profile<-glm(`IBS,HL,PNC`~., data = df[,c("IBS,HL,PNC","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
c3_4564_profile<-glm(`DEP,PNC,ANX`~., data = df[,c("DEP,PNC,ANX","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c4_4564_profile<-glm(`AST,COPD,PNC`~., data = df[,c("AST,COPD,PNC","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c5_4564_profile<-glm(`ALC,PSM,PNC`~., data = df[,c("ALC,PSM,PNC","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

tab_model(c1_4564_profile,c2_4564_profile,c3_4564_profile,c4_4564_profile,c5_4564_profile,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          file = "C:\\MM_proc_data\\Results\\results_train_patid_4564_profile.html")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_patid_4564_profile.pdf")
a<-plot_coefs(c1_4564_profile,c2_4564_profile,c3_4564_profile,c4_4564_profile,c5_4564_profile,scale = TRUE,
              model.names = c("HYP,DIA,PNC","IBS,HL,PNC", "DEP,PNC,ANX","AST,COPD,PNC","ALC,PSM,PNC"),
              legend.title = "MM clusters",ci_level = 0.95)
a
dev.off()



## 65-84 ##
# relabel predictors for meaningful plots ref.#
df<-lcapatid_train_unsup6584[,c("ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
                                "dur","d2013","d2016",
                                "modalcls_1","modalcls_2","modalcls_3","modalcls_4",
                                "modalcls_5","modalcls_6",
                                "gender_2",
                                "imd2015_5_2","imd2015_5_3","imd2015_5_4","imd2015_5_5",
                                "smokstatus_1","smokstatus_3","bmi","age")]
# know the class labels #
lcapatid_train_unsup6584 %>% group_by(modalcls) %>% summarise (top3morb_3=paste(top3morb_3,collapse=""),clsprop=mean(clsprop))



# Careful, label is not ordered, need to match with class label
colnames(df)[7:21]<-c("CHD,AF,DIA","PNC,CHD,DEP", "HL,PSD,IBS","COPD,AST,PNC",
                      "DEP,PNC,ANX","HYP,DIA,CKD",
                      "Female",
                      "IMD_2","IMD_3","IMD_4","IMD_5",
                      "Current_smoker", "Ex_smoker","BMI","age")

# set variable label for response
set_label(df$ctconsult12mons_filterdur) <- "Consultations"
set_label(df$cthospspl_12mons) <- "Hospitalisation"
set_label(df$ctprod12mons4plus_BNF) <- "Repeat prescriptions"

# Note, need to leave out 1 ref. class HL,PSD.IBS
# Here we orderthe covariates from largest to smallest in cluster size #

m1_6584<-glm.nb(ctconsult12mons_filterdur~., data = df[,c("ctconsult12mons_filterdur","HYP,DIA,CKD",
                                                "DEP,PNC,ANX","CHD,AF,DIA",
                                                "COPD,AST,PNC","PNC,CHD,DEP","Female",
                                                "IMD_2","IMD_3","IMD_4","IMD_5",
                                                "Current_smoker", "Ex_smoker","BMI","age")])
m2_6584<-glm.nb(cthospspl_12mons~., data = df[,c("cthospspl_12mons","HYP,DIA,CKD",
                                                 "DEP,PNC,ANX","CHD,AF,DIA",
                                                 "COPD,AST,PNC","PNC,CHD,DEP","Female",
                                                 "IMD_2","IMD_3","IMD_4","IMD_5",
                                                 "Current_smoker", "Ex_smoker","BMI","age")])
m3_6584<-glm.nb(ctprod12mons4plus_BNF~., data = df[,c("ctprod12mons4plus_BNF","HYP,DIA,CKD",
                                                      "DEP,PNC,ANX","CHD,AF,DIA",
                                                      "COPD,AST,PNC","PNC,CHD,DEP","Female",
                                                      "IMD_2","IMD_3","IMD_4","IMD_5",
                                                      "Current_smoker", "Ex_smoker","BMI","age")])
m4_6584<-glm(d2013~., data = df[,c("d2013","HYP,DIA,CKD",
                                   "DEP,PNC,ANX","CHD,AF,DIA",
                                   "COPD,AST,PNC","PNC,CHD,DEP","Female",
                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                   "Current_smoker", "Ex_smoker","BMI","age")])
m5_6584<-glm(d2016~., data = df[,c("d2016","HYP,DIA,CKD",
                                   "DEP,PNC,ANX","CHD,AF,DIA",
                                   "COPD,AST,PNC","PNC,CHD,DEP","Female",
                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                   "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

#summarise GLM exp results #
tab_model(m1_6584,m2_6584,m3_6584,m4_6584,m5_6584,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          dv.labels = c("GP consultation", "Hospitalisation","Repeat prescription",
                        "2-year mortality","5-year mortality"),
          file = "C:\\MM_proc_data\\Results\\expresults_train_patid_6584.html")

#summarise GLM log-odds results #

results_train_patid_6584<-mtable("No. GP visits"=m1_6584,
                                 "No. hospital spells"=m2_6584,
                                 "No. repeat prescription"=m3_6584,
                                 "Mortality 2-yr"=m4_6584,
                                 "Mortality 5-yr"=m5_6584,
                                 summary.stats = TRUE)
capture.output(results_train_patid_6584, file = "C:\\MM_proc_data\\Results\\results_train_patid_6584.txt")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_unsup_service_patid_6584.pdf")
a<-plot_coefs(m1_6584,m2_6584,m3_6584,scale = TRUE,
              model.names = c("# Consultation","# Hospitalisation spells",
                              "# Repeat prescriptions"),
              legend.title = "Service use",ci_level = 0.95)
a
dev.off()

pdf(file="C:\\MM_proc_data\\Results\\glm_unsup_mort_patid_6584.pdf")
b<-plot_summs(m4_6584,m5_6584,scale = TRUE,
              model.names = c("2-year mortality",
                              "5-year mortality"),
              legend.title = "Mortality")
b
dev.off()


# relationship between clusters and demographics (order in terms of class size) #
c1_6584_profile<-glm(`HYP,DIA,CKD`~., data = df[,c("HYP,DIA,CKD",
                                                   "Female",
                                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c2_6584_profile<-glm(`HL,PSD,IBS`~., data = df[,c("HL,PSD,IBS","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
c3_6584_profile<-glm(`DEP,PNC,ANX`~., data = df[,c("DEP,PNC,ANX","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c4_6584_profile<-glm(`CHD,AF,DIA`~., data = df[,c("CHD,AF,DIA","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c5_6584_profile<-glm(`COPD,AST,PNC`~., data = df[,c("COPD,AST,PNC","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c6_6584_profile<-glm(`PNC,CHD,DEP`~., data = df[,c("PNC,CHD,DEP","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                    "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))


tab_model(c1_6584_profile,c2_6584_profile,c3_6584_profile,c4_6584_profile,c5_6584_profile,c6_6584_profile,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          file = "C:\\MM_proc_data\\Results\\results_train_patid_6584_profile.html")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_patid_6584_profile.pdf")
a<-plot_coefs(c1_6584_profile,c2_6584_profile,c3_6584_profile,c4_6584_profile,c5_6584_profile,c6_6584_profile,scale = TRUE,
              model.names = c("HYP,DIA,CKD","HL,PSD,IBS", "DEP,PNC,ANX","CHD,AF,DIA","COPD,AST,PNC","PNC,CHD,DEP"),
              legend.title = "MM clusters",ci_level = 0.95)
a
dev.off()

## 85plus ## 
# relabel predictors for meaningful plots ref.#
df<-lcapatid_train_unsup85plus[,c("ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
                                "dur","d2013","d2016",
                                "modalcls_1","modalcls_2","modalcls_3","modalcls_4",
                                "gender_2",
                                "imd2015_5_2","imd2015_5_3","imd2015_5_4","imd2015_5_5",
                                "smokstatus_1","smokstatus_3","bmi","age")]
# know the class labels #
lcapatid_train_unsup85plus %>% group_by(modalcls) %>% summarise (top3morb_3=paste(top3morb_3,collapse=""),clsprop=mean(clsprop))

# Careful, label is not ordered, need to match with class label
colnames(df)[7:19]<-c("HF,CHD,AF","AST,COPD,PNC", "HYP,HL,DIA","PNC,DEP,CSP",
                      "Female",
                      "IMD_2","IMD_3","IMD_4","IMD_5",
                      "Current_smoker", "Ex_smoker","BMI","age")

# set variable label for response
set_label(df$ctconsult12mons_filterdur) <- "Consultations"
set_label(df$cthospspl_12mons) <- "Hospitalisation"
set_label(df$ctprod12mons4plus_BNF) <- "Repeat prescriptions"

# Note, need to leave out 1 ref. class HYP,HL,DIA
# Here we orderthe covariates from largest to smallest in cluster size #

m1_85plus<-glm.nb(ctconsult12mons_filterdur~., data = df[,c("ctconsult12mons_filterdur",
                                                  "PNC,DEP,CSP","HF,CHD,AF",
                                                  "AST,COPD,PNC","Female",
                                                  "IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
m2_85plus<-glm.nb(cthospspl_12mons~., data = df[,c("cthospspl_12mons",
                                                   "PNC,DEP,CSP","HF,CHD,AF",
                                                   "AST,COPD,PNC","Female",
                                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
m3_85plus<-glm.nb(ctprod12mons4plus_BNF~., data = df[,c("ctprod12mons4plus_BNF",
                                                        "PNC,DEP,CSP","HF,CHD,AF",
                                                        "AST,COPD,PNC","Female",
                                                        "IMD_2","IMD_3","IMD_4","IMD_5",
                                                        "Current_smoker", "Ex_smoker","BMI","age")])
m4_85plus<-glm(d2013~., data = df[,c("d2013",
                                     "PNC,DEP,CSP","HF,CHD,AF",
                                     "AST,COPD,PNC","Female",
                                     "IMD_2","IMD_3","IMD_4","IMD_5",
                                     "Current_smoker", "Ex_smoker","BMI","age")])
m5_85plus<-glm(d2016~., data = df[,c("d2016",
                                     "PNC,DEP,CSP","HF,CHD,AF",
                                     "AST,COPD,PNC","Female",
                                     "IMD_2","IMD_3","IMD_4","IMD_5",
                                     "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
## rm(list = ls(pattern = tbresults_))

#summarise GLM exp results #
tab_model(m1_85plus,m2_85plus,m3_85plus,m4_85plus,m5_85plus,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          dv.labels = c("GP consultation", "Hospitalisation","Repeat prescription",
                        "2-year mortality","5-year mortality"),
          file = "C:\\MM_proc_data\\Results\\expresults_train_patid_85plus.html")

#summarise GLM log-odds results #


results_train_patid_85plus<-mtable("No. GP visits"=m1_85plus,
                                 "No. hospital spells"=m2_85plus,
                                 "No. repeat prescriptions"=m3_85plus,
                                 "Mortality 2-yr"=m4_85plus,
                                 "Mortality 5-yr"=m5_85plus,
                                 summary.stats = TRUE)
capture.output(results_train_patid_85plus, file = "C:\\MM_proc_data\\Results\\results_train_patid_85plus.txt")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_unsup_service_patid_85plus.pdf")
a<-plot_coefs(m1_85plus,m2_85plus,m3_85plus,scale = TRUE,
              model.names = c("# Consultation","# Hospitalisation spells",
                              "# Repeat prescriptions"),
              legend.title = "Service use",ci_level = 0.95)
a
dev.off()

pdf(file="C:\\MM_proc_data\\Results\\glm_unsup_mort_patid_85plus.pdf")
b<-plot_summs(m4_85plus,m5_85plus,scale = TRUE,
              model.names = c("2-year mortality",
                              "5-year mortality"),
              legend.title = "Mortality")
b
dev.off()

# relationship between clusters and demographics (order in terms of class size) #
c1_85plus_profile<-glm(`HYP,HL,DIA`~., data = df[,c("HYP,HL,DIA",
                                                     "Female",
                                                     "IMD_2","IMD_3","IMD_4","IMD_5",
                                                     "Current_smoker", "Ex_smoker","BMI","age")])
c2_85plus_profile<-glm(`PNC,DEP,CSP`~., data = df[,c("PNC,DEP,CSP","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
c3_85plus_profile<-glm(`HF,CHD,AF`~., data = df[,c("HF,CHD,AF","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c4_85plus_profile<-glm(`AST,COPD,PNC`~., data = df[,c("AST,COPD,PNC","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])


# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

tab_model(c1_85plus_profile,c2_85plus_profile,c3_85plus_profile,c4_85plus_profile,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          file = "C:\\MM_proc_data\\Results\\results_train_patid_85plus_profile.html")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_patid_85plus_profile.pdf")
a<-plot_coefs(c1_85plus_profile,c2_85plus_profile,c3_85plus_profile,c4_85plus_profile,scale = TRUE,
              model.names = c("HYP,HL,DIA","PNC,DEP,CSP", "HF,CHD,AF","AST,COPD,PNC"),
              legend.title = "MM clusters",ci_level = 0.95)
a
dev.off()

save.image(file='C:\\MM_proc_data\\final_patid_strata_final_filterdur.RData')                                                                                                                                                                    

#############################################################################################################################
# TEST data # 
# Show GLM results for unsupervised data (test)#

## 18-44 ##
# relabel predictors for meaningful plots ref. = clslabel3 for testing (IBS)#
df<-lcapatid_test_unsup1844[,c("ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
                                "dur","d2013","d2016",
                                "modalcls_1","modalcls_2","modalcls_3","modalcls_4",
                                "modalcls_5",
                                "gender_2",
                                "imd2015_5_2","imd2015_5_3","imd2015_5_4","imd2015_5_5",
                               "smokstatus_1","smokstatus_3","bmi","age")]

# know the class labels #
lcapatid_test_unsup1844 %>% group_by(modalcls) %>% summarise (top3morb_3=paste(top3morb_3,collapse=""),clsprop=mean(clsprop))


# Careful, label is not ordered, need to match with class label
colnames(df)[7:20]<-c("PSM,ALC,PNC","IBS,DEP,AST", "HYP,DIA,DEP","AST,HL,PNC",
                      "DEP,ANX,PNC",
                      "Female",
                      "IMD_2","IMD_3","IMD_4","IMD_5",
                      "Current_smoker", "Ex_smoker","BMI","age"
                        )

# set variable label for response
set_label(df$ctconsult12mons_filterdur) <- "Consultations"
set_label(df$cthospspl_12mons) <- "Hospitalisation"
set_label(df$ctprod12mons4plus_BNF) <- "Repeat prescription"

# Note, need to leave out 1 ref. class IBS,DEP,AST.
# Here we orderthe covariates from largest to smallest in cluster size #

m1_1844<-glm.nb(ctconsult12mons_filterdur~., data = df[,c("ctconsult12mons_filterdur","DEP,ANX,PNC","AST,HL,PNC",
                                                "HYP,DIA,DEP","PSM,ALC,PNC","Female",
                                                "IMD_2","IMD_3","IMD_4","IMD_5",
                                                "Current_smoker", "Ex_smoker","BMI","age")])
m2_1844<-glm.nb(cthospspl_12mons~., data = df[,c("cthospspl_12mons","DEP,ANX,PNC","AST,HL,PNC",
                                                 "HYP,DIA,DEP","PSM,ALC,PNC","Female",
                                                 "IMD_2","IMD_3","IMD_4","IMD_5",
                                                 "Current_smoker", "Ex_smoker","BMI","age")])
m3_1844<-glm.nb(ctprod12mons4plus_BNF~., data = df[,c("ctprod12mons4plus_BNF","DEP,ANX,PNC","AST,HL,PNC",
                                                      "HYP,DIA,DEP","PSM,ALC,PNC","Female",
                                                      "IMD_2","IMD_3","IMD_4","IMD_5",
                                                      "Current_smoker", "Ex_smoker","BMI","age")])
m4_1844<-glm(d2013~., data = df[,c("d2013","DEP,ANX,PNC","AST,HL,PNC",
                                   "HYP,DIA,DEP","PSM,ALC,PNC","Female",
                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                   "Current_smoker", "Ex_smoker","BMI","age")])
m5_1844<-glm(d2016~., data = df[,c("d2016","DEP,ANX,PNC","AST,HL,PNC",
                                   "HYP,DIA,DEP","PSM,ALC,PNC","Female",
                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                   "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
# rm(list = ls(pattern = "tbresults_"))
# summarise GLM exp results #
tab_model(m1_1844,m2_1844,m3_1844,m4_1844,m5_1844,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          dv.labels = c("GP consultation", "Hospitalisation","Repeat prescription",
                        "2-year mortality","5-year mortality"),
          file = "C:\\MM_proc_data\\Results\\expresults_test_patid_1844.html")

#summarise GLM log-odds results #

results_test_patid_1844<-mtable("No. GP visits"=m1_1844,
                                 "No. hospital spells"=m2_1844,
                                 "No. repeat prescriptions"=m3_1844,
                                 "Mortality 2-yr"=m4_1844,
                                 "Mortality 5-yr"=m5_1844,
                                 summary.stats = TRUE)
capture.output(results_test_patid_1844, file = "C:\\MM_proc_data\\Results\\results_test_patid_1844.txt")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_unsuptest_service_patid_1844.pdf")
a<-plot_coefs(m1_1844,m2_1844,m3_1844,scale = TRUE,
              model.names = c("# Consultation","# Hospitalisation spells",
                              "# Repeat prescriptions"),
              legend.title = "Service use",ci_level = 0.95)
a
dev.off()

pdf(file="C:\\MM_proc_data\\Results\\glm_unsuptest_mort_patid_1844.pdf")
b<-plot_summs(m4_1844,m5_1844,scale = TRUE,
              model.names = c("2-year mortality",
                              "5-year mortality"),
              legend.title = "Mortality")
b
dev.off()

# relationship between clusters and demographics (order in terms of class size) #
c1_1844_profile<-glm(`DEP,ANX,PNC`~., data = df[,c("DEP,ANX,PNC",
                                                    "Female",
                                                    "IMD_2","IMD_3","IMD_4","IMD_5",
                                                    "Current_smoker", "Ex_smoker","BMI","age")])
c2_1844_profile<-glm(`AST,HL,PNC`~., data = df[,c("AST,HL,PNC","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                     "Current_smoker", "Ex_smoker","BMI","age")])
c3_1844_profile<-glm(`IBS,DEP,AST`~., data = df[,c("IBS,DEP,AST","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                     "Current_smoker", "Ex_smoker","BMI","age")])
c4_1844_profile<-glm(`HYP,DIA,DEP`~., data = df[,c("HYP,DIA,DEP","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                      "Current_smoker", "Ex_smoker","BMI","age")])
c5_1844_profile<-glm(`PSM,ALC,PNC`~., data = df[,c("PSM,ALC,PNC","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

tab_model(c1_1844_profile,c2_1844_profile,c3_1844_profile,c4_1844_profile,c5_1844_profile,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          file = "C:\\MM_proc_data\\Results\\results_test_patid_1844_profile.html")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_test_patid_1844_profile.pdf")
a<-plot_coefs(c1_1844_profile,c2_1844_profile,c3_1844_profile,c4_1844_profile,c5_1844_profile,scale = TRUE,
              model.names = c("DEP,ANX,PNC","AST,HL,PNC", "IBS,DE,AST","HYP,DIA,DEP","PSM,ALC,PNC"),
              legend.title = "MM clusters",ci_level = 0.95)
a
dev.off()


## 45-64 ##
df<-lcapatid_test_unsup4564[,c("ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
                                "dur","d2013","d2016",
                                "modalcls_1","modalcls_2","modalcls_3","modalcls_4",
                                "modalcls_5",
                                "gender_2",
                                "imd2015_5_2","imd2015_5_3","imd2015_5_4","imd2015_5_5",
                               "smokstatus_1","smokstatus_3","bmi","age")]
# know the class labels #
lcapatid_test_unsup4564 %>% group_by(modalcls) %>% summarise (top3morb_3=paste(top3morb_3,collapse=""),clsprop=mean(clsprop))

# Careful, label is not ordered, need to match with class label
colnames(df)[7:20]<-c("IBS,HL,AST","DEP,PNC,ANX", "ALC,PSM,PNC","PNC,CHD,DIA",
                      "HYP,DIA,CHD",
                      "Female",
                      "IMD_2","IMD_3","IMD_4","IMD_5",
                      "Current_smoker", "Ex_smoker","BMI","age")

# set variable label for response
set_label(df$ctconsult12mons_filterdur) <- "Consultations"
set_label(df$cthospspl_12mons) <- "Hospitalisation"
set_label(df$ctprod12mons4plus_BNF) <- "Repeat prescription"

# Note, need to leave out 1 ref. class IBS,HL,AST
# Here we orderthe covariates from largest to smallest in cluster size #

m1_4564<-glm.nb(ctconsult12mons_filterdur~., data = df[,c("ctconsult12mons_filterdur","HYP,DIA,CHD",
                                                "DEP,PNC,ANX","PNC,CHD,DIA","ALC,PSM,PNC","Female",
                                                "IMD_2","IMD_3","IMD_4","IMD_5",
                                                "Current_smoker", "Ex_smoker","BMI","age")])
m2_4564<-glm.nb(cthospspl_12mons~., data = df[,c("cthospspl_12mons","HYP,DIA,CHD",
                                                 "DEP,PNC,ANX","PNC,CHD,DIA","ALC,PSM,PNC","Female",
                                                 "IMD_2","IMD_3","IMD_4","IMD_5",
                                                 "Current_smoker", "Ex_smoker","BMI","age")])
m3_4564<-glm.nb(ctprod12mons4plus_BNF~., data = df[,c("ctprod12mons4plus_BNF","HYP,DIA,CHD",
                                                      "DEP,PNC,ANX","PNC,CHD,DIA","ALC,PSM,PNC","Female",
                                                      "IMD_2","IMD_3","IMD_4","IMD_5",
                                                      "Current_smoker", "Ex_smoker","BMI","age")])
m4_4564<-glm(d2013~., data = df[,c("d2013","HYP,DIA,CHD",
                                   "DEP,PNC,ANX","PNC,CHD,DIA","ALC,PSM,PNC","Female",
                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                   "Current_smoker", "Ex_smoker","BMI","age")])
m5_4564<-glm(d2016~., data = df[,c("d2016","HYP,DIA,CHD",
                                   "DEP,PNC,ANX","PNC,CHD,DIA","ALC,PSM,PNC","Female",
                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                   "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

#summarise GLM exp results #
tab_model(m1_4564,m2_4564,m3_4564,m4_4564,m5_4564,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          dv.labels = c("GP consultation", "Hospitalisation","Repeat prescription",
                        "2-year mortality","5-year mortality"),
          file = "C:\\MM_proc_data\\Results\\expresults_test_patid_4564.html")

#summarise GLM log-odds results #

results_test_patid_4564<-mtable("No. GP visits"=m1_4564,
                                 "No. hospital spells"=m2_4564,
                                 "No. repeat prescriptions"=m3_4564,
                                 "Mortality 2-yr"=m4_4564,
                                 "Mortality 5-yr"=m5_4564,
                                 summary.stats = TRUE)
capture.output(results_test_patid_4564, file = "C:\\MM_proc_data\\Results\\results_test_patid_4564.txt")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_unsuptest_service_patid_4564.pdf")
a<-plot_coefs(m1_4564,m2_4564,m3_4564,scale = TRUE,
              model.names = c("# Consultation","# Hospitalisation spells",
                              "# Repeat prescriptions"),
              legend.title = "Service use",ci_level = 0.95)
a
dev.off()

pdf(file="C:\\MM_proc_data\\Results\\glm_unsuptest_mort_patid_4564.pdf")
b<-plot_summs(m4_4564,m5_4564,scale = TRUE,
              model.names = c("2-year mortality",
                              "5-year mortality"),
              legend.title = "Mortality")
b
dev.off()

# relationship between clusters and demographics (order in terms of class size) #
c1_4564_profile<-glm(`HYP,DIA,CHD`~., data = df[,c("HYP,DIA,CHD",
                                                   "Female",
                                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c2_4564_profile<-glm(`IBS,HL,AST`~., data = df[,c("IBS,HL,AST","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
c3_4564_profile<-glm(`DEP,PNC,ANX`~., data = df[,c("DEP,PNC,ANX","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
c4_4564_profile<-glm(`PNC,CHD,DIA`~., data = df[,c("PNC,CHD,DIA","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c5_4564_profile<-glm(`ALC,PSM,PNC`~., data = df[,c("ALC,PSM,PNC","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))
tab_model(c1_4564_profile,c2_4564_profile,c3_4564_profile,c4_4564_profile,c5_4564_profile,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          file = "C:\\MM_proc_data\\Results\\results_test_patid_4564_profile.html")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_test_patid_4564_profile.pdf")
a<-plot_coefs(c1_4564_profile,c2_4564_profile,c3_4564_profile,c4_4564_profile,c5_4564_profile,scale = TRUE,
              model.names = c("HYP,DIA,CHD","IBS,HL,AST", "DEP,PNC,ANX","PNC,CHD,DIA","ALC,PSM,PNC"),
              legend.title = "MM clusters",ci_level = 0.95)
a
dev.off()


## 65-84 ##
# relabel predictors for meaningful plots ref.#
df<-lcapatid_test_unsup6584[,c("ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
                                "dur","d2013","d2016",
                                "modalcls_1","modalcls_2","modalcls_3","modalcls_4",
                                "modalcls_5","modalcls_6",
                                "gender_2",
                                "imd2015_5_2","imd2015_5_3","imd2015_5_4","imd2015_5_5",
                               "smokstatus_1","smokstatus_3","bmi","age")]

# know the class labels #
lcapatid_test_unsup6584 %>% group_by(modalcls) %>% summarise (top3morb_3=paste(top3morb_3,collapse=""),clsprop=mean(clsprop))

# Careful, label is not ordered, need to match with class label
colnames(df)[7:21]<-c("AST,COPD,PNC","DEP,PNC,ANX", "PNC,DEP,CHD","HL,IBS,DIV",
                      "CHD,AF,DIA","HYP,DIA,CKD",
                      "Female",
                      "IMD_2","IMD_3","IMD_4","IMD_5",
                      "Current_smoker", "Ex_smoker","BMI","age")

# set variable label for response
set_label(df$ctconsult12mons_filterdur) <- "Consultations"
set_label(df$cthospspl_12mons) <- "Hospitalisation"
set_label(df$ctprod12mons4plus_BNF) <- "Repeat prescription"

# Note, need to leave out 1 ref. class HL,IBS,DIV
# Here we orderthe covariates from largest to smallest in cluster size #

m1_6584<-glm.nb(ctconsult12mons_filterdur~., data = df[,c("ctconsult12mons_filterdur","HYP,DIA,CKD",
                                                "CHD,AF,DIA",
                                                "DEP,PNC,ANX","AST,COPD,PNC","PNC,DEP,CHD","Female",
                                                "IMD_2","IMD_3","IMD_4","IMD_5",
                                                "Current_smoker", "Ex_smoker","BMI","age")])
m2_6584<-glm.nb(cthospspl_12mons~., data = df[,c("cthospspl_12mons","HYP,DIA,CKD",
                                                 "CHD,AF,DIA",
                                                 "DEP,PNC,ANX","AST,COPD,PNC","PNC,DEP,CHD","Female",
                                                 "IMD_2","IMD_3","IMD_4","IMD_5",
                                                 "Current_smoker", "Ex_smoker","BMI","age")])
m3_6584<-glm.nb(ctprod12mons4plus_BNF~., data = df[,c("ctprod12mons4plus_BNF","HYP,DIA,CKD",
                                                      "CHD,AF,DIA",
                                                      "DEP,PNC,ANX","AST,COPD,PNC","PNC,DEP,CHD","Female",
                                                      "IMD_2","IMD_3","IMD_4","IMD_5",
                                                      "Current_smoker", "Ex_smoker","BMI","age")])
m4_6584<-glm(d2013~., data = df[,c("d2013","HYP,DIA,CKD",
                                   "CHD,AF,DIA",
                                   "DEP,PNC,ANX","AST,COPD,PNC","PNC,DEP,CHD","Female",
                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                   "Current_smoker", "Ex_smoker","BMI","age")])
m5_6584<-glm(d2016~., data = df[,c("d2016","HYP,DIA,CKD",
                                   "CHD,AF,DIA",
                                   "DEP,PNC,ANX","AST,COPD,PNC","PNC,DEP,CHD","Female",
                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                   "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

#summarise GLM exp results #
tab_model(m1_6584,m2_6584,m3_6584,m4_6584,m5_6584,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          dv.labels = c("GP consultation", "Hospitalisation","Repeat prescription",
                        "2-year mortality","5-year mortality"),
          file = "C:\\MM_proc_data\\Results\\expresults_test_patid_6584.html")

#summarise GLM log-odds results #
results_test_patid_6584<-mtable("No. GP visits"=m1_6584,
                                 "No. hospital spells"=m2_6584,
                                 "No. repeat prescriptions"=m3_6584,
                                 "Mortality 2-yr"=m4_6584,
                                 "Mortality 5-yr"=m5_6584,
                                 summary.stats = TRUE)
capture.output(results_test_patid_6584, file = "C:\\MM_proc_data\\Results\\results_test_patid_6584.txt")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_unsuptest_service_patid_6584.pdf")
a<-plot_coefs(m1_6584,m2_6584,m3_6584,scale = TRUE,
              model.names = c("# Consultation","# Hospitalisation spells",
                              "# Repeat prescriptions"),
              legend.title = "Service use",ci_level = 0.95)
a
dev.off()

pdf(file="C:\\MM_proc_data\\Results\\glm_unsuptest_mort_patid_6584.pdf")
b<-plot_summs(m4_6584,m5_6584,scale = TRUE,
              model.names = c("2-year mortality",
                              "5-year mortality"),
              legend.title = "Mortality")
b
dev.off()

# relationship between clusters and demographics (order in terms of class size) #
c1_6584_profile<-glm(`HYP,DIA,CKD`~., data = df[,c("HYP,DIA,CKD",
                                                   "Female",
                                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c2_6584_profile<-glm(`HL,IBS,DIV`~., data = df[,c("HL,IBS,DIV","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
c3_6584_profile<-glm(`CHD,AF,DIA`~., data = df[,c("CHD,AF,DIA","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
c4_6584_profile<-glm(`DEP,PNC,ANX`~., data = df[,c("DEP,PNC,ANX","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c5_6584_profile<-glm(`AST,COPD,PNC`~., data = df[,c("AST,COPD,PNC","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
c6_6584_profile<-glm(`PNC,DEP,CHD`~., data = df[,c("PNC,DEP,CHD","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                    "Current_smoker", "Ex_smoker","BMI","age")])
# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

tab_model(c1_6584_profile,c2_6584_profile,c3_6584_profile,c4_6584_profile,c5_6584_profile,c6_6584_profile,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          file = "C:\\MM_proc_data\\Results\\results_test_patid_6584_profile.html")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_test_patid_6584_profile.pdf")
a<-plot_coefs(c1_6584_profile,c2_6584_profile,c3_6584_profile,c4_6584_profile,c5_6584_profile,c6_6584_profile,scale = TRUE,
              model.names = c("HYP,DIA,CKD","HL,IBS,DIV", "CHD,AF,DIA","DEP,PNC,ANX",
                              "AST,COPD,PNC","PNC,DEP,CHD"),
              legend.title = "MM clusters",ci_level = 0.95)
a
dev.off()

## 85plus ##
# relabel predictors for meaningful plots ref.#
df<-lcapatid_test_unsup85plus[,c("ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
                                  "dur","d2013","d2016",
                                  "modalcls_1","modalcls_2","modalcls_3","modalcls_4",
                                  "gender_2",
                                  "imd2015_5_2","imd2015_5_3","imd2015_5_4","imd2015_5_5",
                                 "smokstatus_1","smokstatus_3","bmi","age")]

lcapatid_test_unsup85plus %>% group_by(modalcls) %>% summarise (top3morb_3=paste(top3morb_3,collapse=""),clsprop=mean(clsprop))

# Careful, label is not ordered, need to match with class label
colnames(df)[7:19]<-c("PNC,DEP,CSP","COPD,AST,PSD", "DEP,PNC,CSP","HYP,CHD,HL",
                      "Female",
                      "IMD_2","IMD_3","IMD_4","IMD_5",
                      "Current_smoker", "Ex_smoker","BMI","age")

# set variable label for response
set_label(df$ctconsult12mons_filterdur) <- "Consultations"
set_label(df$cthospspl_12mons) <- "Hospitalisation"
set_label(df$ctprod12mons4plus_BNF) <- "Repeat prescription"

# Note, need to leave out 1 ref. class HYP,CHD,HL
# Here we orderthe covariates from largest to smallest in cluster size #

m1_85plus<-glm.nb(ctconsult12mons_filterdur~., data = df[,c("ctconsult12mons_filterdur",
                                                  "DEP,PNC,CSP",
                                                  "PNC,DEP,CSP","COPD,AST,PSD","Female",
                                                  "IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
m2_85plus<-glm.nb(cthospspl_12mons~., data = df[,c("cthospspl_12mons",
                                                   "DEP,PNC,CSP",
                                                   "PNC,DEP,CSP","COPD,AST,PSD","Female",
                                                   "IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])
m3_85plus<-glm.nb(ctprod12mons4plus_BNF~., data = df[,c("ctprod12mons4plus_BNF",
                                                        "DEP,PNC,CSP",
                                                        "PNC,DEP,CSP","COPD,AST,PSD","Female",
                                                        "IMD_2","IMD_3","IMD_4","IMD_5",
                                                        "Current_smoker", "Ex_smoker","BMI","age")])
m4_85plus<-glm(d2013~., data = df[,c("d2013",
                                     "DEP,PNC,CSP",
                                     "PNC,DEP,CSP","COPD,AST,PSD","Female",
                                     "IMD_2","IMD_3","IMD_4","IMD_5",
                                     "Current_smoker", "Ex_smoker","BMI","age")])
m5_85plus<-glm(d2016~., data = df[,c("d2016",
                                     "DEP,PNC,CSP",
                                     "PNC,DEP,CSP","COPD,AST,PSD","Female",
                                     "IMD_2","IMD_3","IMD_4","IMD_5",
                                     "Current_smoker", "Ex_smoker","BMI","age")])

# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

#summarise GLM exp results #
tab_model(m1_85plus,m2_85plus,m3_85plus,m4_85plus,m5_85plus,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          dv.labels = c("GP consultation", "Hospitalisation","Repeat prescription",
                        "2-year mortality","5-year mortality"),
          file = "C:\\MM_proc_data\\Results\\expresults_test_patid_85plus.html")

#summarise GLM log-odds results #

results_test_patid_85plus<-mtable("No. GP visits"=m1_85plus,
                                   "No. hospital spells"=m2_85plus,
                                   "No. repeat prescription"=m3_85plus,
                                   "Mortality 2-yr"=m4_85plus,
                                   "Mortality 5-yr"=m5_85plus,
                                   summary.stats = TRUE)
capture.output(results_test_patid_85plus, file = "C:\\MM_proc_data\\Results\\results_test_patid_85plus.txt")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_unsuptest_service_patid_85plus.pdf")
a<-plot_coefs(m1_85plus,m2_85plus,m3_85plus,scale = TRUE,
              model.names = c("# Consultation","# Hospitalisation spells",
                              "# Repeat prescriptions"),
              legend.title = "Service use",ci_level = 0.95)
a
dev.off()

pdf(file="C:\\MM_proc_data\\Results\\glm_unsuptest_mort_patid_85plus.pdf")
b<-plot_summs(m4_85plus,m5_85plus,scale = TRUE,
              model.names = c("2-year mortality",
                              "5-year mortality"),
              legend.title = "Mortality")
b
dev.off()

# relationship between clusters and demographics (order in terms of class size) #
c1_85plus_profile<-glm(`HYP,CHD,HL`~., data = df[,c("HYP,CHD,HL",
                                                     "Female",
                                                     "IMD_2","IMD_3","IMD_4","IMD_5",
                                                     "Current_smoker", "Ex_smoker","BMI","age")])
c2_85plus_profile<-glm(`DEP,PNC,CSP`~., data = df[,c("DEP,PNC,CSP","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
c3_85plus_profile<-glm(`PNC,DEP,CSP`~., data = df[,c("PNC,DEP,CSP","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                  "Current_smoker", "Ex_smoker","BMI","age")])
c4_85plus_profile<-glm(`COPD,AST,PSD`~., data = df[,c("COPD,AST,PSD","Female","IMD_2","IMD_3","IMD_4","IMD_5",
                                                   "Current_smoker", "Ex_smoker","BMI","age")])


# Pull results together into 1 table  
#rm(list = ls(pattern = "tbresults_"))

tab_model(c1_85plus_profile,c2_85plus_profile,c3_85plus_profile,c4_85plus_profile,transform = "exp",collapse.ci = TRUE,
          p.style = "a", p.threshold = c(0.1, 0.05,0.01),
          file = "C:\\MM_proc_data\\Results\\results_test_patid_85plus_profile.html")

# graph all estimated coefficients #
pdf(file="C:\\MM_proc_data\\Results\\glm_test_patid_85plus_profile.pdf")
a<-plot_coefs(c1_85plus_profile,c2_85plus_profile,c3_85plus_profile,c4_85plus_profile,scale = TRUE,
              model.names = c("HYP,CHD,HL","DEP,PNC,CSP", "PNC,DEP,CSP","COPD,AST,PSD"),
              legend.title = "MM clusters",ci_level = 0.95)
a
dev.off()
#####################################################################################################

# Clusterd heatmap (now heatmap is in excel)#

# Unsupervised : training# 

# 18-44 #

a<-lcmodel_mm_train_unsup1844
a<-merge(a,classsize_unsup1844[,c("top3morb_3","clslabel")],by="clslabel")

a$Var2<-paste(round(100*a$clsprop,2),"%",";",a$top3morb_3)

# convert data from long to wide
a<-a[order(-a$clsprop,-a$value),] # order by highest to lowest on class size and conditional prob.
a<- reshape(a[,c("Morbidityfull","value","Var2")], 
            timevar = "Var2",
            idvar = c("Morbidityfull"),
            direction = "wide")
a2 <- a[,-1]
rownames(a2) <- a[,1]
a<-a2
names(a) <- substring(names(a), 7) # now a's colnames is very long
heatdata_unsup<-a # store map info into heatdata; plotted graph is based on "a"
colnames(a) <- paste('C',1:5,sep='') # this is order from highest to lowest on class size

# find extra characteristics for each cluster #
c<-lcapatid_train_unsup1844

  ## smokestatus final order: 1= current, 2= never, 3 = ex

extra <-c[,c("modalcls","clsprop","imd2015_5","num_multimorb","gender","smokstatus",
             "bmi","age",
             "ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
             "d2013","d2016")]  
extra<-extra %>%
  group_by(modalcls) %>%
  summarise("% patients"=mean(clsprop),
            "% greater deprivation"=sum(imd2015_5>3,na.rm = T)/n()*100,
            "Avg. age"=mean(age),
            "Med. age"=median(age),
            "Q1 age"=quantile(age,0.25),
            "Q3 age"=quantile(age,0.75),
            "% female"=sum(gender==2)/n()*100,
            "% current smokers"=sum(smokstatus==1,na.rm = T)/n()*100,
            "Avg. BMI"=mean(bmi,na.rm = T),
            "Med. BMI"=median(bmi,na.rm = T),
            "Q1 BMI"=quantile(bmi,0.25,na.rm = T),
            "Q3 BMI"=quantile(bmi,0.75,na.rm = T),
            "Avg. # morbidities"=mean(num_multimorb),
            "Med. # morbidities"=median(num_multimorb),
            "Q1 # morbidities"=quantile(num_multimorb,0.25),
            "Q3 # morbidities"=quantile(num_multimorb,0.75),
            "Avg. # consultations"=mean(ctconsult12mons_filterdur),
            "Med. # consultations"=median(ctconsult12mons_filterdur),
            "Q1 # consultations"=quantile(ctconsult12mons_filterdur,0.25),
            "Q3 # consultations"=quantile(ctconsult12mons_filterdur,0.75),
            "Avg. # hospital spells"=mean(cthospspl_12mons),
            "Med. # hospital spells"=median(cthospspl_12mons),
            "Q1 # hospital spells"=quantile(cthospspl_12mons,0.25),
            "Q3 # hospital spells"=quantile(cthospspl_12mons,0.75),
            "Avg. # polypharmacy"=mean(ctprod12mons4plus_BNF),
            "Med. # polypharmacy"=median(ctprod12mons4plus_BNF),
            "Q1 # polypharmacy"=quantile(ctprod12mons4plus_BNF,0.25),
            "Q3 # polypharmacy"=quantile(ctprod12mons4plus_BNF,0.75),
            "% deaths in 2 years"=sum(d2013==1)/n()*100,
            "% deaths in 5 years"=sum(d2016==1)/n()*100) 

extra<-as.data.frame(extra[order(-extra$"% patients"),]) # order by highest to lowest on class size
extra<-subset(extra,select=-c(modalcls))


# prepare row and col names for annotation #
rownames(extra) <-c("DEP,ANX,PNC","PNC,HL,HYP","AST,IBS,DEP",
                    "IBS,DEP,HL","PSM,ALC,DEP")
colnames(a)<-c("DEP,ANX,PNC","PNC,HL,HYP","AST,IBS,DEP",
               "IBS,DEP,HL","PSM,ALC,DEP")

# Descriptive stats for class features #

write.table(extra, "C:\\MM_proc_data\\Results\\unsup_clusterfeaturetrain1844.txt",sep="\t")

# output cluster profiles #
a1<- tibble::rownames_to_column(a, "Morbidityfull")
a2<-merge(a1,prev[,c("Morbidityfull","Prop.")],by="Morbidityfull")
a2<-as.data.frame(a2[order(-a2$"Prop."),]) # order by highest to lowest on prop.(prevalence)

profile_train_1844<-a2

write.table(a2, "C:\\MM_proc_data\\Results\\unsup_clusterprofiletrain1844.txt",sep="\t")

# prepare colours for plot
myBlues<- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100)

anno_colors <- list("% deaths in 5 years" =rev(heat.colors(100)), 
                    "% deaths in 2 years"=rev(heat.colors(100)),
                    "Avg. # polypharmacy"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # hospital spells"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # consultations"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "% current smokers"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% greater deprivation"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% patients"=colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))(100))

emf(file="C:\\MM_proc_data\\Results\\map_unsup_patid1844.emf",emfPlus = T)
mapunsup_patid<-pheatmap::pheatmap(a,cluster_rows = F, cluster_cols = F, cellwidth = 15, 
                               fontsize = 6,color=myBlues,annotation_col = extra[,c(1,2,5,8:12)],
                               annotation_colors = anno_colors,display_numbers = T,legend=F)
mapunsup_patid

dev.off()

# 45-64 #

a<-lcmodel_mm_train_unsup4564
a<-merge(a,classsize_unsup4564[,c("top3morb_3","clslabel")],by="clslabel")

a$Var2<-paste(round(100*a$clsprop,2),"%",";",a$top3morb_3)

# convert data from long to wide
a<-a[order(-a$clsprop,-a$value),] # order by highest to lowest on class size and conditional prob.
a<- reshape(a[,c("Morbidityfull","value","Var2")], 
            timevar = "Var2",
            idvar = c("Morbidityfull"),
            direction = "wide")
a2 <- a[,-1]
rownames(a2) <- a[,1]
a<-a2
names(a) <- substring(names(a), 7) # now a's colnames is very long
heatdata_unsup<-a # store map info into heatdata; plotted graph is based on "a"
colnames(a) <- paste('C',1:5,sep='') # this is order from highest to lowest on class size

# find extra characteristics for each cluster #
c<-lcapatid_train_unsup4564

extra <-c[,c("modalcls","clsprop","imd2015_5","num_multimorb","gender","smokstatus",
             "bmi","age",
             "ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
             "d2013","d2016")]  
extra<-extra %>%
  group_by(modalcls) %>%
  summarise("% patients"=mean(clsprop),
            "% greater deprivation"=sum(imd2015_5>3,na.rm = T)/n()*100,
            "Avg. age"=mean(age),
            "Med. age"=median(age),
            "Q1 age"=quantile(age,0.25),
            "Q3 age"=quantile(age,0.75),
            "% female"=sum(gender==2)/n()*100,
            "% current smokers"=sum(smokstatus==1,na.rm = T)/n()*100,
            "Avg. BMI"=mean(bmi,na.rm = T),
            "Med. BMI"=median(bmi,na.rm = T),
            "Q1 BMI"=quantile(bmi,0.25,na.rm = T),
            "Q3 BMI"=quantile(bmi,0.75,na.rm = T),
            "Avg. # morbidities"=mean(num_multimorb),
            "Med. # morbidities"=median(num_multimorb),
            "Q1 # morbidities"=quantile(num_multimorb,0.25),
            "Q3 # morbidities"=quantile(num_multimorb,0.75),
            "Avg. # consultations"=mean(ctconsult12mons_filterdur),
            "Med. # consultations"=median(ctconsult12mons_filterdur),
            "Q1 # consultations"=quantile(ctconsult12mons_filterdur,0.25),
            "Q3 # consultations"=quantile(ctconsult12mons_filterdur,0.75),
            "Avg. # hospital spells"=mean(cthospspl_12mons),
            "Med. # hospital spells"=median(cthospspl_12mons),
            "Q1 # hospital spells"=quantile(cthospspl_12mons,0.25),
            "Q3 # hospital spells"=quantile(cthospspl_12mons,0.75),
            "Avg. # polypharmacy"=mean(ctprod12mons4plus_BNF),
            "Med. # polypharmacy"=median(ctprod12mons4plus_BNF),
            "Q1 # polypharmacy"=quantile(ctprod12mons4plus_BNF,0.25),
            "Q3 # polypharmacy"=quantile(ctprod12mons4plus_BNF,0.75),
            "% deaths in 2 years"=sum(d2013==1)/n()*100,
            "% deaths in 5 years"=sum(d2016==1)/n()*100) 

extra<-as.data.frame(extra[order(-extra$"% patients"),]) # order by highest to lowest on class size
extra<-subset(extra,select=-c(modalcls))


# prepare row and col names for annotation #
rownames(extra) <-c("HYP,DIA,PNC","IBS,HL,PNC","DEP,PNC,ANX",
                    "AST,COPD,PNC","ALC,PSM,PNC")
colnames(a)<-c("HYP,DIA,PNC","IBS,HL,PNC","DEP,PNC,ANX",
               "AST,COPD,PNC","ALC,PSM,PNC")

# Descriptive stats for class features #

write.table(extra, "C:\\MM_proc_data\\Results\\unsup_clusterfeaturetrain4564.txt",sep="\t")

# output cluster profiles #
a1<- tibble::rownames_to_column(a, "Morbidityfull")
a2<-merge(a1,prev[,c("Morbidityfull","Prop.")],by="Morbidityfull")
a2<-as.data.frame(a2[order(-a2$"Prop."),]) # order by highest to lowest on prop.(prevalence)

profile_train_4564<-a2

write.table(a2, "C:\\MM_proc_data\\Results\\unsup_clusterprofiletrain4564.txt",sep="\t")

# prepare colours for plot
myBlues<- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100)

anno_colors <- list("% deaths in 5 years" =rev(heat.colors(100)), 
                    "% deaths in 2 years"=rev(heat.colors(100)),
                    "Avg. # polypharmacy"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # hospital spells"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # consultations"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "% current smokers"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% greater deprivation"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% patients"=colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))(100))


emf(file="C:\\MM_proc_data\\Results\\map_unsup_patid4564.emf",emfPlus = F)
mapunsup_patid<-pheatmap::pheatmap(a,cluster_rows = F, cluster_cols = F, cellwidth = 15, 
                                   fontsize = 6,color=myBlues,annotation_col = extra[,c(1,2,5,8:12)],
                                   annotation_colors = anno_colors,display_numbers = T)
mapunsup_patid

dev.off()

# 65-84 #

a<-lcmodel_mm_train_unsup6584
a<-merge(a,classsize_unsup6584[,c("top3morb_3","clslabel")],by="clslabel")

a$Var2<-paste(round(100*a$clsprop,2),"%",";",a$top3morb_3)

# convert data from long to wide
a<-a[order(-a$clsprop,-a$value),] # order by highest to lowest on class size and conditional prob.
a<- reshape(a[,c("Morbidityfull","value","Var2")], 
            timevar = "Var2",
            idvar = c("Morbidityfull"),
            direction = "wide")
a2 <- a[,-1]
rownames(a2) <- a[,1]
a<-a2
names(a) <- substring(names(a), 7) # now a's colnames is very long
heatdata_unsup<-a # store map info into heatdata; plotted graph is based on "a"
colnames(a) <- paste('C',1:6,sep='') # this is order from highest to lowest on class size

# find extra characteristics for each cluster #
c<-lcapatid_train_unsup6584

extra <-c[,c("modalcls","clsprop","imd2015_5","num_multimorb","gender","smokstatus",
             "bmi","age",
             "ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
             "d2013","d2016")]  
extra<-extra %>%
  group_by(modalcls) %>%
  summarise("% patients"=mean(clsprop),
            "% greater deprivation"=sum(imd2015_5>3,na.rm = T)/n()*100,
            "Avg. age"=mean(age),
            "Med. age"=median(age),
            "Q1 age"=quantile(age,0.25),
            "Q3 age"=quantile(age,0.75),
            "% female"=sum(gender==2)/n()*100,
            "% current smokers"=sum(smokstatus==1,na.rm = T)/n()*100,
            "Avg. BMI"=mean(bmi,na.rm = T),
            "Med. BMI"=median(bmi,na.rm = T),
            "Q1 BMI"=quantile(bmi,0.25,na.rm = T),
            "Q3 BMI"=quantile(bmi,0.75,na.rm = T),
            "Avg. # morbidities"=mean(num_multimorb),
            "Med. # morbidities"=median(num_multimorb),
            "Q1 # morbidities"=quantile(num_multimorb,0.25),
            "Q3 # morbidities"=quantile(num_multimorb,0.75),
            "Avg. # consultations"=mean(ctconsult12mons_filterdur),
            "Med. # consultations"=median(ctconsult12mons_filterdur),
            "Q1 # consultations"=quantile(ctconsult12mons_filterdur,0.25),
            "Q3 # consultations"=quantile(ctconsult12mons_filterdur,0.75),
            "Avg. # hospital spells"=mean(cthospspl_12mons),
            "Med. # hospital spells"=median(cthospspl_12mons),
            "Q1 # hospital spells"=quantile(cthospspl_12mons,0.25),
            "Q3 # hospital spells"=quantile(cthospspl_12mons,0.75),
            "Avg. # polypharmacy"=mean(ctprod12mons4plus_BNF),
            "Med. # polypharmacy"=median(ctprod12mons4plus_BNF),
            "Q1 # polypharmacy"=quantile(ctprod12mons4plus_BNF,0.25),
            "Q3 # polypharmacy"=quantile(ctprod12mons4plus_BNF,0.75),
            "% deaths in 2 years"=sum(d2013==1)/n()*100,
            "% deaths in 5 years"=sum(d2016==1)/n()*100) 


extra<-as.data.frame(extra[order(-extra$"% patients"),]) # order by highest to lowest on class size
extra<-subset(extra,select=-c(modalcls))


# prepare row and col names for annotation #
rownames(extra) <-c("HYP,DIA,CKD","HL,PSD,IBS","DEP,PNC,ANX",
                    "CHD,AF,DIA","COPD,AST,PNC","PNC,CHD,DEP")
colnames(a)<-c("HYP,DIA,CKD","HL,PSD,IBS","DEP,PNC,ANX",
               "CHD,AF,DIA","COPD,AST,PNC","PNC,CHD,DEP")

# Descriptive stats for class features #

write.table(extra, "C:\\MM_proc_data\\Results\\unsup_clusterfeaturetrain6584.txt",sep="\t")

# output cluster profiles #
a1<- tibble::rownames_to_column(a, "Morbidityfull")
a2<-merge(a1,prev[,c("Morbidityfull","Prop.")],by="Morbidityfull")
a2<-as.data.frame(a2[order(-a2$"Prop."),]) # order by highest to lowest on prop.(prevalence)

profile_train_6584<-a2

write.table(a2, "C:\\MM_proc_data\\Results\\unsup_clusterprofiletrain6584.txt",sep="\t")

# prepare colours for plot
myBlues<- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100)

anno_colors <- list("% deaths in 5 years" =rev(heat.colors(100)), 
                    "% deaths in 2 years"=rev(heat.colors(100)),
                    "Avg. # polypharmacy"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # hospital spells"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # consultations"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "% current smokers"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% greater deprivation"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% patients"=colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))(100))


emf(file="C:\\MM_proc_data\\Results\\map_unsup_patid6584.emf",emfPlus = F)
mapunsup_patid<-pheatmap::pheatmap(a,cluster_rows = F, cluster_cols = F, cellwidth = 15, 
                                   fontsize = 6,color=myBlues,annotation_col = extra[,c(1,2,5,8:12)],
                                   annotation_colors = anno_colors,display_numbers = T)
mapunsup_patid

dev.off()

# 85 plus #

a<-lcmodel_mm_train_unsup85plus
a<-merge(a,classsize_unsup85plus[,c("top3morb_3","clslabel")],by="clslabel")

a$Var2<-paste(round(100*a$clsprop,2),"%",";",a$top3morb_3)

# convert data from long to wide
a<-a[order(-a$clsprop,-a$value),] # order by highest to lowest on class size and conditional prob.
a<- reshape(a[,c("Morbidityfull","value","Var2")], 
            timevar = "Var2",
            idvar = c("Morbidityfull"),
            direction = "wide")
a2 <- a[,-1]
rownames(a2) <- a[,1]
a<-a2
names(a) <- substring(names(a), 7) # now a's colnames is very long
heatdata_unsup<-a # store map info into heatdata; plotted graph is based on "a"
colnames(a) <- paste('C',1:4,sep='') # this is order from highest to lowest on class size

# find extra characteristics for each cluster #
c<-lcapatid_train_unsup85plus

extra <-c[,c("modalcls","clsprop","imd2015_5","num_multimorb","gender","smokstatus",
             "bmi","age",
             "ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
             "d2013","d2016")]  
extra<-extra %>%
  group_by(modalcls) %>%
  summarise("% patients"=mean(clsprop),
            "% greater deprivation"=sum(imd2015_5>3,na.rm = T)/n()*100,
            "Avg. age"=mean(age),
            "Med. age"=median(age),
            "Q1 age"=quantile(age,0.25),
            "Q3 age"=quantile(age,0.75),
            "% female"=sum(gender==2)/n()*100,
            "% current smokers"=sum(smokstatus==1,na.rm = T)/n()*100,
            "Avg. BMI"=mean(bmi,na.rm = T),
            "Med. BMI"=median(bmi,na.rm = T),
            "Q1 BMI"=quantile(bmi,0.25,na.rm = T),
            "Q3 BMI"=quantile(bmi,0.75,na.rm = T),
            "Avg. # morbidities"=mean(num_multimorb),
            "Med. # morbidities"=median(num_multimorb),
            "Q1 # morbidities"=quantile(num_multimorb,0.25),
            "Q3 # morbidities"=quantile(num_multimorb,0.75),
            "Avg. # consultations"=mean(ctconsult12mons_filterdur),
            "Med. # consultations"=median(ctconsult12mons_filterdur),
            "Q1 # consultations"=quantile(ctconsult12mons_filterdur,0.25),
            "Q3 # consultations"=quantile(ctconsult12mons_filterdur,0.75),
            "Avg. # hospital spells"=mean(cthospspl_12mons),
            "Med. # hospital spells"=median(cthospspl_12mons),
            "Q1 # hospital spells"=quantile(cthospspl_12mons,0.25),
            "Q3 # hospital spells"=quantile(cthospspl_12mons,0.75),
            "Avg. # polypharmacy"=mean(ctprod12mons4plus_BNF),
            "Med. # polypharmacy"=median(ctprod12mons4plus_BNF),
            "Q1 # polypharmacy"=quantile(ctprod12mons4plus_BNF,0.25),
            "Q3 # polypharmacy"=quantile(ctprod12mons4plus_BNF,0.75),
            "% deaths in 2 years"=sum(d2013==1)/n()*100,
            "% deaths in 5 years"=sum(d2016==1)/n()*100)
extra<-as.data.frame(extra[order(-extra$"% patients"),]) # order by highest to lowest on class size
extra<-subset(extra,select=-c(modalcls))


# prepare row and col names for annotation #
rownames(extra) <-c("HYP,HL,DIA","PNC,DEP,CSP","HF,CHD,AF",
                    "AST,COPD,PNC")
colnames(a)<-c("HYP,HL,DIA","PNC,DEP,CSP","HF,CHD,AF",
               "AST,COPD,PNC")

# Descriptive stats for class features #

write.table(extra, "C:\\MM_proc_data\\Results\\unsup_clusterfeaturetrain85plus.txt",sep="\t")

# output cluster profiles #
a1<- tibble::rownames_to_column(a, "Morbidityfull")
a2<-merge(a1,prev[,c("Morbidityfull","Prop.")],by="Morbidityfull")
a2<-as.data.frame(a2[order(-a2$"Prop."),]) # order by highest to lowest on prop.(prevalence)

profile_train_85plus<-a2

write.table(a2, "C:\\MM_proc_data\\Results\\unsup_clusterprofiletrain85plus.txt",sep="\t")

# prepare colours for plot
myBlues<- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100)

anno_colors <- list("% deaths in 5 years" =rev(heat.colors(100)), 
                    "% deaths in 2 years"=rev(heat.colors(100)),
                    "Avg. # polypharmacy"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # hospital spells"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # consultations"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "% current smokers"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% greater deprivation"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% patients"=colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))(100))


emf(file="C:\\MM_proc_data\\Results\\map_unsup_patid85plus.emf",emfPlus = F)
mapunsup_patid<-pheatmap::pheatmap(a,cluster_rows = F, cluster_cols = F, cellwidth = 15, 
                                   fontsize = 6,color=myBlues,annotation_col = extra[,c(1,2,5,8:12)],
                                   annotation_colors = anno_colors,display_numbers = T)
mapunsup_patid

dev.off()

save.image(file='C:\\MM_proc_data\\final_patid_strata_final_filterdur.RData')                                                                                                                                                                    

#####################################################################################################
# Clusterd heatmap #

# Unsupervised : test# 

# 18-44 #

a<-lcmodel_mm_test_unsup1844
a<-merge(a,classsize_test_unsup1844[,c("top3morb_3","clslabel")],by="clslabel")

a$Var2<-paste(round(100*a$clsprop,2),"%",";",a$top3morb_3)

# convert data from long to wide
a<-a[order(-a$clsprop,-a$value),] # order by highest to lowest on class size and conditional prob.
a<- reshape(a[,c("Morbidityfull","value","Var2")], 
            timevar = "Var2",
            idvar = c("Morbidityfull"),
            direction = "wide")
a2 <- a[,-1]
rownames(a2) <- a[,1]
a<-a2
names(a) <- substring(names(a), 7) # now a's colnames is very long
heatdata_unsup<-a # store map info into heatdata; plotted graph is based on "a"
colnames(a) <- paste('C',1:5,sep='') # this is order from highest to lowest on class size

# find extra characteristics for each cluster #
c<-lcapatid_test_unsup1844

## smokestatus final order: 1= current, 2= never, 3 = ex

extra <-c[,c("modalcls","clsprop","imd2015_5","num_multimorb","gender","smokstatus",
             "bmi","age",
             "ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
             "d2013","d2016")]  
extra<-extra %>%
  group_by(modalcls) %>%
  summarise("% patients"=mean(clsprop),
            "% greater deprivation"=sum(imd2015_5>3,na.rm = T)/n()*100,
            "Med. age"=median(age),
            "% female"=sum(gender==2)/n()*100,
            "% current smokers"=sum(smokstatus==1,na.rm = T)/n()*100,
            "Med. BMI"=median(bmi,na.rm = T),
            "Med. # morbidities"=median(num_multimorb),
            "Avg. # consultations"=mean(ctconsult12mons_filterdur),
            "Avg. # hospital spells"=mean(cthospspl_12mons),
            "Avg. # polypharmacy"=mean(ctprod12mons4plus_BNF),
            "% deaths in 2 years"=sum(d2013==1)/n()*100,
            "% deaths in 5 years"=sum(d2016==1)/n()*100) 

extra<-as.data.frame(extra[order(-extra$"% patients"),]) # order by highest to lowest on class size
extra<-subset(extra,select=-c(modalcls))


# prepare row and col names for annotation #
rownames(extra) <-c("DEP,ANX,PNC","AST,HL,PNC","IBS,DEP,AST",
                    "HYP,DIA,DEP","PSM,ALC,PNC")
colnames(a)<-c("DEP,ANX,PNC","AST,HL,PNC","IBS,DEP,AST",
               "HYP,DIA,DEP","PSM,ALC,PNC")

# Descriptive stats for class features #

write.table(extra, "C:\\MM_proc_data\\Results\\unsup_clusterfeaturetest1844.txt",sep="\t")

# output cluster profiles #
a1<- tibble::rownames_to_column(a, "Morbidityfull")
a2<-merge(a1,prev[,c("Morbidityfull","Prop.")],by="Morbidityfull")
a2<-as.data.frame(a2[order(-a2$"Prop."),]) # order by highest to lowest on prop.(prevalence)

profile_test_1844<-a2

write.table(a2, "C:\\MM_proc_data\\Results\\unsup_clusterprofiletest1844.txt",sep="\t")

# prepare colours for plot
myBlues<- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100)

anno_colors <- list("% deaths in 5 years" =rev(heat.colors(100)), 
                    "% deaths in 2 years"=rev(heat.colors(100)),
                    "Avg. # polypharmacy"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # hospital spells"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # consultations"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "% current smokers"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% greater deprivation"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% patients"=colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))(100))

emf(file="C:\\MM_proc_data\\Results\\map_unsup_test_patid1844.emf",emfPlus = T)
mapunsup_patid<-pheatmap::pheatmap(a,cluster_rows = F, cluster_cols = F, cellwidth = 15, 
                                   fontsize = 6,color=myBlues,annotation_col = extra[,c(1,2,5,8:12)],
                                   annotation_colors = anno_colors,display_numbers = T)
mapunsup_patid

dev.off()

# 45-64 #

a<-lcmodel_mm_test_unsup4564
a<-merge(a,classsize_test_unsup4564[,c("top3morb_3","clslabel")],by="clslabel")

a$Var2<-paste(round(100*a$clsprop,2),"%",";",a$top3morb_3)

# convert data from long to wide
a<-a[order(-a$clsprop,-a$value),] # order by highest to lowest on class size and conditional prob.
a<- reshape(a[,c("Morbidityfull","value","Var2")], 
            timevar = "Var2",
            idvar = c("Morbidityfull"),
            direction = "wide")
a2 <- a[,-1]
rownames(a2) <- a[,1]
a<-a2
names(a) <- substring(names(a), 7) # now a's colnames is very long
heatdata_unsup<-a # store map info into heatdata; plotted graph is based on "a"
colnames(a) <- paste('C',1:5,sep='') # this is order from highest to lowest on class size

# find extra characteristics for each cluster #
c<-lcapatid_test_unsup4564

extra <-c[,c("modalcls","clsprop","imd2015_5","num_multimorb","gender","smokstatus",
             "bmi","age",
             "ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
             "d2013","d2016")]  
extra<-extra %>%
  group_by(modalcls) %>%
  summarise("% patients"=mean(clsprop),
            "% greater deprivation"=sum(imd2015_5>3,na.rm = T)/n()*100,
            "Med. age"=median(age),
            "% female"=sum(gender==2)/n()*100,
            "% current smokers"=sum(smokstatus==1,na.rm = T)/n()*100,
            "Med. BMI"=median(bmi,na.rm = T),
            "Med. # morbidities"=median(num_multimorb),
            "Avg. # consultations"=mean(ctconsult12mons_filterdur),
            "Avg. # hospital spells"=mean(cthospspl_12mons),
            "Avg. # polypharmacy"=mean(ctprod12mons4plus_BNF),
            "% deaths in 2 years"=sum(d2013==1)/n()*100,
            "% deaths in 5 years"=sum(d2016==1)/n()*100) 

extra<-as.data.frame(extra[order(-extra$"% patients"),]) # order by highest to lowest on class size
extra<-subset(extra,select=-c(modalcls))


# prepare row and col names for annotation #
rownames(extra) <-c("HYP,DIA,CHD","IBS,HL,AST","DEP,PNC,ANX",
                    "PNC,CHD,DIA","ALC,PSM,PNC")
colnames(a)<-c("HYP,DIA,CHD","IBS,HL,AST","DEP,PNC,ANX",
               "PNC,CHD,DIA","ALC,PSM,PNC")

# Descriptive stats for class features #

write.table(extra, "C:\\MM_proc_data\\Results\\unsup_clusterfeaturetest4564.txt",sep="\t")

# output cluster profiles #
a1<- tibble::rownames_to_column(a, "Morbidityfull")
a2<-merge(a1,prev[,c("Morbidityfull","Prop.")],by="Morbidityfull")
a2<-as.data.frame(a2[order(-a2$"Prop."),]) # order by highest to lowest on prop.(prevalence)

profile_test_4564<-a2

write.table(a2, "C:\\MM_proc_data\\Results\\unsup_clusterprofiletest4564.txt",sep="\t")

# prepare colours for plot
myBlues<- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100)

anno_colors <- list("% deaths in 5 years" =rev(heat.colors(100)), 
                    "% deaths in 2 years"=rev(heat.colors(100)),
                    "Avg. # polypharmacy"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # hospital spells"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # consultations"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "% current smokers"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% greater deprivation"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% patients"=colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))(100))


emf(file="C:\\MM_proc_data\\Results\\map_unsup_test_patid4564.emf",emfPlus = F)
mapunsup_patid<-pheatmap::pheatmap(a,cluster_rows = F, cluster_cols = F, cellwidth = 15, 
                                   fontsize = 6,color=myBlues,annotation_col = extra[,c(1,2,5,8:12)],
                                   annotation_colors = anno_colors,display_numbers = T)
mapunsup_patid

dev.off()

# 65-84 #

a<-lcmodel_mm_test_unsup6584
a<-merge(a,classsize_test_unsup6584[,c("top3morb_3","clslabel")],by="clslabel")

a$Var2<-paste(round(100*a$clsprop,2),"%",";",a$top3morb_3)

# convert data from long to wide
a<-a[order(-a$clsprop,-a$value),] # order by highest to lowest on class size and conditional prob.
a<- reshape(a[,c("Morbidityfull","value","Var2")], 
            timevar = "Var2",
            idvar = c("Morbidityfull"),
            direction = "wide")
a2 <- a[,-1]
rownames(a2) <- a[,1]
a<-a2
names(a) <- substring(names(a), 7) # now a's colnames is very long
heatdata_unsup<-a # store map info into heatdata; plotted graph is based on "a"
colnames(a) <- paste('C',1:6,sep='') # this is order from highest to lowest on class size

# find extra characteristics for each cluster #
c<-lcapatid_test_unsup6584

extra <-c[,c("modalcls","clsprop","imd2015_5","num_multimorb","gender","smokstatus",
             "bmi","age",
             "ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
             "d2013","d2016")]  
extra<-extra %>%
  group_by(modalcls) %>%
  summarise("% patients"=mean(clsprop),
            "% greater deprivation"=sum(imd2015_5>3,na.rm = T)/n()*100,
            "Med. age"=median(age),
            "% female"=sum(gender==2)/n()*100,
            "% current smokers"=sum(smokstatus==1,na.rm = T)/n()*100,
            "Med. BMI"=median(bmi,na.rm = T),
            "Med. # morbidities"=median(num_multimorb),
            "Avg. # consultations"=mean(ctconsult12mons_filterdur),
            "Avg. # hospital spells"=mean(cthospspl_12mons),
            "Avg. # polypharmacy"=mean(ctprod12mons4plus_BNF),
            "% deaths in 2 years"=sum(d2013==1)/n()*100,
            "% deaths in 5 years"=sum(d2016==1)/n()*100) 


extra<-as.data.frame(extra[order(-extra$"% patients"),]) # order by highest to lowest on class size
extra<-subset(extra,select=-c(modalcls))


# prepare row and col names for annotation #
rownames(extra) <-c("HYP,DIA,CKD","HL,IBS,DIV","CHD,AF,DIA",
                    "DEP,PNC,ANX","AST,COPD,PNC","PNC,DEP,CHD")
colnames(a)<-c("HYP,DIA,CKD","HL,IBS,DIV","CHD,AF,DIA",
               "DEP,PNC,ANX","AST,COPD,PNC","PNC,DEP,CHD")

# Descriptive stats for class features #

write.table(extra, "C:\\MM_proc_data\\Results\\unsup_clusterfeaturetest6584.txt",sep="\t")

# output cluster profiles #
a1<- tibble::rownames_to_column(a, "Morbidityfull")
a2<-merge(a1,prev[,c("Morbidityfull","Prop.")],by="Morbidityfull")
a2<-as.data.frame(a2[order(-a2$"Prop."),]) # order by highest to lowest on prop.(prevalence)

profile_test_6584<-a2

write.table(a2, "C:\\MM_proc_data\\Results\\unsup_clusterprofiletest6584.txt",sep="\t")

# prepare colours for plot
myBlues<- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100)

anno_colors <- list("% deaths in 5 years" =rev(heat.colors(100)), 
                    "% deaths in 2 years"=rev(heat.colors(100)),
                    "Avg. # polypharmacy"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # hospital spells"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # consultations"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "% current smokers"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% greater deprivation"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% patients"=colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))(100))


emf(file="C:\\MM_proc_data\\Results\\map_unsup_test_patid6584.emf",emfPlus = F)
mapunsup_patid<-pheatmap::pheatmap(a,cluster_rows = F, cluster_cols = F, cellwidth = 15, 
                                   fontsize = 6,color=myBlues,annotation_col = extra[,c(1,2,5,8:12)],
                                   annotation_colors = anno_colors,display_numbers = T)
mapunsup_patid

dev.off()

# 85 plus #

a<-lcmodel_mm_test_unsup85plus
a<-merge(a,classsize_test_unsup85plus[,c("top3morb_3","clslabel")],by="clslabel")

a$Var2<-paste(round(100*a$clsprop,2),"%",";",a$top3morb_3)

# convert data from long to wide
a<-a[order(-a$clsprop,-a$value),] # order by highest to lowest on class size and conditional prob.
a<- reshape(a[,c("Morbidityfull","value","Var2")], 
            timevar = "Var2",
            idvar = c("Morbidityfull"),
            direction = "wide")
a2 <- a[,-1]
rownames(a2) <- a[,1]
a<-a2
names(a) <- substring(names(a), 7) # now a's colnames is very long
heatdata_unsup<-a # store map info into heatdata; plotted graph is based on "a"
colnames(a) <- paste('C',1:4,sep='') # this is order from highest to lowest on class size

# find extra characteristics for each cluster #
c<-lcapatid_test_unsup85plus

extra <-c[,c("modalcls","clsprop","imd2015_5","num_multimorb","gender","smokstatus",
             "bmi","age",
             "ctconsult12mons_filterdur","cthospspl_12mons","ctprod12mons4plus_BNF",
             "d2013","d2016")]  
extra<-extra %>%
  group_by(modalcls) %>%
  summarise("% patients"=mean(clsprop),
            "% greater deprivation"=sum(imd2015_5>3,na.rm = T)/n()*100,
            "Med. age"=median(age),
            "% female"=sum(gender==2)/n()*100,
            "% current smokers"=sum(smokstatus==1,na.rm = T)/n()*100,
            "Med. BMI"=median(bmi,na.rm = T),
            "Med. # morbidities"=median(num_multimorb),
            "Avg. # consultations"=mean(ctconsult12mons_filterdur),
            "Avg. # hospital spells"=mean(cthospspl_12mons),
            "Avg. # polypharmacy"=mean(ctprod12mons4plus_BNF),
            "% deaths in 2 years"=sum(d2013==1)/n()*100,
            "% deaths in 5 years"=sum(d2016==1)/n()*100) 
extra<-as.data.frame(extra[order(-extra$"% patients"),]) # order by highest to lowest on class size
extra<-subset(extra,select=-c(modalcls))


# prepare row and col names for annotation #
rownames(extra) <-c("HYP,CHD,HL","DEP,PNC,CSP","PNC,DEP,CSP",
                    "COPD,AST,PSD")
colnames(a)<-c("HYP,CHD,HL","DEP,PNC,CSP","PNC,DEP,CSP",
               "COPD,AST,PSD")

# Descriptive stats for class features #

write.table(extra, "C:\\MM_proc_data\\Results\\unsup_clusterfeaturetest85plus.txt",sep="\t")

# output cluster profiles #
a1<- tibble::rownames_to_column(a, "Morbidityfull")
a2<-merge(a1,prev[,c("Morbidityfull","Prop.")],by="Morbidityfull")
a2<-as.data.frame(a2[order(-a2$"Prop."),]) # order by highest to lowest on prop.(prevalence)

profile_test_85plus<-a2

write.table(a2, "C:\\MM_proc_data\\Results\\unsup_clusterprofiletest85plus.txt",sep="\t")

# prepare colours for plot
myBlues<- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100)

anno_colors <- list("% deaths in 5 years" =rev(heat.colors(100)), 
                    "% deaths in 2 years"=rev(heat.colors(100)),
                    "Avg. # polypharmacy"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # hospital spells"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "Avg. # consultations"=colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(100),
                    "% current smokers"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% greater deprivation"=colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"))(100),
                    "% patients"=colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))(100))


emf(file="C:\\MM_proc_data\\Results\\map_unsup_test_patid85plus.emf",emfPlus = F)
mapunsup_patid<-pheatmap::pheatmap(a,cluster_rows = F, cluster_cols = F, cellwidth = 15, 
                                   fontsize = 6,color=myBlues,annotation_col = extra[,c(1,2,5,8:12)],
                                   annotation_colors = anno_colors,display_numbers = T)
mapunsup_patid

dev.off()

#######################################################################################

# comparison of training and test cluster profile solutions # 
# sqrt(Jensen-shannon divergence statistic) & pearson's correlation coefficient
# for pairwise comparison #

# note, for test_85plus, there is no migrane so need to insert this in the disease list#
top<-profile_test_85plus[1:33,]
b<- profile_test_85plus[1,]
b[1,1]<-"Migraine"
b[1,2:6]<-c(0,0,0,0,0.3911466)
bottom<-profile_test_85plus[34:37,]
c<-rbind(top,b,bottom)
c<-c[order(-c$Prop),]
profile_test_85plus<-c #corrected, with 38 LTCs


# 18-44 #
match_1844<-JSD_pair(profile_train_1844[2:6],profile_test_1844[2:6])
write.table(match_1844, "C:\\MM_proc_data\\Results\\match_1844.txt",sep="\t")

a<-cor(profile_train_1844[2:6],profile_test_1844[2:6])
write.table(a, "C:\\MM_proc_data\\Results\\match_1844.txt",sep="\t",append=T)


# 45-64 #
match_4564<-JSD_pair(profile_train_4564[2:6],profile_test_4564[2:6])
write.table(match_4564, "C:\\MM_proc_data\\Results\\match_4564.txt",sep="\t")

a<-cor(profile_train_4564[2:6],profile_test_4564[2:6])
write.table(a, "C:\\MM_proc_data\\Results\\match_4564.txt",sep="\t",append=T)

# 65-84 #
match_6584<-JSD_pair(profile_train_6584[2:7],profile_test_6584[2:7])
write.table(match_6584, "C:\\MM_proc_data\\Results\\match_6584.txt",sep="\t")

a<-cor(profile_train_6584[2:7],profile_test_6584[2:7])
write.table(a, "C:\\MM_proc_data\\Results\\match_6584.txt",sep="\t",append=T)

# 85+ #
match_85plus<-JSD_pair(profile_train_85plus[2:5],profile_test_85plus[2:5])
write.table(match_85plus, "C:\\MM_proc_data\\Results\\match_85plus.txt",sep="\t")

a<-cor(profile_train_85plus[2:5],profile_test_85plus[2:5])
write.table(a, "C:\\MM_proc_data\\Results\\match_85plus.txt",sep="\t",append=T)

save.image(file='C:\\MM_proc_data\\final_patid_strata_final_filterdur.RData')

#######################################################################################
## END ## 


