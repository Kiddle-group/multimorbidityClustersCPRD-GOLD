# Last updated: 29 May 2019 #
# This file contains functions that are used across R files that accompany the paper 
# "Identification of clusters of multimorbid patients in UK general practice: a latent class analysis."s


dateNumeric <- function(dates,...){ 
  # function to turn dates into numbers, not sure if I still use
  
  as.numeric(as.Date(dates,...))
  
}

find_medcodes <- function(medcodes,lm_time=0){
  # given a list of medcodes, finds them in all tables, ATM focused on COPD diagnosis,
  # later will modify to work with landmark times (lm_time)
  
  # empty df to save results in
  relevant <- data.frame(patid=NA,eventdate=as.Date('1990-01-01'),medcode=NA)
  
  # if no landmark time specified, extract before COPD diagnosis
  if(lm_time==0){
    
    print('from COPD diagnosis')
    
    # for each table, find instances, store them in relevant (ffdf) 
    # and deal with errors if nothing is found
    
    print('Clinical')
    
    a <- tryCatch(
      {
        res <- clinical_before[clinical_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
        relevant <- ffdfappend(relevant,res)
      }, warning = function(w) {}, error = function(e) {}, finally = {}
    )
    
    print('Referral')
    
    a <- tryCatch(
      {
        res <- referral_before[referral_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
        relevant <- ffdfappend(relevant,res)
      }, warning = function(w) {}, error = function(e) {}, finally = {}
    )
    
    print('Test')
    
    a <- tryCatch(
      {
        res <- test_before[test_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
        relevant <- ffdfappend(relevant,res)
      }, warning = function(w) {}, error = function(e) {}, finally = {}
    )
    
    print('Immunisation')
    
    a <- tryCatch(
      {
        res <- immunisation_before[immunisation_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
        relevant <- ffdfappend(relevant,res)
      }, warning = function(w) {}, error = function(e) {}, finally = {}
    )
    
    
  } else {
    
    if(lm_time==2012){
      
      print('2012')
      
      # for each table, find instances, store them in relevant (ffdf) 
      # and deal with errors if nothing is found
      
      print('Clinical')
      
      a <- tryCatch(
        {
          res <- clinical2012_before[clinical2012_before$medcode[1:(dim(clinical2012_before))[1]] %in% medcodes$medcode,c('patid','eventdate','medcode')]
          relevant <- rbind(relevant,res)
        }, warning = function(w) {}, error = function(e) {}, finally = {}
      )
      
      print('Referral')
      
      a <- tryCatch(
        {
          res <- referral2012_before[referral2012_before$medcode[1:(dim(referral2012_before))[1]] %in% medcodes$medcode,c('patid','eventdate','medcode')]
          relevant <- rbind(relevant,res)
        }, warning = function(w) {}, error = function(e) {}, finally = {}
      )
      
      print('Test')
      
      a <- tryCatch(
        {
          res <- test2012_before[
            test2012_before$medcode[1:(dim(test2012_before))[1]] %in% medcodes$medcode,c('patid','eventdate','medcode')]
          relevant <- rbind(relevant,res)
        }, warning = function(w) {}, error = function(e) {}, finally = {}
      )
      
      print('Immunisation')
      
      a <- tryCatch(
        {
          res <- immunisation2012_before[immunisation2012_before$medcode[1:(dim(immunisation2012_before))[1]] %in% medcodes$medcode,c('patid','eventdate','medcode')]
          relevant <- rbind(relevant,res)
        }, warning = function(w) {}, error = function(e) {}, finally = {}
      )
      
      
    }

    
  }
  
  return(relevant)
  
}

silviaScoreM <- function(data){
# function to calculate Cambridge Multimorbidity Score for mortality
  
  score <- numeric(dim(data)[1])
  
  score <- score + (data$hyp == 'present') * -2.09
  score <- score + ((data$anx == 'present') | (data$dep == 'present')) * 7.04
  score <- score + (data$pnc == 'present') * 16.46
  score <- score + (data$hel == 'present') * -3.94
  score <- score + (data$ibs == 'present') * -1.33
  score <- score + (data$ast == 'present') * -2.73
  score <- score + (data$dia == 'present') * 10.23
  score <- score + (data$chd == 'present') * 4.22
  score <- score + (data$ckd == 'present') * 16.61
  score <- score + (data$atr == 'present') * 22.14
  score <- score + (data$con == 'present') * 35.42
  score <- score + (data$str == 'present') * 20.63
  score <- score + (data$copd == 'present') * 42.50
  score <- score + (data$rhe == 'present') * -0.39
  score <- score + (data$can == 'present') * 62
  score <- score + (data$ap == 'present') * 12.72
  score <- score + (data$hf == 'present') * 43.47
  score <- score + (data$dem == 'present') * 124.42
  score <- score + (data$scz == 'present') * 7.2
  score <- score + (data$epi == 'present') * 18.26
  
}

silviaScoreG <- function(data){
  # function to calculate general Cambridge Multimorbidity Score
  
  score <- numeric(dim(data)[1])
  
  score <- score + (data$hyp == 'present') * 0.08
  score <- score + ((data$anx == 'present') | (data$dep == 'present')) * 0.5
  score <- score + (data$pnc == 'present') * 0.92
  score <- score + (data$hel == 'present') * 0.09
  score <- score + (data$ibs == 'present') * 0.21
  score <- score + (data$ast == 'present') * 0.19
  score <- score + (data$dia == 'present') * 0.75
  score <- score + (data$chd == 'present') * 0.49
  score <- score + (data$ckd == 'present') * 0.53
  score <- score + (data$atr == 'present') * 1.34
  score <- score + (data$con == 'present') * 1.12
  score <- score + (data$str == 'present') * 0.80
  score <- score + (data$copd == 'present') * 1.46
  score <- score + (data$rhe == 'present') * 0.43
  score <- score + (data$can == 'present') * 1.53
  score <- score + (data$ap == 'present') * 0.65
  score <- score + (data$hf == 'present') * 1.18
  score <- score + (data$dem == 'present') * 2.50
  score <- score + (data$scz == 'present') * 0.64
  score <- score + (data$epi == 'present') * 0.92
  
}

# is a patient known to be dead by X years post COPD diagnosis?
dead_by <- function(years,censored,time){censored == 0 & years < time}

# how many patients in this subset are dead X years post COPD diagnosis?
num_dead <- function(subset,time){length(which(dead_by(surv_copd$years[subset],surv_copd$censored[subset],time)))}

# add up the number of true instances
sumBool <- function(a){sum(as.numeric(a))}

# count number of NA (missing) data
numNA <- function(a){length(which(is.na(a)))}


# given cox model and new data, predict probability of five years survival
pred_coxph <- function(fit,data){
  
  #fit <- coxph(SurvObj2 ~ age + gender + smoke,data=small)
  
  # multiply coefficients by patients variable values to get prognostic scores
  X_new <- model.matrix(fit$formula,data= data)
  X_b <- fit$coef%*%t(X_new[,-1])
  
  # extract baseline hazard and index for five year prediction
  bh <- basehaz(fit)
  five_years <- dim(bh)[1]
  
  # use cox model to extract probability
  exp(-bh[five_years,1])^exp(X_b)
  
  #plot(bh[,2],exp(-bh[,1]))
  #lines(bh[,2],exp(-bh[,1])^exp(X_b[1]),col='blue')
}

R2BGLiMS_format <- function(data,old=vector('list',1)){
  # I was going to use R2BGLiMS, but didn't in the end. But the format is useful for
  # other things. Essentially standardise all data.
  
  # saves me renaming some old code variables
  small <- data
  
  # prepare output
  out <- vector('list',4)
  names(out) <- c('data','scaling_factors','means','year5_scaled')
  
  # columns that should be converted to integers
  cols_int <- c(12:19,21:58)[-c(2,4,5,6,7,45)]
  
  # convert them to integers
  for (i in cols_int){
    
    # missingness indicators are T/F, not factor
    if(i != 14 & i != 19 & i != 58){
      
      small[,i] <- as.integer(small[,i])-1
      
    } else {
      
      small[,i] <- as.integer(small[,i])
      
    }
    
  }
  
  # deal with factors, creating contrasts
  formula <- as.formula(paste('SurvObj~',paste(colnames(small)[c(12:19,21:58)][-c(17,14,10,27,41,24,25)],collapse='+')))
  
  mat <- model.matrix(formula,data= small)[,-1]
  
  std_small <- cbind(small[,1:10],mat)
  
  
  # scale all variables by dividing by standard deviation, from this or a previous run
  for (i in 11:dim(std_small)[2]){
    
    if (length(old) == 4){
      
      out$scaling_factors[i-10] <- old$scaling_factors[i-10]
      
    } else {
      
      out$scaling_factors[i-10] <- sd(mat[,i-10])
      
    }
    
    std_small[,i] <- mat[,i-10]/out$scaling_factors[i-10]
    
  }
  
  
  
  #colnames(std_small)[17:47]
  
  
  # create pairwise interactions between co-morbidities
  for (i in 1:30){
    
    for (j in (i+1):31){
      
      if (i == 30 & j == 31){}
      else {
        
        str <- paste(colnames(std_small)[20:50][i],colnames(std_small)[20:50][j],sep='_')
        
        std_small[,str] <- std_small[,c(20:50)[i]] * std_small[,c(20:50)[j]]
        
      }
      
    }
    
  }
  
  
  # center variables
  for (i in c(11:514)){
    
    
    
    if (length(old) == 4){
      
      out$means[i-10] <- old$means[i-10]
      
    } else {
      
      out$means[i-10] <- mean(std_small[,i])
      
    }
    
    std_small[,i] <- std_small[,i] - out$means[i-10]
    
  }
  
  
  # only need for R2BGLiMs
  out$year5_scaled <- 5/max(std_small$years)
  std_small$years <- std_small$years/max(std_small$years)# Recommend scaling survival times to between 0 and 1
  
  std_small$dead <- 1-std_small$censored
  
  out$data <- std_small
  
  return(out)
  
}

# function to load BIG therapy data bits by bits into dataframe and to extract info. from each bit

bigTherapyExtract <- function(prods,chunkSize=10^7){
  
  # num rows in therapy
  N <- dim(therapy2012_before)[1]
  
  # we will only this many rows into R at a time
  n <- chunkSize
  
  # intialise ff vector, will use to say if a row matches prodcode list
  ind <- ff(vmode='boolean',length=N)
  
  # for each chunk
  for (i in 1:ceiling(N/n)){
    
    # define indexes of chunk, i.e. 1:n, (n+1):(2*n), etc
    ind2 <- ((i-1)*n + 1):(i*n) 
    
    # last chunk maybe too long
    if(ind2[n] > N){
      
      # if so, remove indexes higher than needed
      ind2 <- ind2[1:which(ind2==N)]
      
    }
    
    # using chunk indexes on ind and therapy2012_before to ensure
    # we don't load too large a vector into RAM
    ind[ind2] <- therapy2012_before$prodcode[ind2] %in% prods$prodcode
    
    # report on progress, including cumulative number of matches
    print(paste('chunk',i,'of',ceiling(N/n)))
    print(sum(ind))
    
  }
  
  # use ind, constructed by chunking, to extract relevant rows
  therapy_relevant <- therapy2012_before[ind,c('patid','eventdate','prodcode')]
  
  return(therapy_relevant)
  
}


# function to add co-morbidity data for codelists with ever recorded inclusion criteria
ever_recorded <- function(codes,...){
  
  # for each eligible subject, see if code is present
  relevant <- as.data.table.ffdf(find_medcodes(codes,...))
  
  # keep earliest event per patient
  setorderv(relevant,c('patid','eventdate'))
  relevant <- relevant[!duplicated(relevant[,'patid']),]
  
  relevant <- relevant[-1,]
  
  #extract first event per patient
  return(relevant)
  
}

find_latest <- function(relevant,ffdf=T){

  if (ffdf){relevant <- as.data.table.ffdf(relevant)}
  
  # keep latest event per patient
  setorder(relevant,patid,-eventdate)
  relevant <- relevant[!duplicated(relevant[,'patid']),]
 
  return(relevant) 
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# move some columns to last 
movetolast <- function(data, move) {
  data[c(setdiff(names(data), move), move)]
}

# LCA entropy (standardised)
entropy<-function (p) sum(-p*log(p))

# LCA bivariate residuals
bvr <- function(fit) {
  stopifnot(class(fit) == "poLCA")
  
  ov_names <- names(fit$predcell)[1:(ncol(fit$predcell) - 2)]
  ov_combn <- combn(ov_names, 2)
  
  get_bvr <- function(ov_pair) {
    form_obs <- as.formula(paste0("observed ~ ", ov_pair[1], " + ", ov_pair[2]))
    form_exp <- as.formula(paste0("expected ~ ", ov_pair[1], " + ", ov_pair[2]))
    
    counts_obs <- xtabs(form_obs, data = fit$predcell)
    counts_exp <- xtabs(form_exp, data = fit$predcell)
    
    bvr <- sum((counts_obs - counts_exp)^2 / counts_exp)
    
    bvr
  }
  
  bvr_pairs <- apply(ov_combn, 2, get_bvr)
  # names(bvr_pairs) <- apply(ov_combn, 2, paste, collapse = "<->")
  attr(bvr_pairs, "class") <- "dist"
  attr(bvr_pairs, "Size") <- length(ov_names)
  attr(bvr_pairs, "Labels") <- ov_names
  attr(bvr_pairs, "Diag") <- FALSE
  attr(bvr_pairs, "Upper") <- FALSE
  
  bvr_pairs
}

# To allocate heatmaps across a page
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

# # create a merge indicator for later use#
# 
# train_prac_mm$source<-"master"
# results_train_consult2$source<-"using"
# 
# a<- merge(x = train_prac_mm, y = results_train_consult2,
#           all = TRUE, by = "patid")
# a$rowSource <- apply(a[c("source.x", "source.y")], 1, 
#                      function(x) paste(na.omit(x), collapse = ""))
# head(a)

# medians for factor and ordered classes#
quantile.ordered <- function(x, probs = seq(0, 1, 0.25)) {
  tab <- table(x)
  cdf <- cumsum(tab / sum(tab))
  idx <- sapply(probs, function(p) min(which(cdf >= p)))
  levels(x)[idx] 
}

quantile.factor <- quantile.ordered
median.ordered <- function(x) quantile(x, 0.5)
median.factor <- median.ordered


detachAllPackages <- function(keep = NULL, keep.basic = TRUE) {
  # function for detaching all attached packages (except basic ones)
  basic.packages <- c("package:stats","package:graphics","package:grDevices",
                      "package:utils","package:datasets","package:methods",
                      "package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1,
                                  TRUE, FALSE)]
  if (!is.null(keep)){
    package.list <- setdiff(package.list, paste("package", keep, sep = ":"))
  }
  if (keep.basic){
    package.list <- setdiff(package.list, basic.packages)
  }
  if (length(package.list) > 0) {
    for (package in package.list) detach(package, character.only = TRUE)
  }
}
######################################################################################################
# # a (not yet, as not neccesary) function to compute classification matrix for validation #
# # est=estimated conditional prob. in the training set.
# # file = lca results.gh5 file
# 
#   a<-cbind(colmorbidlist,as.data.frame(mplus.get.estimated_probabilities("H:\\Mplus_MM\\unsup_patidsplit\\lca-unsup_8.gh5","process1",2,2)))
#   colnames(a) <- c("Morbidity","class 1:","class 2:","class 3:","class 4:","class 5:","class 6:","class 7:","class 8:")
#   a <- merge(a,prev[,c('Morbidity','Morbidityfull')],by='Morbidity',all.x=TRUE)
#   b<-test_pat_mm[,c(colmorbidlist,"patid")]
#   # order b columns (test data) so that comorb align with estimation matrix
#   order_morbid<-as.character(unlist(a[,1]))
#   #order training set class results by clslabel
#   classsize_cv<-classsize_unsup2[order(classsize_unsup2$clslabel),]
#   
#   df_cv<-data.frame(matrix(ncol=8,nrow=dim(b)[1])) # empty dataframe to put values in
#   b<-b[,c(order_morbid,"patid")]
#   bcopy<-b
#   for (n in 1:dim(bcopy)[1]){
#     print (n)
#     prob<-1
#     for (j in 1:8) { # 8 cls
#       prob<-1
#       for (i in 1:dim(a)[1]){ #38 dim
#         prob<-prob*as.numeric(a[i,j+1])^as.numeric(bcopy[n,i]-1)
#       }
#       c<-prob*classsize_cv$clsprop[j]/100
#       df_cv[n,j]<-c
#     }
#     df_cv[n,]<-df_cv[n,]/apply(df_cv[n,],1,sum)
#   }
#   
#   df_cv$modal<-max.col(df_cv,"first")
#   options(scipen = 999) #turn off scientific
#   df_cv
  
  # compute class proportions in the test set (hope that it is simliar to that in the training set)
##################################################################################################################


# Transpose a data frame # 
# transpose
# extra2<-extra %>%
#   t() %>%
#   as.data.frame(stringsAsFactors = F) %>%
#   rownames_to_column("value") %>%
#   `colnames<-`(.[1,]) %>%
#   .[-1,] %>%
#   `rownames<-`(NULL)


# # write latent gold syntax for 3-step misclassification correction # 
# / on-going /
# makeNewSyntax = function(dataDir,
#                          itemNames, #The names of the indicators
#                          mLevels, #A character vector being either "ordinal" or "continuous" to indicate the measurement level of each variable
#                          sets = 16, #number of sets of initial values
#                          iterations = 50){
#   
#   newSyntaxToBe = utils::capture.output(cat(paste("
#                                                   //LG5.0//
#                                                   version = 5.0
#                                                   infile '", dataDir,"'
#                                                   
#                                                   model
#                                                   options
#                                                     step3 modal ml;
#                                                     maxthreads=5;
#                                                   algorithm
#                                                     tolerance=1e-008 emtolerance=0,01 emiterations=250 nriterations=50;
#                                                   startvalues
#                                                     seed=0 sets=", sets," tolerance=1e-005 iterations=", iterations,";
#                                                   bayes
#                                                     categorical=1 variances=1 latent=1 poisson=1;
#                                                   montecarlo
#                                                     seed=0 sets=0 replicates=500 tolerance=1e-008;
#                                                   quadrature  
#                                                     nodes=10;
#                                                   missing  
#                                                     includeall;
#                                                   output
#                                                     parameters=first estimatedvalues=model reorderclasses
#                                                   write='C:\Users\zhuy18\Desktop\pj2\2lv_2distal_re3\high2000\esthigh2000_[[rep]].txt';
#                                                   variables
#                                                   caseweight ", weight,";
#                                                   dependent;
#                                                   latent
#                                                   Cluster nominal 1;
#                                                   equations
#                                                   Cluster <- 1;
#                                                   end model
#                                                   ")))
# 
# 
#   newSyntaxToBe[grep("dependent", newSyntaxToBe)] =
#     utils::capture.output(cat(paste0("   dependent ", paste(itemNames, mLevels, collapse = ", "), ";", sep = "")))
# 
#   newSyntaxToBe[length(newSyntaxToBe)] = paste0("   ", itemNames[1]," <- 1 + Cluster;")
#   for(idxVar in 2:length(itemNames)){
#     newSyntaxToBe[length(newSyntaxToBe) + 1] = paste0("   ", itemNames[idxVar]," <- 1 + Cluster;")
#   }
# 
#   for(idxVarCon in which(mLevels == "continuous")){
#     newSyntaxToBe[length(newSyntaxToBe) + 1] = paste0("   ", itemNames[idxVarCon]," | Cluster;")
#   }
# 
#   newSyntaxToBe[length(newSyntaxToBe) + 1] = "end model"
#   newSyntaxToBe[length(newSyntaxToBe) + 1] = ""
#   return(newSyntaxToBe)
# }
# ************************************************************************************************ #

# legacy functions, not needed
surv_weibull <- function(intercept,log_scale,t){exp(-(exp((log(t)-intercept)/exp(log_scale) )))}

pred_weibull <- function(scale,coef,formula,data,t){
  
  mat <- model.matrix(formula,data= data)
  
  surv_weibull(coef%*%t(mat),log(scale),t)
  
}

pred_coxnet <- function(b,data){
  
  X_new <- model.matrix(fit$formula,data= data)
  
  X_b <- b%*%t(data[,names(b)])
  
  bh <- basehaz(fit)
  
  five_years <- dim(bh)[1]
  
  exp(-bh[five_years,1])^exp(X_b)
  
  #plot(bh[,2],exp(-bh[,1]))
  #lines(bh[,2],exp(-bh[,1])^exp(X_b[1]),col='blue')
}


