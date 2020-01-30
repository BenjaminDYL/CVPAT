CVPAT <- function(MV, CVFolds, Model1, Model2, testtype, BootSamp, boot.Di=FALSE, seed=FALSE, scale = TRUE){
  if (seed==FALSE) {
  }
  else{
    old <- .Random.seed
    on.exit( { .Random.seed <<- old } )
    set.seed(seed)
  }
  
  AllResults <- list(boot.p.values = NULL, losses = NULL, t.stat = NULL, p.value = NULL, conf.int = NULL)
  N <- nrow(MV)
  if (scale == TRUE) {
    MV <- scale(MV, center = TRUE, scale = TRUE)
  }
  #Bootstrap on D_i's
  if (boot.Di == TRUE) {
    Losses <- LossFun(N,MV,CVFolds,Model1$inner,Model1$reflective,Model1$formative,Model2$inner,Model2$reflective,Model2$formative)
    t.stat <- t.test(Losses$LossM2,Losses$LossM1,alternative = testtype, paired=TRUE)$statistic
    if (BootSamp==0) {
      Boot_di <- rep(NaN,5)
      names(Boot_di) <- c("OrgTtest.t", "p.value_perc_Ttest", "tstat_boot_Var", "p.value_var_ttest", "p.value_perc_D")
    }
    if (BootSamp>0) {
      Boot_di<-Bootstrap_di(Losses$LossM1, Losses$LossM2, BootSamp, testtype, N)
    }
    # Output results from bootstrap on D_i's
    AllResults[[1]] <- cbind("p.value.perc.t"=Boot_di["p.value_perc_Ttest"],
                             "p.value.b.v.t"=Boot_di["p.value_var_ttest"],
                             "p.value.perc.D"=Boot_di["p.value_perc_D"])
    AllResults[[2]] <- Losses
    AllResults[[3]] <- cbind("t.stat"=t.stat,"t.stat.b.v"=Boot_di["tstat_boot_Var"])
    AllResults[[4]] <- t.test(Losses$LossM2,Losses$LossM1,alternative = testtype, paired=TRUE)$p.value
    AllResults[[5]] <- t.test(Losses$LossM2,Losses$LossM1,alternative = testtype, paired=TRUE)$conf.int
    
  }
  # Bootstrap on MV
  if (boot.Di == FALSE) {
    Losses <- LossFun(N,MV,CVFolds,Model1$inner,Model1$reflective,Model1$formative,Model2$inner,Model2$reflective,Model2$formative)
    mv.org.ttest<-t.test(Losses$LossM2,Losses$LossM1,mu=0,alternative=testtype, paired=TRUE)$statistic
    Di <- Losses$LossM2 - Losses$LossM1
    D.bar <- mean(Di)
    if (BootSamp==0) {
      boot_MV <- rep(NaN, 5)
      names(boot_MV) <- c("p.value.perc.t", "p.value_var_ttest", "p.value_perc_D",
                          "t.stat.b.v", "prop.non.conv")
    }
    if (BootSamp>0) {
      boot_MV <- Bootstrap_MV(mv.org.ttest, Di, D.bar, N,MV,CVFolds,InnerSpecM1,ReflecSpecM1,FormSpecM1,InnerSpecM2,ReflecSpecM2,FormSpecM2, BootSamp)
    }
    # Save results from bootstrapping on MV
    
    AllResults[[1]] <- cbind("p.value.perc.t" = boot_MV["p.value.perc.t"], "p.value.b.v.t" = boot_MV["p.value_var_ttest"],
                             "p.value.perc.D" = boot_MV["p.value.perc.D"])
    AllResults[[2]] <- Losses
    AllResults[[3]] <- cbind("t.stat"=mv.org.ttest,"t.stat.b.v"=boot_MV["t.stat.b.v"])
    AllResults[[4]] <- t.test(Losses$LossM2,Losses$LossM1,alternative = testtype, paired=TRUE)$p.value
    AllResults[[5]] <- t.test(Losses$LossM2,Losses$LossM1,alternative = testtype, paired=TRUE)$conf.int
    AllResults$conv.fail <- boot_MV["prop.non.conv"]
    
    
  }
  
  AllResults
}

# Helpers -----------------------------------------------------------------
LossFun <- function(N,MV,CVFolds,InnerSpecM1,ReflecSpecM1,FormSpecM1,InnerSpecM2,ReflecSpecM2,FormSpecM2) {
  # Relation between MV and LV
  MVLVrelM1<-ReflecSpecM1+t(FormSpecM1)
  MVLVrelM2<-ReflecSpecM2+t(FormSpecM2)
  # Specifying model in matrix-PLS format
  M1<-list(inner = InnerSpecM1,reflective = ReflecSpecM1,formative = FormSpecM1)
  M2<-list(inner = InnerSpecM2,reflective = ReflecSpecM2,formative = FormSpecM2)
  ## Allocating memory to mean of squared errors and estimated parameters
  LossM1 <- matrix(rep(0,N),nrow = N, ncol = 1)
  LossM2 <- matrix(rep(0,N),nrow = N, ncol = 1)
  endogLVM1<-colnames(ReflecSpecM1)[!colnames(ReflecSpecM1)%in%colnames(ReflecSpecM1)[which(rowSums(InnerSpecM1)==0)]]
  LossM1_sepLV<- setNames(vector("list",length(endogLVM1)),endogLVM1)
  LossM1_sepLV<-lapply(LossM1_sepLV, function(x) x<-matrix(rep(0,N),nrow = N, ncol = 1))
  endogLVM2<-colnames(ReflecSpecM2)[!colnames(ReflecSpecM2)%in%colnames(ReflecSpecM2)[which(rowSums(InnerSpecM2)==0)]]
  LossM2_sepLV<- setNames(vector("list",length(endogLVM2)),endogLVM2)
  LossM2_sepLV<-lapply(LossM2_sepLV, function(x) x<-matrix(rep(0,N),nrow = N, ncol = 1))
  ## Split data into folds
  folds<-cvFolds(N,CVFolds,1)
  folds<-cbind(folds$which,folds$subsets)
  for (i in 1:CVFolds) {
    testIndexes<- folds[folds[,1]==i,2]
    TempMV     <- MV[-testIndexes,]
    ## Estimate PLS
    TempPLSM1 <- matrixpls(cov(TempMV),model = M1)
    TempPLSM2 <- matrixpls(cov(TempMV),model = M2)
    ## Save results from estimation
    EstLoadingsM1<-attr(TempPLSM1, "reflective")
    EstWeightsM1<-t(attr(TempPLSM1, "W"))
    EstInnerM1<-attr(TempPLSM1, "inner")
    EstLoadingsM2<-attr(TempPLSM2, "reflective")
    EstWeightsM2<-t(attr(TempPLSM2, "W"))
    EstInnerM2<-attr(TempPLSM2, "inner")
    ## Calculate prediction errors
    # Calculate proxy of LV from indicators
    FoldSize <- length(testIndexes)
    if (FoldSize>1) {
      LVProxyM1<-t(EstWeightsM1)%*%t(MV[testIndexes,])
      LVProxyM2<-t(EstWeightsM2)%*%t(MV[testIndexes,])
    } else {
      LVProxyM1<-t(EstWeightsM1)%*%MV[testIndexes,]
      LVProxyM2<-t(EstWeightsM2)%*%MV[testIndexes,]
    }
    # Predict endog LV from exog LV
    PredictLVM1<-EstInnerM1%*%LVProxyM1
    PredictLVM2<-EstInnerM2%*%LVProxyM2
    # Predict MV connected to endogenous LV
    PredictMVM1<-EstLoadingsM1%*%PredictLVM1
    PredictMVM2<-EstLoadingsM2%*%PredictLVM2
    ## prediction errors
    PseudoPredErrorsM1<- MV[testIndexes,]-t(PredictMVM1)
    PseudoPredErrorsM2<- MV[testIndexes,]-t(PredictMVM2)
    # Identifying predicition errros connected to endogenous LV.
    PureExogLVM1<-which(rowSums(InnerSpecM1)==0)
    if (length(PureExogLVM1)==1){
      PureExogMVM1<-which(MVLVrelM1[,PureExogLVM1]!=0)
    } else {
      PureExogMVM1<-which(rowSums(MVLVrelM1[,PureExogLVM1])!=0)
    }
    PredErrorsM1<-as.matrix(PseudoPredErrorsM1[,-PureExogMVM1])
    PureExogLVM2<-which(rowSums(InnerSpecM2)==0)
    if (length(PureExogLVM2)==1){
      PureExogMVM2<-which(MVLVrelM2[,PureExogLVM2]!=0)
    } else {
      PureExogMVM2<-which(rowSums(MVLVrelM2[,PureExogLVM2])!=0)
    }
    PredErrorsM2<-as.matrix(PseudoPredErrorsM2[,-PureExogMVM2])
    ## Loss function for entire model
    if (FoldSize>1){
      LossM1[testIndexes,]<-rowMeans(PredErrorsM1^2)
      LossM2[testIndexes,]<-rowMeans(PredErrorsM2^2)
    } else {
      LossM1[testIndexes,]<-mean(PredErrorsM1^2)
      LossM2[testIndexes,]<-mean(PredErrorsM2^2)
    }
    ## Losses for each endogenous LV
    for (j in endogLVM1) {
      if (length(rownames(ReflecSpecM1)[ReflecSpecM1[,j]==1])>1) {
        if (FoldSize>1) {
          LossM1_sepLV[[j]][testIndexes]<-rowMeans(PseudoPredErrorsM1[,(rownames(ReflecSpecM1)[ReflecSpecM1[,j]==1])]^2)
        } else {
          LossM1_sepLV[[j]][testIndexes]<-mean(PseudoPredErrorsM1[,(rownames(ReflecSpecM1)[ReflecSpecM1[,j]==1])]^2)
        }
      }else{
        LossM1_sepLV[[j]][testIndexes]<-PseudoPredErrorsM1[,(rownames(ReflecSpecM1)[ReflecSpecM1[,j]==1])]^2
      }
    }
    
    for (j in endogLVM2) {
      if (length(rownames(ReflecSpecM2)[ReflecSpecM2[,j]==1])>1) {
        if (FoldSize>1) {
          LossM2_sepLV[[j]][testIndexes]<-rowMeans(PseudoPredErrorsM2[,(rownames(ReflecSpecM2)[ReflecSpecM2[,j]==1])]^2)
        } else {
          LossM2_sepLV[[j]][testIndexes]<-mean(PseudoPredErrorsM2[,(rownames(ReflecSpecM2)[ReflecSpecM2[,j]==1])]^2)
        }
      }else{
        LossM2_sepLV[[j]][testIndexes]<-PseudoPredErrorsM2[,(rownames(ReflecSpecM2)[ReflecSpecM2[,j]==1])]^2
        
      }
    }
  }
  list(LossM1 = LossM1, LossM2 = LossM2, LossM1_sepLV = LossM1_sepLV, LossM2_sepLV = LossM2_sepLV )
}
# Bootstrapping the t-test
Bootstrap_di <- function (LossM1,LossM2,BootSamp,testtype, N){
  # Originial test-statistic
  OrgTtest<-t.test(LossM2,LossM1,alternative = testtype, paired=TRUE)$statistic
  # Originial average difference in losses
  OrgDbar<-mean(LossM2-LossM1)
  # Differences in loss functions under the null
  D_0<-LossM2-LossM1-OrgDbar
  # Differences in loss functions
  D<-LossM2-LossM1
  #Allocating memory to bootrap
  BootSample <- matrix(0,ncol=2,nrow=length(D))
  BootDbar <- rep(0,BootSamp)
  m_losses<-cbind(LossM1,LossM2)
  tStat <- rep(0,BootSamp)
  for (b in 1:BootSamp) {
    BootSample <- m_losses[sample((1:length(D)), length(D), replace=TRUE),]
    tStat[b]<-t.test(BootSample[,2],BootSample[,1],mu=mean(D),alternative=testtype, paired=TRUE)$statistic
    BootDbar[b] <- mean(sample(D_0, length(D_0), replace=TRUE))
  }
  SorttStat<-sort(tStat, decreasing = FALSE)
  SortBootDbar<-sort(BootDbar, decreasing = FALSE)
  # Bootstrap variance on Dbar for t-test
  std<-sqrt(var(BootDbar))
  tstat_boot_Var<-OrgDbar/std
  # Calculating p-values
  if (testtype=="two.sided") {
    p.value_perc_Ttest<-(sum(SorttStat>abs(OrgTtest))+sum(SorttStat<=(-abs(OrgTtest))))/BootSamp
    p.value_perc_D<-(sum(SortBootDbar>abs(OrgDbar))+sum(SortBootDbar<=(-abs(OrgDbar))))/BootSamp
    p.value_var_ttest<-2*pt(-abs(tstat_boot_Var),(N-1), lower.tail = TRUE)
  }
  if (testtype=="greater") {
    p.value_perc_Ttest<-1-(head(which(SorttStat>OrgTtest),1)-1)/(BootSamp+1)
    p.value_perc_D<-1-(head(which(SortBootDbar>OrgDbar),1)-1)/(BootSamp+1)
    p.value_var_ttest<-pt(tstat_boot_Var,(N-1), lower.tail = FALSE)
    if (length(which(SorttStat>OrgTtest))==0){
      p.value_perc_Ttest=0
      p.value_perc_D=0
    }
  }
  Results <- c("OrgTtest"=OrgTtest, "p.value_perc_Ttest"=p.value_perc_Ttest,
               "tstat_boot_Var" = tstat_boot_Var, "p.value_var_ttest" = p.value_var_ttest,
               "p.value_perc_D" = p.value_perc_D)
  return(Results)
}
# Bootstrap on MV
Bootstrap_MV <- function(mv.org.ttest, Di, D.bar, N,MV,CVFolds,InnerSpecM1,ReflecSpecM1,FormSpecM1,InnerSpecM2,ReflecSpecM2,FormSpecM2, BootSamp){
  b.t.stat <- rep(0,BootSamp)
  b.D.bar <- rep(0,BootSamp)
  for (b in 1:BootSamp) {
    b.samp<-sample(1:N,N,replace = TRUE)
    b.MV<-MV[b.samp,]
    b.Losses <- tryCatch(LossFun(N,b.MV,CVFolds,Model1$inner,Model1$reflective,Model1$formative,Model2$inner,Model2$reflective,Model2$formative),
                         warning = function(w) {NA},
                         error = function(e) {NA})
    if (is.na(b.Losses[1])) {
      b.D.bar[b] <- NA
      b.t.stat[b]<- NA
    } else {
      b.D.bar[b] <- mean(b.Losses$LossM2 - b.Losses$LossM1) - D.bar
      b.t.stat[b]<-t.test(b.Losses$LossM2,b.Losses$LossM1,mu=D.bar,alternative=testtype, paired=TRUE)$statistic
    }
  }
  prop.non.conv <- sum(is.na(b.D.bar))/length(b.D.bar)
  BootSamp<-length(b.D.bar[!is.na(b.D.bar)])
  b.D.bar<-b.D.bar[!is.na(b.D.bar)]
  b.t.stat<-b.t.stat[!is.na(b.t.stat)]
  # Sorted t-statistics
  t.sort <- sort(b.t.stat, decreasing = FALSE)
  t.stat.b.v <- D.bar/sqrt(var(b.D.bar))
  D.bar.sort <- sort(b.D.bar, decreasing = FALSE)
  if (testtype=="two.sided") {
    # Percentile-t method, p-value
    p.value.perc.t<-(sum(t.sort>abs(mv.org.ttest))+sum(t.sort<=(-abs(mv.org.ttest))))/BootSamp
    # Botstrapped variance of D.bar, p-value
    p.value_var_ttest<-2*pt(-abs(t.stat.b.v),(N-1), lower.tail = TRUE)
    # percentile D.bar, p-value
    p.value_perc_D<-(sum(D.bar.sort>abs(D.bar))+sum(D.bar.sort<=(-abs(D.bar))))/BootSamp
  }
  if (testtype=="greater") {
    # Percentile-t method, p-value
    p.value.perc.t<-1-(head(which(t.sort>mv.org.ttest),1)-1)/(BootSamp+1)
    if (length(which(t.sort>mv.org.ttest))==0){
      p.value.perc.t=0
    }
    # Botstrapped variance of D.bar, p-value
    p.value_var_ttest<-pt(t.stat.b.v,(N-1), lower.tail = FALSE)
    # percentile D.bar, p-value
    p.value_perc_D<-1-(head(which(D.bar.sort>D.bar),1)-1)/(BootSamp+1)
    if (length(which(D.bar.sort>D.bar))==0){
      p.value_perc_D=0
    }
  }
  return(c("p.value.perc.t" = p.value.perc.t, "p.value_var_ttest" = p.value_var_ttest, "p.value.perc.D" = p.value_perc_D,
           "t.stat.b.v" = t.stat.b.v, "prop.non.conv" = prop.non.conv))
}
