# Load required packages --------------------------------------------------
library(matrixpls)
library(cvTools)

# Load data ---------------------------------------------------------------
load(file="corp_rep.rda")

# matrixpls specification for competing models ----------------------------

var_names <- colnames(corp_rep)
# Model specifications, Proposed Model (PM)
CSOR = c(0,0,0,0,0,0,0,0)
ATTR = c(0,0,0,0,0,0,0,0)
PERF = c(0,0,0,0,0,0,0,0)
QUAL = c(0,0,0,0,0,0,0,0)
LIKE = c(1,1,1,1,0,0,0,0)
COMP = c(1,1,1,1,0,0,0,0)
CUSA = c(0,0,0,0,1,1,0,0)
CUSL = c(0,0,0,0,1,1,1,0)
inner_PM  = rbind(CSOR, ATTR, PERF, QUAL, LIKE, COMP, CUSA, CUSL)
colnames(inner_PM) <- rownames(inner_PM)

# Reflective measurement model
reflective<-matrix(0,length(var_names),ncol(inner_PM))
reflective[22:24,5] <-1
reflective[25:27,6] <-1
reflective[28,7]    <-1
reflective[29:31,8] <-1
dimnames(reflective) <- list(var_names,colnames(inner_PM))

# Formative measurement model
formative <- (reflective*0)
formative[1:5,1]   <-1
formative[6:8,2]   <-1
formative[9:13,3]  <-1
formative[14:21,4] <-1
formative <- t(formative)
# Model relations summarized in list for matrixpls package, PM
PM <- list(inner = inner_PM,reflective = reflective,formative = formative)

# Model specifications Alternative Model (AM)
CSOR = c(0,0,0,0,0,0,0,0)
ATTR = c(0,0,0,0,0,0,0,0)
PERF = c(0,0,0,0,0,0,0,0)
QUAL = c(0,0,0,0,0,0,0,0)
LIKE = c(1,1,1,1,0,0,0,0)
COMP = c(1,1,1,1,0,0,0,0)
CUSA = c(1,1,1,1,1,1,0,0)
CUSL = c(1,1,1,1,1,1,1,0)
inner_AM  = rbind(CSOR, ATTR, PERF, QUAL, LIKE, COMP, CUSA, CUSL)
colnames(inner_AM) <- rownames(inner_AM)

# Model relations summarized in list for matrixpls package, AM
AM <- list(inner = inner_AM,reflective = reflective,formative = formative)

# Model comparison with CVPAT ---------------------------------------------
# Load CVPAT function
source("CVPAT.R")
res_CVPAT <- CVPAT(MV=corp_rep,
                   CVFolds = 10,Model1=AM,Model2 = PM,
                   testtype=c("greater"),BootSamp = 2000,boot.Di = T,seed=TRUE,scale=TRUE)
# Average losses for each model (PM has slightly lower loss than AM)
res_CVPAT$losses$avg_losses$avg_losses_M1
res_CVPAT$losses$avg_losses$avg_losses_M2
# Bootstrapped p-values 
res_CVPAT$boot.p.values
# non-bootstrapped p-value
res_CVPAT$p.value
# non-bootstrapped confidence interval
res_CVPAT$conf.int
