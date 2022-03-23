################################################
#### Framework to run 3-step LASSO approach ####
################################################

## This script contains an example of the code used to run a 3-step LASSO framework for feature selection followed by optimization and performance evaluation 

## In this example the data has already been loaded and cleaned
## Brief description of the dataframe and variables used 
## fenland : dataframe containing all the required variables (both the proteome and phenotyp data)
## IGT : name of the variable encoding the phenotype to be predicted (coded as 0/1)
## prots : vector contraining the 4979 protein names as they are called in the fenland dataframe

## Packages required
library(caret)
library(glmnet)
library(doMC)
library(ROSE)
library(pROC)

######################################################
####   Split the data into training and testing   ####
######################################################

## Split data set in training (50%), optimization (25%) and testing sets (25%)
##for reproducibilty set seed
set.seed(3456)
trainIndex <- createDataPartition(fenland$IGT,p=.5,list = FALSE,times = 1) ## use variable to be used as the outcome to generate this split to ensure comparable number of cases between training and testing sets
#Get train set
fenland.Train <- fenland[trainIndex,]
#Get test set (has to further divided to be used as optimization and final test) 
fenland.Test  <- fenland[-trainIndex,]

## Split into optimization and final test
set.seed(3456)
trainIndex2<-createDataPartition(fenland.Test$IGT,p=.5,list = F,times = 1)

fenland.Opti<-fenland.Test[trainIndex2,]
fenland.finalTest<-fenland.Test[-trainIndex2,]

######################################################
####     Run feature selection on training set    ####
######################################################

## Create lists to store results
cf<-list()
las.morb<-list()

## Draw 100 random subsets of the training data (80% each time) 
jj <- lapply(1:100, function(x) sample(nrow(fenland.Train),round(nrow(fenland.Train)*0.80), replace = F)) ## subsampling approach

## kernel for parallel processing
registerDoMC(24)

## Loop to 
for (i in 1:length(jj)){
  print(i)
  las.morb[[i]] <- train(fenland.Train[jj[[i]],prots], fenland.Train[jj[[i]], "IGT"],
                         family="binomial",
                         method = "glmnet",
                         tuneGrid = as.data.frame(expand.grid(alpha=1, lambda=10^-seq(10,.25,-.25))),
                         trControl = trainControl(method="repeatedcv", number=10, repeats = 10, sampling="rose", ## weighting scheme to deal with case imbalance
                                                  allowParallel=T)
  )
  cf[[i]]       <- as.matrix(coef(las.morb[[i]]$finalModel, s = las.morb[[i]]$finalModel$lambdaOpt))
  las.morb[[i]] <- las.morb[[i]]$results
  gc()
}

## Save training results 
save(las.morb, file="../data_input/las.morb.IGT.RData")
save(cf, file="../data_input/cf.IGT.RData")

########################################################
####   Optimization using top informative features  ####
########################################################

# load("../data_input/las.morb.IGT.RData")
# load("../data_input/cf.IGT.RData")

## Generate feature selection ranking
## Mtrix of 0/1 to denote whether feature was selected in the final model during each round
prot.select<-matrix(, nrow = 4980,ncol = 100) ## one more row than proteins for the intercept

for (i in 1:100){
  temp<-cbind(ifelse(cf[[i]][,1]!=0,1,0))
  prot.select[,i]<-temp
  row.names(prot.select)<-row.names(temp)
}

##Count matrix 
lasso.count <- as.data.frame(prot.select)
lasso.count$select <- rowSums(lasso.count[,1:100])
lasso.count <- lasso.count[-1,] ## remove the intercept
lasso.count$MRC_seqid <- row.names(lasso.count)

##Sort in order of proteins most times selected across boots 
lasso.sort <- lasso.count[order(lasso.count$select,decreasing = T),]

## Example of optimization using proteins selected in more than 95% of rounds
prots.8 <-filter(lasso.sort,select>95)$MRC_seqid

## Hyperparameter optimization using proteins selected in more than 95% of rounds
las.morb.8 <- train(fenland.Opti[,prots.8], fenland.Opti[, "IGT"],
                    family="binomial",
                    method = "glmnet",
                    tuneGrid = as.data.frame(expand.grid(alpha=1, lambda=10^-seq(10,.25,-.25))),
                    trControl = trainControl(method="repeatedcv", number=10, repeats = 10, sampling="rose")) ## same weighting scheme to deal with case imbalance

#############################################
####   Test final model in the test set  ####
#############################################

## Get final models and optimal lambdas 
glmnet.8 <- las.morb.8$finalModel
s.8 <- las.morb.8$finalModel$lambdaOpt

## Look at predictors selected in the final optimized model
predictors(las.morb.8)

## Test preformance in the final test 
predict.8 <- as.numeric(predict(glmnet.8,type = "response",newx = as.matrix(fenland.finalTest[,prots.8]),s=s.8))

roc.8 <- roc(response = fenland.finalTest$IGT,predictor =  predict.8)
