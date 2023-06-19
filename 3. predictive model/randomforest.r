library(survminer)
library(randomForest)
library(survival)
library(broom)
library(timeROC)
library(caret)
library(ggpubr)
#devtools::install_github("selva86/InformationValue")
library(InformationValue)
#install.packages("ISLR")
library(ISLR)

#DMP showed BF significance in the discovery cohort
cgid=c('cg04893281','cg05773425','cg20005350','cg03546163','cg22282161','cg09610644','cg12234768','cg22675732','cg20728857','cg16500036','cg14911380','cg03775956','cg24447116','cg20015729','cg20783780','cg02584791','cg02869235','cg01068906','cg18017575','cg07623341','cg03237401')
train = read.table("discovery-death1-fdr0.05-residuals.txt",header = T,sep = "\t",row.names = 1)
train_pheno = read.table("20220414-405-pheno.txt",header = T,sep = "\t",row.names = 1)
train = cbind(train[,cgid],train_pheno$death1,train_pheno$death_time1)
#newtrain = cbind(train[,cgid],train_pheno$death1,train_pheno$death_time1,train_pheno$age,train_pheno$sex)
colnames(train)=c(cgid,"death","death_time")
#colnames(newtrain)=c(cgid,"death","death_time","age","sex")

######################################################random forest##############################################
train$death = as.factor(train$death)
train=train[,1:22]
##randomForest

#set.seed(5533)
#rf = randomForest(death~.,data=train)

set.seed(5533)
min=100
num=0
n<-length(names(train))
for (i in 1:(n-1)){
  mtry_fit<- randomForest(death~., data=train, mtry=i)
  err<-mean(mtry_fit$err.rate)
  print(err)
  if(err<min) {    
    min =err     
    num=i }
}
#min
#num

#select mtry with the lowest err
set.seed(5533)
rf = randomForest(death~.,
                  data=train,
                  mtry = 7,
                  ntree = 1000,
                  importance=TRUE)

varImpPlot(rf,
           sort = T,
           n.var = 21,
           main = "Variable Importance")
