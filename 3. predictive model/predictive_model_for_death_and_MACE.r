##predictive model with 5-fold 1000 times cross validation
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
library(caret)
library(timeROC)

####1. data import
train_pheno = read.table("20220518-405-pheno.txt",header = T,sep = "\t",row.names = 1)
train = read.table("discovery-death1-fdr0.05-residuals.txt",header = T,sep = "\t",row.names = 1)
cgid=c('cg04893281','cg05773425','cg20005350','cg03546163','cg22282161','cg09610644','cg12234768','cg22675732','cg20728857','cg16500036','cg14911380','cg03775956','cg24447116','cg20015729','cg20783780','cg02584791','cg02869235','cg01068906','cg18017575','cg07623341','cg03237401')
newtrain = cbind(train[,cgid],train_pheno$death1,train_pheno$death_time1,train_pheno$age,train_pheno$sex,train_pheno$fibrinogen,train_pheno$HDLC,train_pheno$LVEF,train_pheno$SII,train_pheno$PLR,train_pheno$LVMI,train_pheno$CHOL,train_pheno$WBC,train_pheno$LMR,train_pheno$NLR,train_pheno$LDLC,train_pheno$TRIG)
#colnames(train)=c(cgid,"death","death_time")
colnames(newtrain)=c(cgid,"death","death_time","age","sex","FIB","HDLC","LVEF","SII","PLR","LVMI","CHOL","WBC","LMR","NLR","LDLC","TRIG")

####2. construct function for 5-fold 1000 times cross validation
##function
tRocFuction_cv=function(data=null,formula_roc=null){
  for(i in 1:5000){
    train2<- data[folds[[i]],] #folds[[i]] as train datasets
    test <- data[-folds[[i]],] #the remaining as tedst datasets
    model = coxph(formula_roc,data=train2)
    model_pre<-predict(model,type='risk', newdata=test)
    tROC <-timeROC(T=test$death_time,delta = test$death,marker = model_pre,
                   cause = 1,times = c(1,3,5),ROC=T)
    auc_value<- rbind(auc_value,tROC$AUC)
  }
  return(auc_value)
}

####3. calculate auc value
set.seed(1234)
folds <-createMultiFolds(y=Surv(newtrain$death_time,newtrain$death),k=5,times=1000)

##3.1 traditional risk factors model
auc_value<-as.numeric(c(1,3,5))
auc<-as.numeric(c(1,3,5,7))
traditional_factor = c("age","sex","HDLC","LDLC","TRIG","CHOL","FIB","LVEF","LVMI","WBC","NLR","PLR","SII","LMR")

for (j in traditional_factor) {
    model = tRocFuction_cv(data = newtrain,formula_roc = Surv(death_time,death)~get(j))
    model = model[-1,]
    model_auc = colMeans(model)
    model_auc = cbind (model_auc,j)
    auc = rbind(auc,model_auc)
}
auc = auc[-1,]
colnames(auc) = c("1year","3years","5years","model_name")

##3.2 10 cg model
auc_value<-as.numeric(c(1,3,5))
cg_model = tRocFuction_cv(data = newtrain,formula_roc = Surv(death_time,death)~cg07623341+cg03237401+cg03546163+cg20728857+cg01068906+cg16500036+cg20015729+cg20005350+cg03775956+cg20783780)
cg_model = cg_model[-1,]
cg_model_auc = colMeans(cg_model)
cg_model_auc = append(cg_model_auc,"cgmodel")

##3.3 sex + age model
auc_value<-as.numeric(c(1,3,5))
sexage_model = tRocFuction_cv(data = newtrain,formula_roc = Surv(death_time,death)~sex+age)
sexage_model = sexage_model[-1,]
sexage_model_auc = colMeans(sexage_model)
sexage_model_auc = append(sexage_model_auc,"sexage")

##3.4 sex + age + cg model
auc_value<-as.numeric(c(1,3,5))
cgsexage_model = tRocFuction_cv(data = newtrain,formula_roc = Surv(death_time,death)~sex+age+cg07623341+cg03237401+cg03546163+cg20728857+cg01068906+cg16500036+cg20015729+cg20005350+cg03775956+cg20783780)
cgsexage_model = cgsexage_model[-1,]
cgsexage_model_AUC = colMeans(cgsexage_model)
cgsexage_model_AUC = append(cgsexage_model_AUC,"cgsexage")

##3.5 ensemble model
auc_value<-as.numeric(c(1,3,5))
ensemble_model = tRocFuction_cv(data = newtrain,formula_roc = Surv(death_time,death)~HDLC+age+sex+FIB+LVEF+cg07623341+cg03237401+cg03546163+cg20728857+cg01068906+cg16500036+cg20015729+cg20005350+cg03775956+cg20783780)
ensemble_model = ensemble_model[-1,]
ensemble_model_auc = colMeans(ensemble_model)
ensemble_model_auc = append(ensemble_model_auc,"ensemble")

####4. merge all results and save file
death_auc = as.data.frame(rbind(auc,cg_model_auc,sexage_model_auc,cgsexage_model_AUC,ensemble_model_auc))
death_auc$`1year` = as.numeric(death_auc$`1year`)
death_auc$`3years` = as.numeric(death_auc$`3years`)
death_auc$`5years` = as.numeric(death_auc$`5years`)

write.table(death_auc,"preditive_model_for_death_CV.txt",quote = F,sep = "\t")


### 5-fold 1000 times cross validation for MACE
##1. data import
train_pheno = read.table("20220518-405-pheno.txt",header = T,sep = "\t",row.names = 1)
train = read.table("discovery-death1-fdr0.05-residuals.txt",header = T,sep = "\t",row.names = 1)
cgid=c('cg05773425','cg03546163','cg20005350','cg16500036','cg18017575','cg12234768','cg13412579','cg24646285')
newtrain = cbind(train[,cgid],train_pheno$MACE,train_pheno$MACE_time,train_pheno$age,train_pheno$sex,train_pheno$fibrinogen,train_pheno$HDLC,train_pheno$LVEF,train_pheno$SII,train_pheno$PLR,train_pheno$LVMI,train_pheno$LDLC,train_pheno$HDLC,train_pheno$WBC,train_pheno$LMR,train_pheno$NLR,train_pheno$CHOL,train_pheno$TRIG)
#colnames(train)=c(cgid,"death","death_time")
colnames(newtrain)=c(cgid,"MACE","MACE_time","age","sex","FIB","HDLC","LVEF","SII","PLR","LVMI","LDLC","HDLC","WBC","LMR","NLR","CHOL","TRIG")

##2. function
tRocFuction_cv=function(data=null,formula_roc=null){
  for(i in 1:5000){
    train2<- data[folds[[i]],]
    test <- data[-folds[[i]],]
    model = coxph(formula_roc,data=train2)
    model_pre<-predict(model,type='risk', newdata=test)
    tROC <-timeROC(T=test$MACE_time,delta = test$MACE,marker = model_pre,
                   cause = 1,times = c(1,3,5),ROC=T)
    auc_value<- rbind(auc_value,tROC$AUC)
  }
  return(auc_value)
}
##3. construct model
set.seed(1234)
folds <-createMultiFolds(y=Surv(newtrain$MACE_time,newtrain$MACE),k=5,times=1000)

##3.1 traditional risk factors model
auc_value<-as.numeric(c(1,3,5))
auc<-as.numeric(c(1,3,5,7))
traditional_factor = c("age","sex","HDLC","LDLC","TRIG","CHOL","FIB","LVEF","LVMI","WBC","NLR","PLR","SII","LMR")

for (j in traditional_factor) {
    model = tRocFuction_cv(data = newtrain,formula_roc = Surv(MACE_time,MACE)~get(j))
    model = model[-1,]
    model_auc = colMeans(model)
    model_auc = cbind (model_auc,j)
    auc = rbind(auc,model_auc)
}
auc = auc[-1,]
colnames(auc) = c("1year","3years","5years","model_name")

##10 cg model
auc_value<-as.numeric(c(1,3,5))
cg_model = tRocFuction_cv(data = newtrain,formula_roc = Surv(MACE_time,MACE)~cg05773425+cg03546163+cg20005350+cg16500036+cg18017575+cg12234768+cg13412579+cg24646285)
cg_model = cg_model[-1,]
cg_model_auc = colMeans(cg_model)
cg_model_auc = append(cg_model_auc,"cgmodel")

##sex + age model
auc_value<-as.numeric(c(1,3,5))
sexage_model = tRocFuction_cv(data = newtrain,formula_roc = Surv(MACE_time,MACE)~sex+age)
sexage_model = sexage_model[-1,]
sexage_model_auc = colMeans(sexage_model)
sexage_model_auc = append(sexage_model_auc,"sexage")

##sex + age + cg model
auc_value<-as.numeric(c(1,3,5))
cgsexage_model = tRocFuction_cv(data = newtrain,formula_roc = Surv(MACE_time,MACE)~sex+age+cg05773425+cg03546163+cg20005350+cg16500036+cg18017575+cg12234768+cg13412579+cg24646285)
cgsexage_model = cgsexage_model[-1,]
cgsexage_model_AUC = colMeans(cgsexage_model)
cgsexage_model_AUC = append(cgsexage_model_AUC,"cgsexage")

##ensemble model
auc_value<-as.numeric(c(1,3,5))
ensemble_model = tRocFuction_cv(data = newtrain,formula_roc = Surv(MACE_time,MACE)~HDLC+age+sex+FIB+LVEF+cg05773425+cg03546163+cg20005350+cg16500036+cg18017575+cg12234768+cg13412579+cg24646285)
ensemble_model = ensemble_model[-1,]
ensemble_model_auc = colMeans(ensemble_model)
ensemble_model_auc = append(ensemble_model_auc,"ensemble")

MACE_auc = as.data.frame(rbind(auc,cg_model_auc,sexage_model_auc,cgsexage_model_AUC,ensemble_model_auc))
write.table(MACE_auc,"preditive_model_for_MACE_CV.txt",quote = F,sep = "\t")
