#setwd("/home/qinmin/850K/discovery/517data/")
library(survival)
library(plyr)

argv<-commandArgs(TRUE)
beta=read.table(argv[1],header=T,sep="\t",row.names = 1)
pheno = read.table("/home/qinmin/Project/EWAS_OF_CAD/405_discovery/death/20220414-405-pheno.txt",header = T,sep="\t",row.names=1)
#cell = read.table("/home/qinmin/850K/discovery/517data/20220404_517_Cell_type.txt",header = T,sep="\t",row.names=1)
dat = cbind(pheno,scale(beta))
clinical_covar = c("age","sex","PCIorNot","smoking2","arrhythmia","HF","HyperT","Hyperlip","BB","ACEI","CCB","PPI","clop","Aspirin")
#pc_covar = paste0("PC",1:10,sep="")
#cell_covar = c("CD8T","CD4T","NK","Bcell","Mono","Gran")
#covar_name = c(clinical_covar,pc_covar,cell_covar)

##construct function
Multi_cox_model<- 
  function(x,covar,data){
    FML<-as.formula(paste("Surv(death_time1,death1)~ ", paste(c(x,covar), collapse= "+")))
    cox<-coxph(FML,data=data)
    cox1<-summary(cox)
    beta = round(cox1$coefficients[1,1],3)
    HR<-round(cox1$coefficients[1,2],2)
    SE<-round(cox1$coefficients[1,3],2)
    CI <- paste0(round(cox1$conf.int[1,3:4],2),collapse = '-')
    HRCI<-paste0(HR,'(',CI,')')
    P<-cox1$coefficients[1,5]
    Multi_cox_model <- data.frame('Characteristics'=x,
                                'beta'=beta,
                                'se'=SE,
                                'HR(95%CI)' = HRCI,
                                'P' = P)#[-1,]                     
    return(Multi_cox_model)
  }

iv_name = colnames(beta)
Multi_cox <- lapply(iv_name, function (x) Multi_cox_model(x = x, covar = clinical_covar,data = dat))

##merge the result
Multi_cox<- ldply(Multi_cox,data.frame)
write.table(Multi_cox,file = argv[2],quote = F,sep = "\t",row.names=F)
