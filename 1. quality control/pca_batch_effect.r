library(minfi)
library(psych)
library(pheatmap)
library(plyr)
library(circlize)
library(ComplexHeatmap)
library(tidyr)

###1. get the positive control probe intensity
load("rgSet.Rdata")
a = getProbeInfo(rgSet,type = c("Control"))
a = a[a$Type!="NEGATIVE",]
id=as.character(a$Address)
control_matrix=rgSet[id,]
##get the methylated probe intensity
Meth = getGreen(control_matrix)
##get the unmethylated probe intensity
UnMeth = getRed(control_matrix)
Intensity=t(Meth+UnMeth)
write.table(Intensity,"positive_control_probe_intensity.txt",quote=F,sep="\t")

###2. perform the PCA for positive probe intensity matrix
#Intensity=t(Intensity)
#c.pca=prcomp(Intensity,scale=T)
#out=c.pca$x
#write.table(out,"positive_control_probe_pca.txt",quote=F,sep="\t")
#summary(c.pca)

###2. pca analysis for 405 CAD patients
pheno = read.csv("20180124_jda_517s_MethylationEPIC_Sample_Sheet.csv",header = T,sep = ",",row.names = 1)
batch1=pheno[pheno$Operter=="PYRJ",]
Intensity=read.table("positive_control_probe_intensity.txt",header = T,sep = "\t",row.names = 1)
Intensity=t(Intensity)
inten1=Intensity[rownames(batch1),]
c.pca=prcomp(inten1,scale=T)
out=c.pca$x
write.table(out,"positive_control_probe_pca_405.txt",quote=F,sep="\t")

##PCA plot
#library(factoextra)
#pdf("dis_PCA_plot.pdf",width = )
#fviz_pca_ind(c.pca, label="none", habillage=batch1$Sample_Group,
#             addEllipses=TRUE, ellipse.level=0.95,
#             palette = c("#999999", "#E69F00"))

##3. execute kruskal.test between batch and PCA of control probe
pca=read.table("positive_control_probe_pca_405.txt",header = T,sep = "\t",row.names = 1)
batch = read.csv("20180124_jda_517s_MethylationEPIC_Sample_Sheet.csv",header = T,sep = ",")
batch$Sentrix_ID = as.character(batch$Sentrix_ID)
batch1 = batch[batch$Operter=="PYRJ",]
pc20 = pca[,1:20]
dat = cbind(pc20,batch1)

Kruskal_test_model<- function(iv,data,dv){
  #
  FML<-as.formula(paste(iv,paste("~",dv,sep = "")))
  #KrusKal test
  res = kruskal.test(FML,data=data)
  #
  chi2 = res$statistic
  #ֵ
  P<-res$p.value
  #
  Kruskal_test_model <- data.frame('Characteristics'=iv,
                                   'chi2'=chi2,
                                   'P' = P)#[-1,]
  #                     
  return(Kruskal_test_model)
}

pc = paste0("PC",1:20,sep="")
covar = c("Sentrix_ID","Sentrix_Position","Array_Plate","Well")

out = c(1:3)
for (type in covar) {
  Uni_kw = lapply(pc, function(x) Kruskal_test_model(iv = x,data = dat,dv = type))
  a = ldply(Uni_kw,data.frame)
  a$type = type
  out = rbind(out,a[,c(1,3:4)])
}
out = out[-1,]
out2 = spread(out,key="Characteristics",value = "P",drop = FALSE)
rownames(out2) = out2[,1]
out2=out2[,-1]
out2=-log10(out2)

range(out2)

col_fun <- colorRamp2(seq(from = 0, to = 50, length = 9), 
                      RColorBrewer::brewer.pal(name = "YlGn", n = 9))
Heatmap(out2,
        col = col_fun,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_order = as.character(unique(out$Characteristics)),
        row_order = as.character(unique(out$type)),
        border = "black"
        )
