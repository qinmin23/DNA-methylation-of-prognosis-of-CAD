#####################################################################Figure 2A circos plot
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(ArchR)
library(grDevices)
library(stringr)

dat = read.table("115_DMP_for_curcosplot.txt",header = T,sep = "\t")
dat$index = 1
dat[is.na(dat)] = 0
dat = dat[order(dat$CHR,dat$POS),]
#dat$death_beta = ifelse(dat$death>0,1,ifelse(dat$death<0,-1,0))
#dat$MACE_beta = ifelse(dat$MACE>0,1,ifelse(dat$MACE<0,-1,0))
dat$index2 = c(seq(1,13),seq(1,7),seq(1,5),seq(1,5),seq(1,7),seq(1,12),seq(1,6),seq(1,3),seq(1,5),seq(1,6),seq(1,2),seq(1,13),seq(1,3),seq(1,4),seq(1,4),seq(1,7),seq(1,6),1,1,1,seq(1,4))
dat$label = paste(dat$CpG,dat$Gene,sep = " ")
dat$chr2=paste("chr",dat$CHR,sep = "")
aa=dat[,c(3,5)]
aa=scale(aa)
rownames(aa)=dat$CpG

#barcols  <- c("firebrick","steelblue","goldenrod")

#split = unique(str_sort(dat$chr2, numeric = TRUE))
#split = factor(rownames(aa))
split = factor(str_sort(dat$chr2, numeric = TRUE),levels = unique(as.character(dat$chr2)))

##chromosome color
#cols <- ArchR::paletteDiscrete(unique(factor(dat$CHR)))
#cols = sample(colors(TRUE), 21)
color_pal <- c("#FF3200", "#E9A17C", "#E9E4A6", "#1BB6AF", "#0076BB", "#172869")
more_colors <- (grDevices::colorRampPalette(color_pal))(21)
scales::show_col(colours = more_colors)

#colline = rep(rev(viridis(21)),times=table(dat$CHR))
colline = rep(more_colors,times=table(dat$chr2))

col_fun1 = colorRamp2(c(-2,0,2),c("#1b9e77","white","#7570b3"))
#col_fun = c(c("#1b9e77","white","#7570b3"))

circos.par(cell.padding = c(0.01, 0, 0.01, 0),
           gap.after = c(rep(1, 20), 20),
           start.degree = 0,
           track.height = 0.25,
           track.margin = c(0.01,0.01),
           circle.margin = 0.45)

circos.heatmap.initialize(aa,cluster = FALSE,split = split)

#circos.initialize(factors = dat$CpG, 
#                  x = dat$index, 
#                  xlim = c(0,2))

circos.trackPlotRegion(factors = dat$chr2, 
                       y = abs(dat$MACE), 
                       bg.border = more_colors,
                       #bg.col = c(rep(c("white","#f0f0f0"),10),"white"),
                       panel.fun = function(x, y) {
                         circos.yaxis(side = "left", 
                                      at = seq(0, 10, by = 2),
                                      sector.index = get.all.sector.index()[1], 
                                      labels.cex = 0.4, 
                                      labels.niceFacing = FALSE)
                       })

circos.trackLines(factors = dat$chr2, 
                  x = dat$index2-0.5, 
                  y = abs(dat$MACE), 
                  pch=20, 
                  cex=2, 
                  col = colline,
                  type="h", 
                  lwd = 4)
col=viridis(10)

circos.heatmap(as.matrix(aa[,2]),col = col_fun1, track.height = 0.05,cell_width = 0.25,rownames.side = "none")

circos.trackText(factors = dat$chr2,x=dat$index2-1,y=rep(9,115),labels = dat$label,cex = 0.6,facing = "clockwise",col = colline,adj=rep(c(0,1),115))

circos.trackPlotRegion(factors = dat$chr2, 
                       y = abs(dat$death), 
                       bg.border = more_colors,
                       #bg.col = c(rep(c("white","#f0f0f0"),10),"white"),
                       panel.fun = function(x, y) {
                         
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         
                         #text direction (dd) and adjusmtents (aa)
                         theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                         dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                         aa = c(1, 0.5)
                         if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                         
                         #plot labels
                         circos.text(x=mean(xlim), y=-7.5, labels=name, facing = dd, cex=0.6,  adj = aa)
                         
                         #plot main sector
                         #circos.axis(labels=FALSE, major.tick=FALSE)
                         circos.yaxis(side = "left", at = seq(0, 10, by = 2),
                                      sector.index = get.all.sector.index()[1], labels.cex = 0.4, labels.niceFacing = FALSE)
                       }
                       )


circos.trackLines(factors = dat$chr2, 
                  x = dat$index2-0.5, 
                  y = abs(dat$death), 
                  pch=20, 
                  cex=2, 
                  col = colline,
                  type="h", 
                  lwd = 4)

#insert heatmap
circos.heatmap(aa[,1], col = col_fun1,track.height = 0.05)
circos.clear()

lgd_lines = Legend(at = unique(dat$CHR),
                   type = "lines",
                   background = "white",
                   legend_gp = gpar(col = more_colors, lwd = 2),
                   title_position = "lefttop",
                   labels_gp = gpar(fontsize = 14, lex = 4),
                   title = "Chromosome")
lgd_list_vertical = packLegend(lgd_lines)


pushViewport(viewport(x = unit(10, "mm"), y = unit(4, "mm"),
                      #width = grobWidth(lgd_list_vertical),
                      #height = grobHeight(lgd_list_vertical),
                      just = c("left", "top")))

grid.draw(lgd_list_vertical)
upViewport()

#######################################################################Figure 2B and 2C enrichment plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)

GWASCatalog_dataset <- read.table(file ="DMP_GWASCatalog_enrich.txt",header = TRUE, sep = "\t")
GOBP_dataset <- read.table(file ="DMP_GOBP_enrich.txt",header = TRUE, sep = "\t")

#sort data by Pvalue
GWASCatalog_dataset <- arrange(GWASCatalog_dataset,desc(GWASCatalog_dataset[,5]))
GOBP_dataset <- arrange(GOBP_dataset,desc(GOBP_dataset[,5]))

GWASCatalog_dataset$Term <- factor(GWASCatalog_dataset$Term,levels = rev(GWASCatalog_dataset$Term))
GOBP_dataset$Term <- factor(GOBP_dataset$Term,levels = rev(GOBP_dataset$Term))

mytheme <- theme(axis.title=element_text(face="bold", size=14,colour = 'black'), 
                 axis.text=element_text(face="bold", size=14,colour = 'black'), 
                 axis.line = element_line(size=0.5, colour = 'black'), 
                 panel.background = element_rect(color='black'),
                 legend.key = element_blank()
)

#
p <- ggplot(GOBP_dataset,aes(x=GeneRatio,y=Term,colour=-1*log10(P.value),size=Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradientn(colours = rev(rocket(20)))+
  theme_bw()+
  ylab("GO BP Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))
  
plot <- p+mytheme
plot

ggsave(plot,filename = "GWASCatalog.pdf",width = 12,height = 6,dpi=300)
ggsave(plot,filename = "GOBP.pdf",width = 12,height = 7,dpi=300)

#######################################################################Figure 2D heatmap
library(pheatmap)
library(RColorBrewer)
library(viridis)

coef=read.table("pheno_DMP_coef.txt",header = T,sep = "\t",row.names = 1)
coef=coef[,c(1,13,2,14,3,15,4,16,5,17,6,18,7,19,8,20,9,21,10,22,11,23,12,24)]
gene = read.table("115_DMP_for_curcosplot.txt",header = T,sep = "\t")
p = read.table("pheno_DMP_P.txt",header = T,sep = "\t",row.names = 1)
#p = p[,c(1,3,5,7,9,11,13,15,17,19,21,23,2,4,6,8,10,12,14,16,18,20,22,24)]

#beta=a[,c(1,3,5,7)]

annotation_col=data.frame(cohort=factor(rep(c("discovery cohort","validation cohort"),12)))
rownames(annotation_col)=colnames(coef)

#p=a[,c(2,4,6,8)]
if (!is.null(p)){
  ssmt <- p< 0.01
  p[ssmt] <-'**'
  smt <- p >0.01& p <0.05
  p[smt] <- '*'
  p[!ssmt&!smt]<- ''
} else {
  pmt <- F
}

#ann_colors = list(cohort = c('discovery cohort' = "#D7301F",
#                             'validation cohort'= "#225EA8"))

annotation_row = data.frame(class=factor(c(rep("WBC",2),
                                           rep("LMR",2),
                                           rep("NLR",2),
                                           rep("PLR",2),
                                           rep("FIB",2),
                                           rep("SII",2),
                                           rep("LDLC",2),
                                           rep("HDLC",2),
                                           rep("CHOL",2),
                                           rep("TRIG",2),
                                           rep("LVEF",2),
                                           rep("LVMI",2)
                                           )))
rownames(annotation_row)=colnames(coef)
ann_colors = list(cohort = c('discovery cohort' = "#F65C2F",
                             'validation cohort'= "#047FB9"))
class=c('WBC'="#FF3200",
        "LMR"= "#F65C2F", 
        'NLR'="#EE875F", 
        'PLR'="#E9AB82", 
        'FIB'="#E9C592", 
        'SII'="#E9DEA2", 
        'LDLC'= "#A9D5A8", 
        'HDLC'= "#5AC4AC",
        'CHOL'= "#18B1AF",
        'TRIG' = "#0E98B4",
        'LVEF' = "#047FB9",
        'LVMI' = "#0563A8")

#pdf("DMP_phenotype_heatmap2.pdf",height = 8,width = 20)
pheatmap(t(scale(coef)),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         #annotation_col = annotation_row,
         annotation_row = annotation_col,
         border_color = "NA",
         cellheight  = 10,
         annotation_colors = ann_colors,
         display_numbers = t(p),
         color = colorRampPalette(c("#ff7f00", "white", "lightseagreen"))(100)
         )
#dev.off()
