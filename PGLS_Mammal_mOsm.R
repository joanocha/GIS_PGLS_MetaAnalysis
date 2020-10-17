#useful packages
library(ggpubr)
library(Hmisc)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(plotly)
library(psych)
library(ape)
library(nlme)
library(geiger)
library(phytools)
library(phylobase)
library(plyr) 
library(reshape2)

## DATASETS AND SUBSETS
metaUrine <- read.delim('METAdata.txt', na.strings=c("", "NA"), header=TRUE, row.names=1)
dehydrated <-metaUrine[which(metaUrine$control=='D'),]
bmr<-metaUrine[which(metaUrine$BMR.m_O2.g.h!='NA'),]
arvore<-read.tree("RAxML_4098mam1out.newick")
prun_list<-name.check(arvore, metaUrine, data.names= metaUrine$tip.label)
dehydrated_list<-name.check(arvore, dehydrated, data.names= dehydrated$tip.label)
bmr_list<-name.check(arvore, bmr, data.names= bmr$tip.label)
arvore_cortada<-drop.tip(arvore, prun_list$tree_not_data)
metaUrine <-metaUrine[arvore_cortada$tip.label,] #ensure meta urine and arvore cortada are order the same way
arvore_d<-drop.tip(arvore, dehydrated_list$tree_not_data)
dehydrated<-dehydrated[arvore_d$tip.label,] #ensure meta urine and arvore cortada are order the same way
arvore_bmr<-drop.tip(arvore, bmr_list$tree_not_data)
bmr<-bmr[arvore_bmr$tip.label,] #ensure meta urine and arvore cortada are order the same way


#NORMALITY TEST
osm<-ggqqplot(metaUrine$Max.mosm.Kg,  title = "Maximum urine osmolality (mOsm/Kg)") 
logosm<-ggqqplot(log10(metaUrine$Max.mosm.Kg),  title = "Log10(Maximum urine osmolality)") 
mass<-ggqqplot(metaUrine$body.mass.Kg, title = "Body mass (Kg)") 
logmass<-ggqqplot(log10(metaUrine$body.mass.Kg),  title = "Log10(Body mass)") 
BMR<-ggqqplot(metaUrine$BMR.m_O2.g.h,  title = "Basal Metabolic Rate (ml O2 per g.h)") 
logBMR<-ggqqplot(log10(metaUrine$BMR.m_O2.g.h),  title = "Log10(Basal Metabolic Rate)") 
qqplot<-ggarrange(osm,logosm, BMR, logBMR, mass, logmass, legend = "top", labels = c("A", "B", "C", "D", "E", "F"), hjust=0, vjust = 1, nrow = 3, ncol = 2, align = "h" )
ggsave("qqplot.pdf", qqplot, "pdf", dpi=300, width=12, height=7, units="in")

## LOG
metaUrine$Max.mosm.Kg.log <- log10(metaUrine$Max.mosm.Kg)
metaUrine$body.mass.Kg.log<- log10(metaUrine$body.mass.Kg)
metaUrine$BMR.m_O2.g.h.log<-log10(metaUrine$BMR.m_O2.g.h)
dehydrated$Max.mosm.Kg.log <- log10(dehydrated$Max.mosm.Kg)
dehydrated$body.mass.Kg.log<- log10(dehydrated$body.mass.Kg)
bmr$Max.mosm.Kg.log <- log10(bmr$Max.mosm.Kg)
bmr$body.mass.Kg.log<- log10(bmr$body.mass.Kg)
bmr$body.mass.g.log<- log10(bmr$body.mass.g)
bmr$BMR.m_O2.g.h.log<-log10(bmr$BMR.m_O2.g.h)

# CORRELATION
# mOSM/AI/BM
log.bodymass.mOsm<-ggscatter(metaUrine, x = "body.mass.Kg.log", y = "Max.mosm.Kg.log",
                             add = "reg.line",                                
                             conf.int = TRUE,                                  
                             add.params = list(color = "blue", fill = "lightgray"))+ stat_cor(method = "spearman", label.x = 0, label.y = 4.0) + # Add correlation coefficient 
xlab("Log 10 Body mass (Kg)") + ylab("Log 10 Maximum urine osmolality (mOsm/Kg)")
                        
AI.logmOsm<-ggscatter(metaUrine, x = "AI", y = "Max.mosm.Kg.log", 
                      add = "reg.line",                                
                      conf.int = TRUE,                                  
                      add.params = list(color = "blue", fill = "lightgray"))+ stat_cor(method = "spearman", label.x = 0.5, label.y = 4.0) + # Add correlation coefficient 
  xlab("Mean Aridity Index") + ylab("Log 10 Maximum urine osmolality (mOsm/Kg)")

AI.log.bodymass<-ggscatter(metaUrine, x = "AI", y = "body.mass.Kg.log", 
                           add = "reg.line",                                
                           conf.int = TRUE,                                  
                           add.params = list(color = "blue", fill = "lightgray"))+ stat_cor(method = "spearman", label.x = 0.5, label.y = 4.0) + # Add correlation coefficient
  xlab("Mean Aridity Index") + ylab("Log 10 Body mass (Kg)")

urine_ai<-ggarrange(log.bodymass.mOsm,AI.logmOsm,AI.log.bodymass, labels = c("A", "B", "C"), ncol=3, nrow=1)
ggsave("urine_ai.tiff", urine_ai, "tiff", dpi=900, width=12, height=5, units="in")


###  BMR/BIOS
BIO1.MassAdjBMR<-ggscatter(metaUrine, x = "mean_bio1", y = "BMR.m_O2.g.h.log", 
                           add = "reg.line", conf.int = TRUE, 
                           cor.coef = TRUE, cor.method = "spearman",label.x = -8, label.y = -8.0,
                           xlab = "Mean Annual Average Temperature", ylab = "Log mass adjusted BMR (ml O2 per g per h)")
BIO5.MassAdjBMR<-ggscatter(metaUrine, x = "mean_bio5", y = "BMR.m_O2.g.h.log", 
                           add = "reg.line", conf.int = TRUE, 
                           cor.coef = TRUE, cor.method = "spearman",label.x = -8, label.y = -8.0,
                           xlab = "Mean Maximum Temperture warmest month", ylab = "Log mass adjusted BMR (ml O2 per g per h)")
BIO6.MassAdjBMR<-ggscatter(metaUrine, x = "mean_bio6", y = "BMR.m_O2.g.h.log", 
                           add = "reg.line", conf.int = TRUE, 
                           cor.coef = TRUE, cor.method = "spearman",label.x = -8, label.y = -8.0,
                           xlab = "Mean Minimum Temperture coldest month", ylab = "Log mass adjusted BMR (ml O2 per g per h)")
BIO7.MassAdjBMR<-ggscatter(metaUrine, x = "mean_bio7", y = "BMR.m_O2.g.h.log", 
                           add = "reg.line", conf.int = TRUE, 
                           cor.coef = TRUE, cor.method = "spearman",label.x = -8, label.y = -8.0,
                           xlab = "Mean Temperature Annual Range", ylab = "Log mass adjusted BMR (ml O2 per g per h)")
bio_bmr<-ggarrange(BIO1.MassAdjBMR,BIO5.MassAdjBMR,BIO6.MassAdjBMR,BIO7.MassAdjBMR, labels = c("A", "B", "C", "D"), ncol=2, nrow=2)
ggsave("bios_bmr.tiff", bio_bmr, "tiff", dpi=900, width=12, height=12, units="in")

### mOSm/BIOs
BIO1.logmOsm<-ggscatter(metaUrine, x = "mean_bio1", y = "Max.mosm.Kg.log", 
                           add = "reg.line", conf.int = TRUE, 
                           cor.coef = TRUE, cor.method = "spearman",
                           xlab = "Mean Annual Average Temperature", ylab = "Log 10 Maximum urine osmolality (mOsm/Kg)")
BIO5.logmOsm<-ggscatter(metaUrine, x = "mean_bio5", y = "Max.mosm.Kg.log", 
                           add = "reg.line", conf.int = TRUE, 
                           cor.coef = TRUE, cor.method = "spearman",
                           xlab = "Mean Maximum Temperture warmest month", ylab = "Log 10 Maximum urine osmolality (mOsm/Kg)")
BIO6.logmOsm<-ggscatter(metaUrine, x = "mean_bio6", y = "Max.mosm.Kg.log", 
                           add = "reg.line", conf.int = TRUE, 
                           cor.coef = TRUE, cor.method = "spearman",
                           xlab = "Mean Minimum Temperture coldest month", ylab = "Log 10 Maximum urine osmolality (mOsm/Kg)")
BIO7.logmOsm<-ggscatter(metaUrine, x = "mean_bio7", y = "Max.mosm.Kg.log", 
                           add = "reg.line", conf.int = TRUE, 
                           cor.coef = TRUE, cor.method = "spearman",
                           xlab = "Mean Temperature Annual Range", ylab = "Log 10 Maximum urine osmolality (mOsm/Kg)")
bios_mOsm<-ggarrange(BIO1.logmOsm, BIO5.logmOsm, BIO6.logmOsm, BIO7.logmOsm, labels = c("A", "B", "C", "D"), ncol=2, nrow=2)
ggsave("bios_mOsms.tiff", bios_mOsm, "tiff", dpi=900, width=10, height=10, units="in")

### AI/BIOs
BIO1.AI<-ggscatter(metaUrine, x = "mean_bio1", y = "AI", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "spearman",label.x = -8, label.y = -8.0,
                        xlab = "Mean Annual Average Temperature", ylab = "Mean Aridity Index")
BIO5.AI<-ggscatter(metaUrine, x = "mean_bio5", y = "AI", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "spearman",label.x = -8, label.y = -8.0,
                        xlab = "Mean Maximum Temperture warmest month", ylab = "Mean Aridity Index")
BIO6.AI<-ggscatter(metaUrine, x = "mean_bio6", y = "AI", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "spearman",label.x = -8, label.y = -8.0,
                        xlab = "Mean Minimum Temperture coldest month", ylab = "Mean Aridity Index")
BIO7.AI<-ggscatter(metaUrine, x = "mean_bio7", y = "AI", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "spearman",label.x = -8, label.y = -8.0,
                        xlab = "Mean Temperature Annual Range", ylab = "Mean Aridity Index")
bios_AI<-ggarrange(BIO1.AI, BIO5.AI, BIO6.AI, BIO7.AI, labels = c("A", "B", "C", "D"), ncol=2, nrow=2)
ggsave("AI.tiff",bios_AI, "tiff", dpi=900, width=12, height=12, units="in")


#PGLS - example model

# aLL 
model<-gls(Max.mosm.Kg.log~AI+body.mass.Kg.log, data=metaUrine, correlation=corPagel(1,arvore_cortada))
summary(model)
## dataset with BMR
model<-gls(BMR.m_O2.g.h.log~ body.mass.g.log, data=bmr, correlation=corPagel(1,arvore_bmr))
summary(model)
# dataset with dehydrated only
model<-gls(Max.mosm.Kg.log~AI+body.mass.Kg.log, data=dehydrated, correlation=corPagel(1,arvore_d))
summary(model)

## Phylogenetic signal for variables
#example for AI
phylosig(arvore_cortada, metaUrine$AI, method="lambda", test=TRUE) # do it for each variable

#lambda tells the relative importance of phylogeny in predicting the trait


# ANCESTRAL STATE RECONSTRUCTION
arvore_cortada$tip.label<-metaUrine$species
arvore_cortada$Classification<-metaUrine$Classification
arvore_cortada$AI<-metaUrine$AI
plotTree(arvore_cortada,ftype="i",fsize=0.6,lwd=1)
X2<-metaUrine[!duplicated(metaUrine[,'species']),]
row.names(X2)<-unique(X2[,'species'])
AI<-as.matrix(X2)[,14] #choosing the 3rd column of dataset X2
mode(AI)<-'numeric' #should be a named numeric matrix
str(X2) #continuous variables should be labeled as 'numeric'
X2$AI<-as.numeric(X2$AI)
AI<-AI[arvore_cortada$tip.label]
mOsm<-as.matrix(X2)[,12] #choosing the 3rd column of dataset X2
mode(mOsm)<-'numeric' #should be a named numeric matrix
str(X2) #continuous variables should be labeled as 'numeric'
X2$mOsm<-as.numeric(X2$mOsm)
mOsm<-mOsm[arvore_cortada$tip.label]
obj<-contMap(arvore_cortada,AI,plot=FALSE)
osm<-contMap(arvore_cortada,mOsm,plot=FALSE)
grey<-setMap(osm,c( "#FFFFD4","#FED98E", "#FE9929", "#D95F0E", "#993404"))
colormap(colormap = colormaps$viridis, nshades = 6, format = "hex",
         alpha = 1, reverse = FALSE)
arid<-setMap(obj,c("#a50026", "#d73027", "#fee090", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"))
pdf(file = "AncRec_tree.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 13) # The height of the plot in inches
par(mfrow=c(1,2))
plot(arid,lwd=c(3,7),outline=FALSE,xlim=c(-0.1,3), fsize=c(0.7,1),leg.txt="AI")
plot(grey,lwd=c(3,7),outline=FALSE, direction="leftwards",ftype="off",xlim=c(0.1,3), leg.txt="mOsm")
dev.off()


# GRAPHICAL ABSTRACT Mammal Review
plotmetaUrine <- arrange(transform(metaUrine,Classification=factor(Classification,levels=c("Arid", "Semi-Arid", "Dry sub-humid", "Humid"))),Classification)
osm_AI<-ggplot(plotmetaUrine, aes(x=AI, y=Max.mosm.Kg, color = Classification, label=Name)) + xlab("Mean annual aridity index (AI)") + ylab("Maximum urine osmolality (mOsm/Kg)") + 
  geom_point(size=1) + 
  geom_text(aes(label=ifelse(
    Max.mosm.Kg ==9374 |
      Max.mosm.Kg ==1880|
      Max.mosm.Kg ==4650| 
      Max.mosm.Kg ==2362|
      Max.mosm.Kg ==1390|
      Max.mosm.Kg ==4470| 
      Max.mosm.Kg ==7767| 
      Max.mosm.Kg ==6500|
      Max.mosm.Kg ==3170|
      Max.mosm.Kg ==2608|
      Max.mosm.Kg ==9370|
      Max.mosm.Kg ==8773|
      Max.mosm.Kg ==4022 | 
      Max.mosm.Kg ==3027 | 
      Max.mosm.Kg ==837| 
      Max.mosm.Kg ==537,
    as.character(species),'')),hjust=0, vjust=0, size=2.9) +
  scale_color_brewer(palette ="RdYlBu", direction=1) +
  theme_cowplot()
plotmetaUrine <- arrange(transform(metaUrine,Classification=factor(Classification,levels=c("Arid", "Semi-Arid", "Dry sub-humid", "Humid"))),Classification)


### PLOT for main figure in Mammal Review paper
plotmetaUrine <- arrange(transform(metaUrine,Classification=factor(Classification,levels=c("Arid", "Semi-Arid", "Dry sub-humid", "Humid"))),Classification)
sp <- ggscatter(plotmetaUrine, x = "AI", y = "Max.mosm.Kg",
                color = "Classification", palette = "jco",
                size = 1.5, alpha = 0.6, ggtheme = theme_cowplot()) +
  scale_color_brewer(palette ="RdYlBu", direction=1) + 
  xlab("Mean Aridity index") + ylab("Maximum urine osmolality (mOsm/Kg)")            
# Marginal boxplot of x (top panel) and y (right panel)
metaUrine <- arrange(transform(metaUrine,Classification=factor(Classification,levels=c("Arid", "Semi-Arid", "Dry sub-humid", "Humid"))),Classification)
xplot <- ggboxplot(plotmetaUrine, x = "Category", y = "Max.mosm.Kg", 
                   color = "Classification", fill = "Classification", palette = "RdYlBu",
                   alpha = 0.5, ggtheme = theme_cowplot()) + xlab("") + ylab("Maximum urine osmolality (mOsm/Kg)") 
library(cowplot)
myplot<-plot_grid(xplot, sp, ncol = 1,labels = c("A", "B"))
ggsave("myplot.tiff", myplot, "tiff", dpi=600, width=10, height=9, units="in")



### BOX PLOT for Trends in Ecology and Evolution
#all 
plotmetaUrine <- arrange(transform(metaUrine,Classification=factor(Classification,levels=c("Arid", "Semi-Arid", "Dry sub-humid", "Humid"))),Classification)
ggplot(plotmetaUrine, aes(x=Category, y=Max.mosm.Kg, fill=Classification)) + xlab("") + ylab("Maximum urine osmolality (mOsm/Kg)") + 
  geom_boxplot() + 
  scale_x_discrete(limits=c("Evader", "Evaporator", "Endurer")) +
  scale_fill_manual(values=c("#FC4E07", "#C4961A", "#C3D7A4", "#00AFBB","#4E84C4")) +
  theme_classic2()
# dehydrated only
plotmetadehydrated <- arrange(transform(dehydrated,Classification=factor(Classification,levels=c("Arid", "Semi-Arid", "Dry sub-humid", "Humid"))),Classification)
treeReview<-ggplot(plotmetadehydrated, aes(x=Category, y=Max.mosm.Kg, fill=Classification)) + xlab("") + ylab("Maximum urine osmolality (mOsm/Kg)") + 
  geom_boxplot() + 
  scale_x_discrete(limits=c("Evader", "Evaporator", "Endurer")) +
  scale_fill_manual(values=c("#FC4E07", "#C4961A", "#C3D7A4", "#00AFBB","#4E84C4")) +
  theme_classic2()
ggsave("treeReview.tiff", treeReview, "tiff", dpi=900, width=5, height=4, units="in")


