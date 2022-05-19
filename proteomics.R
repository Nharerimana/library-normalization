# title: QC and Normalization of Brain proteomics profiles
# author: Nadia Harerimana

library(plyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(GGally)
library(sva)
library(dplyr)
library(variancePartition)
library(boot)
library(corrplot)
library(lattice)
library(readr)
library(ggrepel)

abundance=read_delim('rushtmt_multiconsensus50batch_Proteins.txt',delim ='\t')
true_data=colnames(abundance)
abundace_name=c(grep('(Normalized)',true_data))
abundace_true=abundance[,abundace_name]
test=as.data.frame(colnames(abundace_true))
abundace_true=log2(abundace_true)

#### remove GIS, and creat GIS dataset for analysis ####
subject=colnames(abundace_true)
length(grep('126',subject))
length(grep('131',subject))
gis=c(grep('126',subject),grep('131',subject))
protein=abundace_true[,-gis]
GIS=abundace_true[,gis]

#### format proteomicsid
colnames(protein)=gsub(".*F","",colnames(protein))
colnames(protein)=gsub(",.*","",colnames(protein))
colnames(protein)=gsub(": ","_",colnames(protein))
colnames(protein)=paste('b',colnames(protein),sep='')
protein$protein=paste(abundance$`Gene Symbol`, abundance$Accession,sep='|')
protein=as.data.frame(protein)
rownames(protein)=protein$protein
protein=protein[,-401]
test=as.data.frame(colnames(protein))

#### format proteomicsid for GIS
colnames(GIS)=gsub(".*F","",colnames(GIS))
colnames(GIS)=gsub(",.*","",colnames(GIS))
colnames(GIS)=gsub(": ","_",colnames(GIS))
colnames(GIS)=paste('b',colnames(GIS),sep='')
GIS$protein=paste(abundance$`Gene Symbol`, abundance$Accession,sep='|')
GIS=as.data.frame(GIS)
rownames(GIS)=GIS$protein
GIS=GIS[,-101]
GIS_126=GIS[,c(1:50)]
GIS_131=GIS[,c(51:100)]


### drop proteins fall outside the 95% CI of regression lines
protein_new=protein
batch_proteinNA=NULL

test=apply(protein,1,function(x) sum(is.na(x)==T))
sum(test==400)


for(i in 1:50){
  gis_126=as.data.frame(GIS_126[,i])
  colnames(gis_126)=colnames(GIS_126)[i]
  rownames(gis_126)=rownames(GIS_126)
  gis_131=as.data.frame(GIS_131[,i])
  colnames(gis_131)=colnames(GIS_131)[i]
  rownames(gis_131)=rownames(GIS_131)
  data_i=cbind(gis_126,gis_131)
  colnames(data_i)=c('gis126','gis131')
  lm.model=lm(gis126~gis131,data=data_i)

  newx = as.data.frame(data_i$gis131)
  colnames(newx)='gis131'
  conf_interval <- as.data.frame(predict(lm.model, newdata=newx, interval="predict",
                           level = 0.95))
  data_i=cbind(data_i,conf_interval)
  data_i$outCI=ifelse(is.na(data_i$gis126)==T,NA,ifelse(data_i$gis126 > data_i$upr |
                                                          data_i$gis126 < data_i$lwr,1,0))
  protein_i=protein[,grep(paste('b',i,'_',sep=''),colnames(protein))]
  for(j in 1:dim(data_i)[1]){
    if( is.na(data_i$outCI[j])==F & data_i$outCI[j]==0)
      {protein_i[j,]=protein_i[j,]}else{protein_i[j,]=NA}
  }

  batch_proteinNA=c(batch_proteinNA,sum(data_i$outCI==0,na.rm=T))
  protein_new[,grep(paste('b',i,'_',sep=''),colnames(protein))]=protein_i

}
batch_proteinNA=batch_proteinNA/12691

protein_new$allNA=apply(protein_new,1,function(x) sum(is.na(x)==T))
sum(protein_new$allNA==400)

### get back to original N400 abundance and set NA in that with NA in protein_new
abundace_true=abundance[,abundace_name]
protein_3=abundace_true[,-gis]
colnames(protein_3)=gsub(".*F","",colnames(protein_3))
colnames(protein_3)=gsub(",.*","",colnames(protein_3))
colnames(protein_3)=gsub(": ","_",colnames(protein_3))
colnames(protein_3)=paste('b',colnames(protein_3),sep='')
## b1_127N, not b01_127N-------use ROSMAP_Traits.csv , not NEW_ROSMAP_Traits.csv !!!!!!!
protein_3$protein=paste(abundance$`Gene Symbol`, abundance$Accession,sep='|')
protein_3=as.data.frame(protein_3)
rownames(protein_3)=protein_3$protein
protein_3=protein_3[,-401]
test=as.data.frame(colnames(protein_3))
rownames(protein_3)

for(i in 1:dim(protein_3)[1]) {
  for(j in 1:dim(protein_3)[2]) {
    if(is.na(protein_new[i,j])==T) {
      protein_3[i,j] = NA
    } else {
      protein_3[i,j] = protein_3[i,j]
    }
  }
}

test=apply(protein_3,1,function(x) sum(is.na(x)==T))
sum(test==400)

### Save a copy of the Q/C'd data prior to imposing a missingness threshold ###

#write.csv(protein_3, "protein_3.csv")
# protein_3 <- read.csv("protein_3.csv")

#### use protein with more than 50% data ###
protein_3$no_NA=apply(protein_3,1,function(x)sum(is.na(x)==T))
# hist(protein_new$no_NA)

protein_3=protein_3[protein_3$no_NA < 0.5*(dim(protein_3)[2]-1),]
protein_3=protein_3[,-(dim(protein_3)[2])]

### Save a copy of the Q/C'd data after imposing 50% missingness threshold ###

#write.csv(protein_3, "protein_3.csv")

### For each subject: divided by sum ###
protein_3=sweep(protein_3,2,colSums(protein_3,na.rm = T),`/`)


### take log 2
protein_3=log2(protein_3)

protein_3=protein_3[,order(colnames(protein_3))]

#write.csv(protein_3,'protein_outlier_removed_N400_log2proportion.csv')

### phenotype: batch, MS, sex, age at death, PMI, study ###
pheno1=read_excel('/Users/nadiaharerimana/Box Sync/4Nadia/ROSMAP_Phenotype/Phe_082018/dataset_529_basic_08-12-2018.xlsx')
pheno1$projid=stri_replace_all_regex(pheno1$projid, "\\b0*(\\d+)\\b", "$1")
pheno1=pheno1[,c('projid','study','age_death','msex','cogdx')]
pheno2=read.csv('~/Box Sync/4Nadia/ROSMAP_Phenotype/Phe_082018/ROS_MAP_TRAITS.csv',header = T,stringsAsFactors = F)
pheno=merge(pheno2[,c('projid','proteomicsid','Batch','PMI','amyloid','tangles','Braak','CERAD')],pheno1,by='projid',all = T)
pheno=pheno[is.na(pheno$proteomicsid)==F,]
pheno=pheno[order(pheno$proteomicsid),]
pheno$diagnosis <- ifelse(pheno$cogdx == 1, 0,1)
table(pheno$diagnosis)

### Braak: 0= 0-3; 1= 4-6
pheno$Braak=ifelse(0<= pheno$Braak & pheno$Braak <= 3,0,ifelse(4 <= pheno$Braak & pheno$Braak <= 6,1,NA))
### Cerad: 0= 1-2; 1= 3-4
pheno$CERAD=ifelse(1 <= pheno$CERAD & pheno$CERAD <=2,0,ifelse(3 <= pheno$CERAD & pheno$CERAD <= 4,1,NA))
### study: ROS =1, MAP =0
pheno$study=ifelse(pheno$study=='ROS ',1,0)
### ms3: batch =1,2,3,5,11 -----ms=1; ms2: others -----ms=0
pheno$ms=ifelse(as.character(pheno$Batch) == '1'|as.character(pheno$Batch) =='2' |as.character(pheno$Batch)=='3'|as.character(pheno$Batch)=='5'|as.character(pheno$Batch)=='11',1,0)


###### regress out batch, MS, sex, age at death, PMI ######
res_log2abundance_after_regress=as.data.frame(matrix(NA,nrow=dim(protein_3)[1],ncol=dim(protein_3)[2]))
rownames(res_log2abundance_after_regress)=rownames(protein_3)
colnames(res_log2abundance_after_regress)=colnames(protein_3)

for(i in 1:dim(protein_3)[1]){
  test_protein=protein_3[i,]
  test_protein=as.data.frame(t(test_protein))
  colnames(test_protein)='log2abundance'
  test_protein$proteomicsid=rownames(test_protein)
  test_protein=merge(test_protein,pheno,by='proteomicsid',all.x = T)
  test_protein=test_protein[order(test_protein$proteomicsid),]

  lm_i=lm(log2abundance~as.factor(Batch)+PMI+age_death+msex+ms,data=test_protein,na.action=na.exclude)
  residual_i=residuals(lm_i)
  res_log2abundance_after_regress[i,]=residual_i
}

#### check ourliers by PCA ####
pca_data=as.data.frame(res_log2abundance_after_regress[complete.cases(res_log2abundance_after_regress),])
pcaObj = prcomp(t(pca_data), scale = T)
summary(pcaObj)
pcpoints = data.frame(pcaObj$x)[,1:2]
pcpoints$proteomicsid=rownames(pcpoints)
pcpoints=merge(pcpoints,pheno,all.x = T,by='proteomicsid')
pcpoints$pc1_outlier=ifelse(abs(pcpoints$PC1) > (mean(pcpoints$PC1)+4*sd(pcpoints$PC1)),pcpoints$proteomicsid,NA)
pcpoints$pc2_outlier=ifelse(abs(pcpoints$PC2) > (mean(pcpoints$PC2)+4*sd(pcpoints$PC2)),pcpoints$proteomicsid,NA)
pcpoints$outlier=ifelse(is.na(pcpoints$pc1_outlier)==F | is.na(pcpoints$pc2_outlier)==F, pcpoints$proteomicsid,NA)

# png("PCA_round1.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(Batch),label=pcpoints$outlier) +
  theme(legend.position="right")+ geom_text_repel(size=3)
# dev.off()
# drop #
people_to_drop=pcpoints$outlier[is.na(pcpoints$outlier)==F]

## re-run PCA w/o that subject
pca_data=as.data.frame(res_log2abundance_after_regress[complete.cases(res_log2abundance_after_regress),!colnames(res_log2abundance_after_regress) %in% people_to_drop])
pcaObj = prcomp(t(pca_data), scale = T)
summary(pcaObj)
pcpoints = data.frame(pcaObj$x)[,1:2]
pcpoints$proteomicsid=rownames(pcpoints)
pcpoints=merge(pcpoints,pheno,all.x = T,by='proteomicsid')
pcpoints$pc1_outlier=ifelse(abs(pcpoints$PC1) > (mean(pcpoints$PC1)+4*sd(pcpoints$PC1)),pcpoints$proteomicsid,NA)
pcpoints$pc2_outlier=ifelse(abs(pcpoints$PC2) > (mean(pcpoints$PC2)+4*sd(pcpoints$PC2)),pcpoints$proteomicsid,NA)
pcpoints$outlier=ifelse(is.na(pcpoints$pc1_outlier)==F | is.na(pcpoints$pc2_outlier)==F, pcpoints$proteomicsid,NA)

# png("PCA_round2.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(Batch),label=pcpoints$outlier) +
  theme(legend.position="right")+ geom_text_repel(size=3)
# dev.off()

# drop #
people_to_drop=c(people_to_drop,pcpoints$outlier[is.na(pcpoints$outlier)==F])

## re-run PCA
pca_data=as.data.frame(res_log2abundance_after_regress[complete.cases(res_log2abundance_after_regress),!colnames(res_log2abundance_after_regress) %in% people_to_drop])
pcaObj = prcomp(t(pca_data), scale = T)
summary(pcaObj)
pcpoints = data.frame(pcaObj$x)[,1:2]
pcpoints$proteomicsid=rownames(pcpoints)
pcpoints=merge(pcpoints,pheno,all.x = T,by='proteomicsid')
pcpoints$pc1_outlier=ifelse(abs(pcpoints$PC1) > (mean(pcpoints$PC1)+4*sd(pcpoints$PC1)),pcpoints$proteomicsid,NA)
pcpoints$pc2_outlier=ifelse(abs(pcpoints$PC2) > (mean(pcpoints$PC2)+4*sd(pcpoints$PC2)),pcpoints$proteomicsid,NA)
pcpoints$outlier=ifelse(is.na(pcpoints$pc1_outlier)==F | is.na(pcpoints$pc2_outlier)==F, pcpoints$proteomicsid,NA)

# png("PCA_round3.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(Batch),label=pcpoints$outlier) +
  theme(legend.position="right")+ geom_text_repel(size=3)
# dev.off()


## go back to log2(abundance), drop 9 samples, save it as a new raw log2 data ##
protein_3=protein_3[,which(colnames(protein_3) %in% colnames(pca_data))]
write.csv(protein_3,'protein_outlier_removed_N391_log2proportion.csv')

######## protein QC phase after removing bad samples: per protein basis  #####################

#### regress out effects of batch, MS2/MS3, sex, PMI, age_death, and study
new_res_log2abundance_after_regress=as.data.frame(matrix(NA,nrow=dim(protein_3)[1],ncol=dim(protein_3)[2]))
rownames(new_res_log2abundance_after_regress)=rownames(protein_3)
colnames(new_res_log2abundance_after_regress)=colnames(protein_3)
for(i in 1:dim(protein_3)[1]){
  test_protein=protein_3[i,]
  test_protein=as.data.frame(t(test_protein))
  colnames(test_protein)='log2abundance'
  test_protein$proteomicsid=rownames(test_protein)
  test_protein=merge(test_protein,pheno,by='proteomicsid',all.x = T)
  test_protein=test_protein[order(test_protein$proteomicsid),]

  lm_i=lm(log2abundance~as.factor(Batch)+PMI+age_death+msex+ms,data=test_protein,na.action=na.exclude)
  residual_i=residuals(lm_i)
  new_res_log2abundance_after_regress[i,]=residual_i
}
write.csv(new_res_log2abundance_after_regress,'n391_residual_log2_batchMSsexPMIage.csv')
#rm(new_res_log2abundance_after_regress)


############################################################################################

################# SVA adjustement to protect AD status #####################
protein <- read.csv("~/Box Sync/4Nadia/ROSMAP_TMT_proteomics/Wen_proteomics_QC/n391_residual_log2_batchMSsexPMIage.csv",stringsAsFactors = F, header = T)
rownames(protein) <- protein$X
protein <- protein[-1]
pheno <- pheno[match(colnames(protein), pheno$proteomicsid),]
colnames(protein) <- pheno$projid

# Now reate the full model matrix - including both the adjustment variables and the variable of interest (AD status)
mod = model.matrix(~as.factor(diagnosis), data=pheno)

# Now create the null model contains only the adjustment variables (batch)
mod0 = model.matrix(~1,data=pheno)

# Identifies the number of latent factors that need to be estimated
protein[is.na(protein)] = 0
protein <- as.matrix(protein)
n.sv = sva::num.sv(protein,mod,method="be", seed = 123456, B = 30)
n.sv # Identified 41 SVs
# Next apply the sva function to estimate the surrogate variables:

svobj = sva(protein,mod,mod0,n.sv=n.sv)$sv
svobj = data.frame(svobj)
colnames(svobj) = paste0('SV',1:dim(svobj)[2])
rownames(svobj) = pheno$projid
svobj10 <- svobj[,1:10]

# Use linear regression to residualize the first 10 SVs in the proteomics counts
protein_after_SVsregress=as.data.frame(matrix(NA,nrow=dim(protein)[1],ncol=dim(protein)[2]))
rownames(protein_after_SVsregress)=rownames(protein)
colnames(protein_after_SVsregress)=colnames(protein)
protein <- as.data.frame(protein)
svobj10$projid <- rownames(svobj10)
for(i in 1:dim(protein)[1]){
  protein_i=protein[i,]
  protein_i=as.data.frame(t(protein_i))
  colnames(protein_i)='proteinid'
  protein_i$projid=rownames(protein_i)
  data_i=merge(svobj10,protein_i,by='projid',all.x = T)
  lm_i=lm(proteinid~SV1+SV2+SV3+SV3+SV4+SV5+SV6+SV7+SV8+SV9+SV10,data= data_i)
  residual_i=residuals(lm_i)
  protein_after_SVsregress[i,]=residual_i
}
list.files()
write.csv(protein_after_SVsregress, file = "~/Box Sync/4Nadia/ROSMAP_proteomics/Wen_proteomics_QC/n392_residual_log2_batchMSsexPMIageSV10.csv")
#

}
############################################################################################
# Plot 1: Input: log2 ratio data (before regressing out sex, batch, PMI, age).
# Variables to include are: batch, MS2/MS3, sex, age_death, PMI, study, cogdx, sqrt(amyloid), sqrt(tangles)

#log2_protein=read.csv('protein_outlier_removed_N391_log2proportion.csv',header = T,stringsAsFactors = F)
#rownames(log2_protein)=log2_protein$X
#log2_protein=log2_protein[,-1]

#info=pheno
#rownames(info)=info$proteomicsid
#info=info[,-1]
#info=info[complete.cases(info),]
#info$Batch=as.factor(info$Batch)

#exprs=log2_protein
#exprs=exprs[,rownames(info)]
#exprs=exprs[complete.cases(exprs),]
#exprs=as.matrix(exprs)
# form = ~ age_death + sqrt_amyloid+ sqrt_tangles + PMI +msex +study + (1|Batch)+cogdx + ms
# form = ~ age_death + sqrt_amyloid+ sqrt_tangles + PMI +msex + (1|Batch)+cogdx + ms
#form = ~ age_death + sqrt_amyloid+ sqrt_tangles + PMI +msex + (1|Batch)+cogdx
#varPart <- fitExtractVarPartModel(exprs, form, info )
#vp <- sortCols( varPart )
#png("./plots/Var_before_regression_1.png", units="in", width=10, height=10, res=300)
#varPar_ggplot(vp)
#dev.off()

# Plot 2: Input: residual data (after regressing out sex, batch, PMI, age)
new_res_log2abundance_after_regress=read.csv('n391_residual_log2_batchMSsexPMIageStudy.csv',header = T,stringsAsFactors = F)
rownames(new_res_log2abundance_after_regress)=new_res_log2abundance_after_regress$X
new_res_log2abundance_after_regress=new_res_log2abundance_after_regress[,-1]
info=pheno
rownames(info)=info$projid
info=info[,-1]
info=info[complete.cases(info),]
info$Batch=as.factor(info$Batch)
info$study <- ifelse(info$study=="ROS",0, 1)
exprs=protein
exprs=exprs[,rownames(info)]
exprs[is.na(exprs)] <- 0
exprs=as.matrix(exprs)


# form = ~ age_death + sqrt_amyloid+ sqrt_tangles + PMI +msex +study + (1|Batch)+cogdx + ms
# form = ~ age_death + sqrt_amyloid+ sqrt_tangles + PMI +msex + (1|Batch)+cogdx + ms
form = ~ age_death + sqrt_amyloid+ sqrt_tangles +PMI +msex + (1|Batch)+ cogdx +study
varPart <- fitExtractVarPartModel(exprs, form, info )
vp <- sortCols( varPart )
png("./plots/Var_after_regression_1.png", units="in", width=10, height=10, res=300)
varPar_ggplot(vp)
dev.off()

# Plot 3: Input: residual data (after regressing out sex, batch, PMI, age), add hidden factor 1
form = ~ age_death + sqrt_amyloid+ sqrt_tangles + PMI +msex + (1|Batch)+cogdx+hiddenFactors.1
varPart <- fitExtractVarPartModel(exprs, form, info )
vp <- sortCols( varPart )
png("./plots/Var_after_regression_hiddenfactor_1.png", units="in", width=10, height=10, res=300)
varPar_ggplot(vp)
dev.off()

# Plot 4: Input: residual data (after regressing out sex, batch, PMI, age), add hidden factor 1 2 3
form = ~ age_death + sqrt_amyloid+ sqrt_tangles + PMI +msex + (1|Batch)+cogdx+hiddenFactors.1+hiddenFactors.2+hiddenFactors.3
varPart <- fitExtractVarPartModel(exprs, form, info )
vp <- sortCols( varPart )
png("./plots/Var_after_regression_hiddenfactor_123.png", units="in", width=10, height=10, res=300)
varPar_ggplot(vp)
dev.off()


###################################### PC plots ##################################################

## plot1. before regression
## use pheno for N=391!!!!!!!!!
pheno=pheno[pheno$proteomicsid %in% colnames(log2_protein),]
pca_data=as.data.frame(log2_protein[complete.cases(log2_protein),])
pcaObj = prcomp(t(pca_data), scale = T)
pcpoints = data.frame(pcaObj$x)[,1:2]
pcpoints$proteomicsid=rownames(pcpoints)
pcpoints=merge(pcpoints,pheno,all.x = T,by='proteomicsid')

png("./plots/PC_before_regression_batch.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(Batch)) +theme(legend.position="right")
dev.off()

png("./plots/PC_before_regression_MS.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(ms)) +theme(legend.position="right")
dev.off()

png("./plots/PC_before_regression_cogdx.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(cogdx)) +theme(legend.position="right")
dev.off()

png("./plots/PC_before_regression_sex.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(msex)) +theme(legend.position="right")
dev.off()


## plot2. after regression
pca_data=as.data.frame(new_res_log2abundance_after_regress[complete.cases(new_res_log2abundance_after_regress),])
pcaObj = prcomp(t(pca_data), scale = T)
pcpoints = data.frame(pcaObj$x)[,1:2]
pcpoints$proteomicsid=rownames(pcpoints)
pcpoints=merge(pcpoints,pheno,all.x = T,by='proteomicsid')

png("./plots/PC_after_regression_batch.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(Batch)) +theme(legend.position="right")
dev.off()

png("./plots/PC_after_regression_MS.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(ms)) +theme(legend.position="right")
dev.off()

png("./plots/PC_after_regression_cogdx.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(cogdx)) +theme(legend.position="right")
dev.off()

png("./plots/PC_after_regression_sex.png", units="in", width=10, height=10, res=300)
qplot(x=PC1, y=PC2, data=pcpoints, colour=as.factor(msex)) +theme(legend.position="right")
dev.off()



################################################################################################

################ correlation analysis ##################

# use pheno N=365!!!
shapiro.test(pheno$hiddenFactors.1)
shapiro.test(pheno$hiddenFactors.2)
shapiro.test(pheno$hiddenFactors.3)
shapiro.test(pheno$hiddenFactors.4)
shapiro.test(pheno$hiddenFactors.5)
shapiro.test(pheno$hiddenFactors.6)
shapiro.test(pheno$hiddenFactors.7)
shapiro.test(pheno$hiddenFactors.8)
shapiro.test(pheno$hiddenFactors.9)
shapiro.test(pheno$hiddenFactors.10)

pheno_for_cor=pheno[,c(2,4,10,13,15:27)]

correlation=as.data.frame(matrix(NA,nrow=6,ncol=10))
colnames(correlation)=colnames(pheno_for_cor)[8:17]
rownames(correlation)=colnames(pheno_for_cor)[2:7]
for(i in 8:17){
  for(j in 2:7){
  data_i=pheno_for_cor[,c(i,j)]
  correlation_i=cor(data_i[complete.cases(data_i),],method='spearman')[1,2]
  correlation[j-1,i-7]=correlation_i}
}

write.csv(correlation,'hiddenfactors_correlation.csv')

pvalue=NULL
for(i in 18:27){
fit=lm(pheno[,i]~as.factor(pheno$Batch))
pvalue_i=anova(fit)$`Pr(>F)`[1]
pvalue=c(pvalue,pvalue_i)
}

##########################################################################################

################ generate sharma cell type that are profiled in ROSMAP ###################

rm(list=ls())
load('/Users/twingo/Box Sync/4Wen/Enrichment_analysis/proteomics_cogdec/sharma_celltype.RData')

protein=read.csv('n391_residual_log2_batchMSsexPMIageStudy.csv')
colnames(protein)[1]='protein'
protein$protein=sub("\\|.*", "", protein$protein)
protein_name=as.data.frame(unique(protein$protein))

sharma_list[["Isolated Astrocytes"]]=sharma_list[["Isolated Astrocytes"]][sharma_list[["Isolated Astrocytes"]] %in% protein_name$`unique(protein$protein)`]
sharma_list[["Isolated Microglia"]]=sharma_list[["Isolated Microglia"]][sharma_list[["Isolated Microglia"]] %in% protein_name$`unique(protein$protein)`]
sharma_list[["Isolated Neurons"]]=sharma_list[["Isolated Neurons"]][sharma_list[["Isolated Neurons"]] %in% protein_name$`unique(protein$protein)`]
sharma_list[["Isolated Oligodendrocytes"]]=sharma_list[["Isolated Oligodendrocytes"]][sharma_list[["Isolated Oligodendrocytes"]] %in% protein_name$`unique(protein$protein)`]

save(sharma_list,file='sharma_celltype_inROSMAP.RData')
