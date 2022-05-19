#title: Covariate and DE analysis of reprocessed counts method adapted from Thanner Perumal
#author: Nadia Harerimana
# Date: March 9th 2020
# Load libraries
devtools::install_github('th1vairam/CovariateAnalysis@dev')
library(CovariateAnalysis) 
library(data.table)
library(tidyverse)
library(dplyr)
library(psych)
library(limma)
library(edgeR)
library(biomaRt)
library(RColorBrewer)
library(preprocessCore)
library(BiocManager)
library(cqn)
library(variancePartition)
library(githubr) # get the package from devtools::install_github('brian-bot/githubr')
library(knitr)
library(iterators)
library(tibble)
library(parallel)
library(doParallel)
library(foreach)
library(ComplexHeatmap)
library(grid)
library(vcd)
library(stringi)


# Load count matrix 
COUNT <- read.table("ROSMAP_Regressed_Counts.tsv")
sampleid <- COUNT[1,]
sampleid <- as.character(sampleid)
colnames(COUNT) <- sampleid
COUNT <- COUNT[-1,]
ensembl_gene_id <- COUNT[,1]
rownames(COUNT) <- ensembl_gene_id
COUNT <- COUNT[,-1]
sampleid <- colnames(COUNT)
ensembl_gene_id <- rownames(COUNT)

# Load technical covariates
metadata <- read.table("ROSMAP_SVA_Covariates.tsv", header = T)

### Divide metadata by brain location 
metadata_ACC <- metadata %>%
      dplyr::filter(tissue == "ACC")
length(unique(metadata_ACC$individualID)) == nrow(metadata_ACC)
n_occur <- data.frame(table(metadata_ACC$individualID))
n_occur[n_occur$Freq > 1,]
ACC_ID <- metadata_ACC$SampleID
Samples_to_drop_ACC=c(grep('Sample_SM-AYOSR',ACC_ID),grep('Sample_SM-AYOTP',ACC_ID))
metadata_ACC <-metadata_ACC[-Samples_to_drop_ACC,]
metadata_ACC$sampleIdentifier <- paste("Sample", metadata_ACC$individualID, metadata_ACC$tissue, sep = "_") 
metadata_ACC <- metadata_ACC[,c("sampleIdentifier", "SampleID","individualID","Diagnosis","race",
                                "spanish", "cogdx", "apoe4","sex","batch", "notes", "pmi" , "rin",
                                "rin2", "age_death" ,"pct_pf_reads_aligned", "pct_intronic_bases" ,
                                "pct_intergenic_bases", "pct_coding_bases", "tissue", "Tissue.Diagnosis",
                                "Tissue.APOE4", "Tissue.Diagnosis.SEX" )]


metadata_DLPC<- metadata %>%
  dplyr::filter(tissue == "DLPFC")
length(unique(metadata_DLPC$individualID)) == nrow(metadata_DLPC)
n_occur <- data.frame(table(metadata_DLPC$individualID))
n_occur[n_occur$Freq > 1,]
DLPC_ID <- metadata_DLPC$SampleID
Samples_to_drop_DLPC=c(grep('Sample_SM-2T717',DLPC_ID),grep('Sample_SM-AYOSS',DLPC_ID),
                      grep('Sample_SM-AYO4G',DLPC_ID),grep('RISK_217',DLPC_ID),
                      grep('RISK_218',DLPC_ID),grep('Sample_SM-2T6Z5',DLPC_ID),
                      grep('Sample_R6284240-PCC',DLPC_ID),grep('RISK_15',DLPC_ID))
metadata_DLPC <-metadata_DLPC[-Samples_to_drop_DLPC,]
metadata_DLPC$sampleIdentifier <- paste("Sample", metadata_DLPC$individualID, metadata_DLPC$tissue, sep = "_") 
metadata_DLPC<- metadata_DLPC[,c("sampleIdentifier", "SampleID","individualID","Diagnosis","race",
                                "spanish", "cogdx", "apoe4","sex","batch", "notes", "pmi" , "rin",
                                "rin2", "age_death" ,"pct_pf_reads_aligned", "pct_intronic_bases" ,
                                "pct_intergenic_bases", "pct_coding_bases", "tissue", "Tissue.Diagnosis",
                                "Tissue.APOE4", "Tissue.Diagnosis.SEX" )]


metadata_PCC<- metadata %>%
  dplyr::filter(tissue == "PCC")
length(unique(metadata_PCC$individualID)) == nrow(metadata_PCC)
metadata_PCC$sampleIdentifier <- paste("Sample", metadata_PCC$individualID, metadata_PCC$tissue, sep = "_") 
metadata_PCC<- metadata_PCC[,c("sampleIdentifier", "SampleID","individualID","Diagnosis","race",
                                 "spanish", "cogdx", "apoe4","sex","batch", "notes", "pmi" , "rin",
                                 "rin2", "age_death" ,"pct_pf_reads_aligned", "pct_intronic_bases" ,
                                 "pct_intergenic_bases", "pct_coding_bases", "tissue", "Tissue.Diagnosis",
                                 "Tissue.APOE4", "Tissue.Diagnosis.SEX" )]

metadata <- rbind(metadata_ACC, metadata_DLPC, metadata_PCC)  

# Pick higher quality RIN batch for sample 492_120515
metadata <- metadata %>%
  dplyr::group_by(sampleIdentifier) %>%
  dplyr::top_n(1, rin)

### Data preprocessing
# Remove samples with no cogdx, RIN, PMI scores and age_death
    metadata<- metadata %>%
      ungroup %>%
        dplyr::filter(!is.na(age_death), 
                  !is.na(cogdx), 
                  !is.na(pmi), 
                  !is.na(rin),
                  !is.na(sampleIdentifier),
                  !is.na(Tissue.Diagnosis),
                  rin >= 5)

 # Fix missing batch, if any
  levels(metadata$batch)[levels(metadata$batch) == ''] = 'NoBatch'       

  # Match metadata to expression data
  indToRetain = intersect(metadata$SampleID, colnames(COUNT))
  removedIDs = setdiff(colnames(COUNT), metadata$SampleID)
  COUNT = COUNT[,indToRetain]
  rownames(metadata) = metadata$sampleIdentifier
  colnames(COUNT) <- metadata$sampleIdentifier
  COUNT <- as.matrix(COUNT)
  class(COUNT) <- "numeric"
  
  
  ### Covariate clustering
  #Determine relationship between metadata
  primaryVariable <- c("cogdx", "Tissue.Diagnosis", "Tissue.APOE4")
  FactorCovariates <- c( "batch", "sex", "race", "Tissue.Diagnosis", "Tissue.APOE4", "cogdx")
  ContCovariates <- c("rin","age_death", "pmi","pct_coding_bases","pct_intergenic_bases", "pct_pf_reads_aligned","pct_intronic_bases")

  # Find inter relation between factor covariates
  covariates = metadata[,c(FactorCovariates,ContCovariates),drop=F]
  covariates <- data.frame(lapply( covariates,function(x){
    x <- sapply(x,function(y){str_replace_all(as.character(y),'\\+','')})}))
  rownames(covariates) <- metadata$sampleIdentifier
 
  # Convert factor covariates to factors
  covariates[,FactorCovariates] = lapply(covariates[,FactorCovariates], factor)
  covariates[,ContCovariates] = lapply(covariates[,ContCovariates], as.character)
  covariates[,ContCovariates] = lapply(covariates[,ContCovariates], as.numeric)
  
  #Correlation/association between covariates at an FDR <= 0.1
  png("~/Emory University/WingoLab - Files/Nadia Harerimana/ROSMAP RNAseq/Correlationassociation_btn_covariates_FDR<=0.1.png", units="in", width=10, height=10, res=300)
  COVARIATES.CORRELATION = CovariateAnalysis::getAssociationStatistics(covariates, PVAL = 0.05)
  draw(COVARIATES.CORRELATION$plot, heatmap_legend_side = 'left', padding  = unit(c(18,2,2,18), 'mm'))
  dev.off()   
  
  
  ### Explore metadata
  my.theme = theme(legend.position = 'top', axis.text.x = element_text(angle = 90, hjust = 1))
  # RIN
  p = list()
  p[[1]] = ggplot(covariates, aes(x = Tissue.Diagnosis, y = rin)) + geom_boxplot()
  p[[1]] = p[[1]] + ggtitle('RIN') +  my.theme
  # AgeAtDeath
  p[[2]] = ggplot(covariates, aes(x = Tissue.Diagnosis, y = age_death)) + geom_boxplot()
  p[[2]] = p[[2]] + ggtitle('AgeOfDeath') +  my.theme
  # PMI
  p[[3]] = ggplot(covariates, aes(x = Tissue.Diagnosis, y = pmi)) + geom_boxplot()
  p[[3]] = p[[3]] + ggtitle('PMI') +  my.theme
  # Intergenic bases
  p[[4]] = ggplot(covariates, aes(x = Tissue.Diagnosis, y = pct_intergenic_bases)) + geom_boxplot()
  p[[4]] = p[[4]] + ggtitle('Fraction Intergenic Bases') +  my.theme
  # Intronic bases
  p[[5]] = ggplot(covariates, aes(x = Tissue.Diagnosis, y = pct_intronic_bases)) + geom_boxplot()
  p[[5]] = p[[5]] + ggtitle('Fraction Intronic Bases') +  my.theme
  # Ribosomal bases
  p[[6]] = ggplot(covariates, aes(x = Tissue.Diagnosis, y = pct_coding_bases)) + geom_boxplot()
  p[[6]] = p[[6]] + ggtitle('Fraction Coding Bases') +  my.theme
  
  png("~/Diagnosis_vs_continouscovariates.png", units="in", width=10, height=6, res=300)
  multiplot(plotlist = p, cols = 3)
  dev.off()
  

  # Get gene specific parameters from synapse
  GENE.PARAM <- read.delim("~/geneParameters.tsv")
  GENE.LEN = dplyr::select(GENE.PARAM, ensembl_gene_id, gene.length) %>% 
    unique() 
  rownames(GENE.LEN) = GENE.LEN$ensembl_gene_id
  GENE.GC = dplyr::select(GENE.PARAM, ensembl_gene_id, percentage_gc_content) %>% 
    unique() 
  rownames(GENE.GC) = GENE.GC$ensembl_gene_id 
  
 ### Filter genes
  #* Remove genes that have less than 1 cpm counts in at least 50% of samples per Tissue x Diagnosis
  # Remove genes with missing gene length and percentage GC content

  genesToAnalyze = covariates %>%
    rownameToFirstColumn('SampleID') %>%
    dlply(.(Tissue.Diagnosis), .fun = function(mtd, count){
      processed.counts = getGeneFilteredGeneExprMatrix(count[,mtd$SampleID],
                                                       MIN_GENE_CPM=1, 
                                                       MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.5)
      processed.counts$filteredExprMatrix$genes
    }, COUNT)
  
  
  genesToAnalyze = unlist(genesToAnalyze) %>% 
    unique() %>%
    intersect(GENE.GC$ensembl_gene_id[!is.na(GENE.GC$percentage_gc_content)]) %>%
    intersect(GENE.LEN$ensembl_gene_id[!is.na(GENE.LEN$gene.length)]) %>%
    setdiff(c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"))
  PROCESSED_COUNTS = getGeneFilteredGeneExprMatrix(COUNT[genesToAnalyze, ], 
                                                   MIN_GENE_CPM=0, 
          MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0)
  
  dim(PROCESSED_COUNTS$filteredExprMatrix) # 18019  2008
  # Check gene biotype
  ## Define biomart object
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "dec2016.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  ## Query biomart
  Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                         filters = "ensembl_gene_id", 
                         values = PROCESSED_COUNTS$filteredExprMatrix$genes,
                         mart = mart)
  
  summary(factor(Ensemble2HGNC$gene_biotype)) %>%
    rownameToFirstColumn('Biotype') %>%
    dplyr::rename(fraction = DF) %>%
    dplyr::mutate(fraction = fraction/dim(PROCESSED_COUNTS$filteredExprMatrix$genes)[1]) %>%
    dplyr::filter(fraction >= 0.01) %>%
    kable
    
  dim(PROCESSED_COUNTS$filteredExprMatrix)[1] #18019
  dim(PROCESSED_COUNTS$filteredExprMatrix)[2] #2008
  
  ## Library Normalisation
  #Library normalisation is performed using cqn (conditional quantile normalisation)
  
  # Compute offset for gene length and gc content
  CQN.GENE_EXPRESSION = cqn(PROCESSED_COUNTS$filteredExprMatrix$counts, 
                            x = GENE.GC[PROCESSED_COUNTS$filteredExprMatrix$genes$genes, 'percentage_gc_content'],
                            lengths = GENE.LEN[PROCESSED_COUNTS$filteredExprMatrix$genes$genes, 'gene.length'],
                            lengthMethod = "smooth", 
                            verbose = FALSE)
  
  CQN.GENE_EXPRESSION$counts[1:5,1:5]
  CQN.GENE_EXPRESSION$glm.offset[1:5,1:5]
  CQN.GENE_EXPRESSION$func1[1:5,1:5]
  CQN.GENE_EXPRESSION$y[1:5,1:5]

  CQN.GENE_EXPRESSION$E = CQN.GENE_EXPRESSION$y + CQN.GENE_EXPRESSION$offset
  dim(CQN.GENE_EXPRESSION$E)# 18019  2008
  
  ### Outlier Analysis
  #### Sample outliers
  #Outlier analysis is performed before library normalisation with raw cpm counts
  
  # Find principal components of expression to plot
  PC <- prcomp(CQN.GENE_EXPRESSION$E, scale.=T, center = T)
  
  # Plot first 2 PCs
  plotdata <- data.frame(sampleIdentifier=rownames(PC$rotation), 
                         PC1=PC$rotation[,1], 
                         PC2=PC$rotation[,2]) 
  
  plotdata <- left_join(plotdata, rownameToFirstColumn(covariates, 'sampleIdentifier')) %>%
    dplyr::mutate(label = sampleIdentifier) %>%
    tidyr::separate(Tissue.Diagnosis, c('Tissue', 'Diagnosis'), sep = '\\.')
  
  plotdata$pc1_outlier=ifelse(abs(plotdata$PC1) > (mean(plotdata$PC1)+2*sd(plotdata$PC1)),plotdata$sampleIdentifier,NA)
  plotdata$pc2_outlier=ifelse(abs(plotdata$PC2) > (mean(plotdata$PC2)+2*sd(plotdata$PC2)),plotdata$sampleIdentifier,NA)
  plotdata$outlier=ifelse(is.na(plotdata$pc1_outlier)==F | is.na(plotdata$pc2_outlier)==F, plotdata$sampleIdentifier,NA)
  
  p <- ggplot(plotdata, aes(x=PC1, y=PC2))
  p <- p + geom_point(aes(color=batch, size=rin, shape = Diagnosis))
  p <- p + theme_bw() + theme(legend.position="right")
  p <- p + geom_text(aes(label= plotdata$outlier), size=4, hjust=0)
  p
  
  png("~/plot_first_2PCs_rawcounts_2SD.png", units="in", width=10, height=8, res=300)
  p
  dev.off()
  
  indToRemove=plotdata$outlier[is.na(plotdata$outlier)==F] # No outlier

  # Plot abberent distribution of logcpm counts
  tmp = covariates %>%
    tidyr::separate(Tissue.Diagnosis, c('Tissue', 'Diagnosis'), sep = '\\.') %>%
    group_by(Tissue, Diagnosis) %>%
    dplyr::summarise(count = n()) %>%
    spread(Diagnosis, count)

#### Gene outliers
#Assign NA values to genes that are above and below 3 std deviation of its distribution

# Set gene counts in specific samples that are deviating 3 sd from other samples to 3SD limit
LOG.CPM = apply(CQN.GENE_EXPRESSION$E, 1, function(x){
  mn = mean(x, na.rm = T)
  std.dev = sd(x, na.rm = T)
  
  x[x < (mn-3*std.dev)] = NA
  x[x > (mn+3*std.dev)] = NA
  return(x)
}) %>% t
CQN.GENE_EXPRESSION$E = LOG.CPM
CQN.GENE_EXPRESSION$E.no.na = CQN.GENE_EXPRESSION$E
CQN.GENE_EXPRESSION$E.no.na[is.na(CQN.GENE_EXPRESSION$E.no.na)] = 0
LIB.SIZE = colSums(PROCESSED_COUNTS$filteredExprMatrix$counts)
NEW.COUNTS = (2^LOG.CPM) * t(replicate(dim(LOG.CPM)[1], LIB.SIZE))/1e6



### Sample clustering
#PCA based clustering of samples

# Find principal components of expression to plot
PC <- prcomp(CQN.GENE_EXPRESSION$E.no.na, scale.=T, center = T)
# Plot first 2 PCs
plotdata <- data.frame(sampleIdentifier=rownames(PC$rotation), 
                       PC1=PC$rotation[,1], 
                       PC2=PC$rotation[,2])


  plotdata <- left_join(plotdata, rownameToFirstColumn(covariates, 'sampleIdentifier')) %>%
  tidyr::separate(Tissue.Diagnosis, c('Region','Diagnosis'), sep = '\\.') %>% 
  dplyr::mutate(Region = factor(Region), batch = factor(batch))

plotdata$pc1_outlier=ifelse(abs(plotdata$PC1) > (mean(plotdata$PC1)+2*sd(plotdata$PC1)),plotdata$sampleIdentifier,NA)
plotdata$pc2_outlier=ifelse(abs(plotdata$PC2) > (mean(plotdata$PC2)+2*sd(plotdata$PC2)),plotdata$sampleIdentifier,NA)
plotdata$outlier=ifelse(is.na(plotdata$pc1_outlier)==F | is.na(plotdata$pc2_outlier)==F, plotdata$sampleIdentifier,NA)


p <- ggplot(plotdata, aes(x=PC1, y=PC2))
p <- p + geom_point(aes(color=Diagnosis, size=rin, shape = sex))
p <- p + theme_bw() + theme(legend.position="right") + facet_grid(.~Region, scales = 'free_y')
p
png("~/plot_first_2PCs_cqn_2SD.png", units="in", width=10, height=8, res=300)
p
dev.off()


#Tree based classification of samples

# Eucledian tree based analysis
png("~/tree.based.classifications.png", units="in", width=10, height=8, res=300)

covariates.tmp = (covariates[,c(FactorCovariates), drop = F])
covariates.tmp[is.na(covariates.tmp)] = 0
tree = hclust(as.dist(t(CQN.GENE_EXPRESSION$E.no.na)))
cols = WGCNA::labels2colors(covariates.tmp);
WGCNA::plotDendroAndColors(tree, 
                           colors = cols, 
                           dendroLabels = FALSE, 
                           abHeight = 0.80, 
                           main = "Sample dendrogram",
                           groupLabels = colnames(covariates.tmp))

dev.off()
#Coexpression of genes 
png("~/Coexpression.of.genes.png", units="in", width=10, height=8, res=300)
cr = cor(t(CQN.GENE_EXPRESSION$E.no.na))
hist(cr, main = 'Distribution of correlation between genes', xlab = 'Correlation')
dev.off()
  
### Significant Covariates
#Correlation between pca of unadjusted mRNA expression and covariates are used to find significant covariates

# Find correlation between PC's of gene expression with covariates
preAdjustedSigCovars = runPCAandPlotCorrelations(CQN.GENE_EXPRESSION$E.no.na, 
                                                 covariates,
                                                 'NULL design(voom-normalized)', 
                                                 isKeyPlot=TRUE, 
                                                 MIN_PVE_PCT_PC = 1)


# Significant covariates to adjust at FDR 0.1 are `r preAdjustedSigCovars$significantCovars`

png("~/preadj.covariates.png", units="in", width=12, height=8, res=300)
preAdjustedSigCovars[["PC_res"]][[2]]$plotData
dev.off()


### Normalisation (iterative design)
#Since many covariates are correlated, re-normalising and re-adjusting COUNTS with an iterative design matrix.
#1. Using a mixed effect model where random effect is chosen as sampleIdentifier
#2. Adding batch and Sex a priori to variable selection
#3. Primary variable of interest Tissue.Diagnosis is excluded from the pool of available covariates for selection

# Primary variable of interest
postAdjustCovars = c("batch", "sex", "cogdx", "Tissue.Diagnosis", "Tissue.APOE4");
# Assign residual covariates
residualCovars = setdiff(preAdjustedSigCovars$significantCovars, c(postAdjustCovars, primaryVariable))
residualSigCovars = preAdjustedSigCovars
covariatesEffects = preAdjustedSigCovars$Effects.significantCovars[residualCovars]
postAdjustCovars = c(postAdjustCovars, names(which.max(covariatesEffects))) %>% unique()

loopCount = 0 
while(length(residualSigCovars$significantCovars)!=0 && loopCount <= 20){
  writeLines(paste('Using following covariates in the model:',
                   paste(postAdjustCovars, collapse=', '),
                   'as fixed effects'))
  
  # Post adjusted design matrix
  DM1 = getDesignMatrix(covariates[,postAdjustCovars,drop=F],Intercept = F)
  DM1$design = DM1$design[,linColumnFinder(DM1$design)$indepCols]
  
  # Estimate voom weights for dispersion control
  cnts = NEW.COUNTS
  cnts[is.na(cnts)] = 0
  VOOM.GENE_EXPRESSION = voom(cnts, 
                              design=DM1$design, 
                              plot=F)
  
  # Fit linear model using new weights and new design
  VOOM.ADJUSTED.FIT = lmFit(CQN.GENE_EXPRESSION$E,
                            design = DM1$design,
                            weights = VOOM.GENE_EXPRESSION$weights)
  
  
  # Residuals after normalisation
  RESIDUAL.GENE_EXPRESSION = residuals.MArrayLM(VOOM.ADJUSTED.FIT,
                                                CQN.GENE_EXPRESSION$E)

  # Residual covariates to choose from
residCovars <- setdiff(c(FactorCovariates,ContCovariates),
                       c(postAdjustCovars, primaryVariable, 'sampleIdentifier'))

  residCovars 
  # "race" "educ"
  
  # Find PC of residual gene expression and significant covariates that are highly correlated with PCs
  expr = RESIDUAL.GENE_EXPRESSION; expr[is.na(expr)] = 0; expr = expr[rowSums(expr)!=0,]
  expr[is.na(expr)] = 0
  residualSigCovars = runPCAandPlotCorrelations(expr, 
                                                covariates[, residCovars, drop=F], 
                                                'adjusted design(voom-normalized)',
                                                isKeyPlot=TRUE)
  
  # Add postadjusted covariates (if any)
  covariatesEffects = residualSigCovars$Effects.significantCovars[residCovars]
  
  postAdjustCovars = c(postAdjustCovars, names(which.max(covariatesEffects)))
  loopCount = loopCount + 1
  # print(loopCount)
  }

modelStr <- paste(paste(gsub('_','\\_',postAdjustCovars), collapse=', '),
                  'as fixed effects')
tmp <- paste('Using following covariates in the final model:', modelStr)
#"Using following covariates in the final model: cogdx, batch, sex, pct_intronic_bases, rin, pct_coding_bases, pct_intergenic_bases, pct_pf_reads_aligned, pmi, age_death, race as fixed effects


### Sanity check

# Find PC of residual gene expression and significant covariates that are highly correlated with PCs

residualSigCovars = runPCAandPlotCorrelations(expr, 
                                              covariates,
                                              'adjusted design(voom-normalized)',
                                              isKeyPlot=TRUE)

png("~/postadj.residualscovariates_and_disease.png", units="in", width=10, height=10, res=300)
residualSigCovars[["PC_res"]][[2]]$plotData
dev.off()


#PCA of residual data
# Find principal components of expression to plot
PC <- prcomp(expr, scale.=T, center = T)
# Plot first 2 PCs
plotdata <- data.frame(sampleIdentifier=rownames(PC$rotation), 
                       PC1=PC$rotation[,1], 
                       PC2=PC$rotation[,2])

plotdata <- left_join(plotdata, rownameToFirstColumn(covariates, 'sampleIdentifier')) %>%
  tidyr::separate(Tissue.Diagnosis, c('Region','Diagnosis'), sep = '\\.') %>% 
  dplyr::mutate(Region = factor(Region))

p <- ggplot(plotdata, aes(x=PC1, y=PC2))
p <- p + geom_point(aes(color=Diagnosis, size=rin, shape = sex))
p <- p + theme_bw() + theme(legend.position="right") + facet_grid(.~Region, scales = 'free_y')

png("~/plot_first_2PCs_residual_normalized_exprdata2.png", units="in", width=10, height=10, res=300)
p
dev.off()

#Tree based clustering of residual data

# Eucledian tree based analysis
png("~/tree.based.classifications.postadj.png", units="in", width=10, height=8, res=300)

covariates.tmp = data.matrix(covariates[,FactorCovariates])
covariates.tmp[is.na(covariates.tmp)] = 0
tree = hclust(as.dist(t(expr)))
cols = WGCNA::labels2colors(covariates.tmp);
WGCNA::plotDendroAndColors(tree, 
                           colors = cols, 
                           dendroLabels = FALSE, 
                           abHeight = 0.80, 
                           main = "Sample dendrogram",
                           groupLabels = colnames(covariates.tmp))

dev.off()


dim(RESIDUAL.GENE_EXPRESSION)

# Store residual gene expression profiles
RESIDUAL.GENE_EXPRESSION %>%
  rownameToFirstColumn('ensembl_gene_id') %>%
  write.table(file = 'ROSMAP_ResidualExpression_tissuediseasestatus.batch.sex.rin.pmi.age.race.fractionbases.tsv', sep = '\t', row.names=F, quote=F)

# Store metadata
metadata %>%
  write.table(file = 'ROSMAP_covariates.tsv', sep = '\t', row.names=F, quote=F)

# Store filtered counts
COUNT %>%
  rownameToFirstColumn('ensembl_gene_id') %>%
  write.table(file = 'ROSMAP_RawCounts.tsv', sep = '\t', row.names=F, quote=F)


CQN.GENE_EXPRESSION$y %>%
  rownameToFirstColumn('ensembl_gene_id') %>%
  write.table(file = 'ROSMAP_logCPM.tsv', sep = '\t', row.names=F, quote=F)


# Store cqn offsets
CQN.GENE_EXPRESSION$offset %>%
  rownameToFirstColumn('ensembl_gene_id') %>%
  write.table(file = 'ROSMAP_offset.tsv', sep = '\t', row.names=F, quote=F)


## Variance partioning plot 
varPar_ggplot=function(vp){
  vp2=(vp*100)
  names=colnames(vp2)
  t=apply(vp2, 2, FUN=median)
  t=round(t,2)
  t[t == 0.00] <- "<0.01"
  textMatrix = paste(names, "\n(", t, "%)\n", sep = "");
  melted=reshape::melt(vp2)
  ggplot(melted, aes(x=variable, y=value, color=variable)) +
    geom_violin(aes(x = variable, y = value), scale = "width") +
    geom_boxplot(aes(x = variable, y = value), width = .1,outlier.size = 0.7)+
    theme_minimal()+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.3, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.3, linetype = "solid", colour = "black"),
          axis.text.x = element_text(size=8, angle=45,colour="black",hjust=1),
          axis.text.y = element_text(size=8,colour="black"))+
    scale_x_discrete(labels=textMatrix)
}

exprs=as.matrix(RESIDUAL.GENE_EXPRESSION)
exprs[1:5,1:5]
phenotype<- metadata
phenotype$sex <- ifelse(phenotype$sex=="female", 0, 1)
phenotype$Diagnosis <- ifelse(phenotype$Diagnosis=="CONTROL",0, 1)
phenotype$Diagnosis <- as.integer(phenotype$Diagnosis)
phenotype$batch<- as.factor(phenotype$batch)
str(phenotype)
form = ~age_death + pmi + sex + (1|batch) + rin + race + pmi + apoe4+  Diagnosis+ pct_pf_reads_aligned + pct_intronic_bases + pct_intergenic_bases + pct_coding_bases 
#Replace mRNAs that have NA as counts with 0
exprs[is.na(exprs)] <- 0
exprs <- exprs[rowSums(exprs) !=0,]
exprs[1:5,1:5]
varPart <- variancePartition::fitExtractVarPartModel(exprs, form, phenotype)
vp <- variancePartition::sortCols( varPart ) 
png("../varPar_diseaseandcovariatesregressedout_2008samples.png", units="in", width=10, height=10, res=600)
varPar_ggplot(vp)
dev.off()
