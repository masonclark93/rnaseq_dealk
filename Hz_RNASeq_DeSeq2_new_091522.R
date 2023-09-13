##### September 13th 2021 #####
##### Dealkylation RNA-Seq experiment ###
#### Author: Mason Clark ###

### About the experiment:
# This in an R script for determinig differential gene expression
# of RNA-Seq data collected from Helicoverpa zea midgut and carcass tissue
# of 48H reared on Cholesterol (control), Sitosterol, Campesterol, and Stigmasterol

### The purpose of the analysis:
# Ultimately, we are interested in identifying genes associated with the process
# by which insects, such as carerpillars, convert plant sterols into cholesterol
# through dealkylation of C24 in the sterol fatty acid side chain

### The analyzes will be conducted in R using DESeq2

##### LIBRARIES, IMPORTING DATA, AND DESeqDataSetFromHTSeqCount DATA OBJECT ####

#libraries
library("DESeq2")
library("ggplot2")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("tidyverse")
library("RColorBrewer")
library("wesanderson")
library("ggsci")
library("plotly")
library("tibble")
library("AnnotationForge")
library("RSQLite")
library("clusterProfiler")
library("enrichplot")
library("ggupset")
library('ggridges')
library("org.Hzea.eg.db")
library("GenomicFeatures")
library("GenomicRanges")
library("ggbio")
library("Gviz")
#library("pathview")


#importing raw counts data from HTSeq counts
sampleFiles <- list.files(path="/Users/masonclark/Desktop/RNA-Seq Dealkylation Assay/Data/RNA_seq/sequence_data/my_analysis/new_r_analysis_071322/raw_counts", pattern="*.csv")
sampleNames <- gsub(sampleFiles, pattern=".csv", replacement="")

#directory to raw counts
dir <- ("/Users/masonclark/Dropbox/Mac/Desktop/RNA-Seq Dealkylation Assay/Data/RNA_seq/sequence_data/my_analysis/new_r_analysis_071322/raw_counts")

#Creating sampleData file
treatment = c("Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Chl",	"Cmp",	"Cmp",	"Cmp",	"Cmp",	"Cmp",	"Cmp",	"Cmp",	"Cmp",	"SCh",	"SCh",	"SCh",	"SCh",	"SCh",	"SCh",	"SCh",	"SCh",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Sit",	"Stg",	"Stg",	"Stg",	"Stg",	"Stg",	"Stg",	"Stg",	"Stg")
tissue = c("carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut",	"carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut",	"carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut",	"carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut",	"carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut",	"carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut",	"carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut",	"carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut",	"carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut",	"carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut",	"carcass",	"carcass",	"carcass",	"carcass",	"midgut",	"midgut",	"midgut",	"midgut")
time=c("0h", "0h",	"0h",	"0h",	"0h",	"0h",	"0h",	"0h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"48h",	"48h",	"48h",	"48h",	"48h",	"48h",	"48h",	"48h",	"4h",	"4h",	"4h",	"4h",	"4h",	"4h",	"4h",	"4h",	"8h",	"8h",	"8h",	"8h",	"8h",	"8h",	"8h",	"8h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"48h",	"48h",	"48h",	"48h",	"48h",	"48h",	"48h",	"48h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"4h",	"4h",	"4h",	"4h",	"4h",	"4h",	"4h",	"4h",	"8h",	"8h",	"8h",	"8h",	"8h",	"8h",	"8h",	"8h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h",	"24h")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, treatment=treatment, tissue=tissue, time=time)

#Next we need to create a DESeqDataset Object
#Design is a mathematic representation of the experiment
#Tilde means: "given"
#Here, gene expression is being evaluate given treatment (diet) and time (0-48H)
#The design is not full rank, so I will manually edit the model matrix by removing the zeroes


#Set model 
des <- formula(~treatment + time)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = dir, design=des)

##### PREFILTERING #####
#The DESeq2 vignette suggests removing rows that have approximately no information
#about the expression of a locus for a given sample. They suggest removing row where the raw count value is =< 1

#Only keep rows with counts >1
keep <-  rowSums(counts(ddsHTSeq)) > 1
dds <- ddsHTSeq[keep,]


##### GLOBAL NORMALIZATION, DISPERSION ESTs, and VARIANCE STABLIZING TRANSFORMATION #####

#Data normalization using DESeq2's median of ratios
#This is defined by counts divided by sample-specific size factors 
#determined by median ratio of gene counts relative 
#to geometric mean per gene

#Estimate size factors for each sample's gene counts
dds <- estimateSizeFactors(dds)
plotDispEsts(dds)

#Homoskedasticity - when the expected amount of variance is approximately the same across different mean values
#For RNA-Seq data, the expected variane grows with the mean
#RNA-Seq is heteroskedastic

#DeSeq2 has two transformation methods:
#Here, I will apply variance stabilization transformation for PCA
#because it is a large dataset (n=88)

#VST - dispersion-mean trend (Anders and Huber 2010)
#Applies variance stablization transformation
#transformation increased homoskedasticity for parametric PCA
#blind = FALSE; this indicates that the global variance of the data
#does not contribute to the expected variance-mean trend of exp
#setting to TRUE is an unsupervised transformation


#VST-transform
#Not directly a log-scale transformation
vsd <- vst(dds, blind = FALSE)


##### GLOBAL DATA ANALYSIS ########


#Get sample distances
sampleDists <- dist(t(assay(vsd)))

#Create matrix
sampleDistmatrix <- as.matrix(sampleDists)


#label rownames
rownames(sampleDistmatrix) <- paste(dds$treatment, dds$time, sep = '-')

#No column names
colnames(sampleDistmatrix) <- NULL

#Use R brewer packages for color gradient
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#Generate heatmap using pheatmap package
pheatmap(sampleDistmatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


##### SPLITTING DATASET BY TISSUE TYPE #####

#Subsetting midgut samples 
dds_mid <- dds[ , c(5:8, 13:16, 21:24, 29:32, 37:40, 45:48, 53:56, 61:64, 69:72, 77:80, 85:88)]

dds_car <- dds[ , c(1:4, 9:12, 17:20, 25:28, 33:36, 41:44, 49:52, 57:60, 65:68, 73:76, 81:84)]


###### MIDGUT PAIRWISE ANALYSIS (Wald Test) #####


#Estimate dispersions of midgut tissue
midgut_dispsersions <- estimateDispersions(dds_mid)

#heteroskedastic
#Mean and variance are not independent in RNA-Seq data
#variance is not continuous with mean
#as mean of norm counts increases, dispersion reduces


#Plot the dispersion estimates
plotDispEsts(dds_mid) #shows that as the mean of normalized counts increases 
#in the data, so does the variance

#black points indicate gene-wise est of dispersion
#red indicates fit of model - expected dispersion value
#blue points indicate final: raw estimates for each gene point


#Releveling dds_mid to designate cholesterol as the reference level for the model
dds_mid$treatment <- relevel(dds_mid$treatment, ref = "Chl")


#Combining treatment and time for a simpler initial design
dds_mid$group <- factor(paste0(dds_mid$treatment, dds_mid$time))

#Now restate design as function of group (this is treatment and time, but not an interaction term)
design(dds_mid) <- ~group

#DEGs analysis
#Wald test fitting:
#1 standard maximum likelihood estimates for general linearized model (GLM) coefficients
# coefficients are the same as beta, log2foldchange, effect size

dds_mid <- DESeq(dds_mid, test = "Wald", betaPrior = FALSE)

#loop for getting all results:

#vector of control samples
control <- c("Chl0h", "Chl4h", "Chl8h", "Chl24h", "Chl24h", "Chl24h", "Chl48h")

#vector of treatment smaples
trtment <- c("Sit4h", "Sit4h", "Sit8h", "Sit24h", "Cmp24h", "Stg24h", "SCh48h")


#The "BH" (aka "fdr"): controls for the false discovery rate, the expected proportion 
#of false discoveries amongst the rejected hypotheses. 
#The false discovery rate is a less stringent condition than the family-wise error rate, 
#so these methods are more powerful than the others.

#Loop for all midgut pairwise contrasts
#generates CSVs of pairwise analyses and volcano plots
for (i in 1:length(trtment)){
  contrast <- results(dds_mid, contrast=c('group', paste(trtment[i]), paste(control[i])), alpha=0.05, pAdjustMethod="fdr", parallel=TRUE)
  #Establish cut-off coniditions and add conditions to new condition column test
  padj_thres <- contrast$padj < 0.05 
  log_fold_thres <- abs(contrast$log2FoldChange) > 1.5
  contrast$test <- padj_thres & log_fold_thres
  
  #omit na values for analysis
  contrast.narm <- na.omit(contrast)
  contrast.narm <- as.data.frame(contrast.narm)
  
  #add new column label for gene_ids
  contrast.narm <- rownames_to_column(contrast.narm, var="gene_id")
  
  #Output the DEGs
  file <- paste(trtment[i], "VS", control[i], "midgut_contrast.csv", sep="_")
  write.csv(as.data.frame(contrast), file=file)
  
  #Volcano plot
  output_pdf <- paste(file, ".volcanoplot.pdf", sep="_")
  pdf(output_pdf)
  plot <- ggplot(contrast.narm, aes(x=log2FoldChange, y=-log10(padj), name = gene_id)) + 
    geom_point(aes(color=test)) + 
    scale_color_npg() +
    scale_fill_npg() +
    geom_vline(xintercept=1, color="black", linetype=2) + 
    geom_vline(xintercept=-1, color="black", linetype=2) + 
    geom_hline(yintercept=-log10(0.05), color="black", linetype=2) +
    theme_bw() +
    theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10))
  print(plot)
  dev.off()
}




###### CARCASS PAIRWISE ANALYSIS (Wald Test) #####

#Estimate dispersions of midgut tissue
carcass_dispsersions <- estimateDispersions(dds_car)

#heteroskedastic
#Mean and variance are not independent in RNA-Seq data
#variance is not continuous with mean
#as mean of norm counts increases, dispersion reduces


#Plot the dispersion estimates
plotDispEsts(dds_car) #shows that as the mean of normalized counts increases 
#in the data, so does the variance

#black points indicate gene-wise est of dispersion
#red indicates fit of model - expected dispersion value
#blue points indicate final: raw estimates for each gene point


#Releveling dds_mid to designate cholesterol as the reference level for the model
dds_car$treatment <- relevel(dds_car$treatment, ref = "Chl")


#Combining treatment and time for a simpler initial design
dds_car$group <- factor(paste0(dds_car$treatment, dds_car$time))

#Now restate design as function of group (this is treatment and time, but not an interaction term)
design(dds_car) <- ~group

#DEGs analysis
#Wald test fitting:
#1 standard maximum likelihood estimates for general linearized model (GLM) coefficients
# coefficients are the same as beta, log2foldchange, effect size

dds_car <- DESeq(dds_car, test = "Wald", betaPrior = FALSE)

#The "BH" (aka "fdr"): controls for the false discovery rate, the expected proportion 
#of false discoveries amongst the rejected hypotheses. 
#The false discovery rate is a less stringent condition than the family-wise error rate, 
#so these methods are more powerful than the others.

#Loop for all midgut pairwise contrasts
#generates CSVs of pairwise analyses and volcano plots
for (i in 1:length(trtment)){
  contrast <- results(dds_car, contrast=c('group', paste(trtment[i]), paste(control[i])), alpha=0.05, pAdjustMethod="fdr", parallel=TRUE)
  #Establish cut-off conditions and add conditions to new condition column test
  padj_thres <- contrast$padj < 0.05 
  log_fold_thres <- abs(contrast$log2FoldChange) > 1.5
  contrast$test <- padj_thres & log_fold_thres
  
  #omit na values for analysis
  contrast.narm <- na.omit(contrast)
  contrast.narm <- as.data.frame(contrast.narm)
  
  #add new column label for gene_ids
  contrast.narm <- rownames_to_column(contrast.narm, var="gene_id")
  
  #Output the DEGs
  file <- paste(trtment[i], "VS", control[i], "carcass_contrast.csv", sep="_")
  write.csv(as.data.frame(contrast), file=file)
  
  #Volcano plot
  output_pdf <- paste(file, ".volcanoplot.pdf", sep="_")
  pdf(output_pdf)
  plot <- ggplot(contrast.narm, aes(x=log2FoldChange, y=-log10(padj), name = gene_id)) + 
    geom_point(aes(color=test)) + 
    scale_color_npg() +
    scale_fill_npg() +
    geom_vline(xintercept=1, color="black", linetype=2) + 
    geom_vline(xintercept=-1, color="black", linetype=2) + 
    geom_hline(yintercept=-log10(0.05), color="black", linetype=2) +
    theme_bw() +
    theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10))
  print(plot)
  dev.off()
}

###### CONTROL CHOLESTEROL GROUP PAIRWISE CONTRASTS #######

#This section is for identifying differentially expressed genes in both tissues
#in control group (cholesterol) individuals only.

#All differentially expressed genes in the control groups will be identified
#and removed from the pairwise analyses of the treatment groups because
#these genes are not associated with a treatment effect - only time (i.e. circadian clock and development)

#control list

control1 <- c("Chl0h", "Chl4h", "Chl8h", "Chl24h")
control2 <- c("Chl4h", "Chl8h", "Chl24h", "Chl48h")

#Loop for all midgut control pairwise contrasts
#generates CSVs of pairwise analyses and volcano plots
for (i in 1:length(control2)){
  contrast <- results(dds_mid, contrast=c('group', paste(control2[i]), paste(control1[i])), alpha=0.05, pAdjustMethod="fdr", parallel=TRUE)
  #Establish cut-off coniditions and add conditions to new condition column test
  padj_thres <- contrast$padj < 0.05 
  log_fold_thres <- abs(contrast$log2FoldChange) > 1.5
  contrast$test <- padj_thres & log_fold_thres
  
  #omit na values for analysis
  contrast.narm <- na.omit(contrast)
  contrast.narm <- as.data.frame(contrast.narm)
  
  #add new column label for gene_ids
  contrast.narm <- rownames_to_column(contrast.narm, var="gene_id")
  
  #Output the DEGs
  file <- paste(control2[i], "VS", control1[i], "midgut_controls_contrast.csv", sep="_")
  write.csv(as.data.frame(contrast), file=file)
}

#Loop for all carcass control pairwise contrasts
#generates CSVs of pairwise analyses and volcano plots
for (i in 1:length(control2)){
  contrast <- results(dds_car, contrast=c('group', paste(control2[i]), paste(control1[i])), alpha=0.05, pAdjustMethod="fdr", parallel=TRUE)
  #Establish cut-off coniditions and add conditions to new condition column test
  padj_thres <- contrast$padj < 0.05 
  log_fold_thres <- abs(contrast$log2FoldChange) > 1.5
  contrast$test <- padj_thres & log_fold_thres
  
  #omit na values for analysis
  contrast.narm <- na.omit(contrast)
  contrast.narm <- as.data.frame(contrast.narm)
  
  #add new column label for gene_ids
  contrast.narm <- rownames_to_column(contrast.narm, var="gene_id")
  
  #Output the DEGs
  file <- paste(control2[i], "VS", control1[i], "carcass_controls_contrast.csv", sep="_")
  write.csv(as.data.frame(contrast), file=file)
}

########## CLUSTERPROFILER FOR GSEA AND KEGG ANALYSIS ############

#For AnnotationForge
#My very own Hz annotation database
organism <- org.Hzea.eg.db

#For KeggAnalysis
kegg_organism <- "hze"

#Let's get a list of the midgutfiles
midgutFiles <- list.files(path="~/Desktop/RNA-Seq Dealkylation Assay/Data/RNA_seq/sequence_data/my_analysis/new_r_analysis_071322/midgut_DEGSs_out/csvs", pattern="*.csv")

#Now all of the carcass files
carcassFiles <- list.files(path="~/Desktop/RNA-Seq Dealkylation Assay/Data/RNA_seq/sequence_data/my_analysis/new_r_analysis_071322/carcass_DEGS_out/csvs", pattern="*.csv")

#Parsing files function
file_degs_parse <- function(file){
  #Read in file
  dat <- read.csv(file, header=TRUE)
  
  #GSEA/KEGG list
  gsea_gene_list <- dat$log2FoldChange
  names(gsea_gene_list) <- dat$gene_ID
  gsea_gene_list <- na.omit(gsea_gene_list)
  gsea_gene_list <- sort(gsea_gene_list, decreasing=TRUE)
  
  #GO Term analysis
  #Need sig genes only
  sig_genes <- subset(dat, padj < 0.05)
  genes <- sig_genes$log2FoldChange
  names(genes) <- sig_genes$gene_ID
  genes <- na.omit(genes)
  
  #Upregulated
  upreg_genes_list <- names(genes)[genes > 1.5]
  
  #Downregulated
  downreg_genes_list <- names(genes)[genes < 1.5]
  
  #Turn outputs into lists
  list <- list(gsea_gene_list, upreg_genes_list, downreg_genes_list) 
  return(list)
}

for(i in 1:length(midgutFiles)){
  dat <- file_degs_parse(midgutFiles[i])
  assign(paste(gsub(midgutFiles[i], pattern="contrast.csv", replacement = "", fixed=TRUE),"_parsed", sep=""), dat)
}

for(i in 1:length(carcassFiles)){
  dat <- file_degs_parse(carcassFiles[i])
  assign(paste(gsub(carcassFiles[i], pattern="contrast.csv", replacement = "", fixed=TRUE),"_parsed", sep=""), dat)
}


#Function for generating GSEA
gsea_figures <- function(list){
  gene_list <- list[[1]]
  gse <- gseGO(geneList=gene_list, 
        ont ="ALL", 
        keyType = "GID", 
        minGSSize = 3, 
        maxGSSize = 800, 
        pvalueCutoff = 0.05, 
        verbose = TRUE, 
        OrgDb = org.Hzea.eg.db, 
        pAdjustMethod = "fdr")
  gse <- pairwise_termsim(gse)
  dotplot <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  emapplot <- emapplot(gse, showCategory = 10)
  cnetplot <- cnetplot(gse, foldChange=gene_list, circular = TRUE, colorEdge=TRUE)
  ridgeplot <- ridgeplot(gse) + labs(x = "enrichment distribution")
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "gsea_dotplot"), dotplot, envir=.GlobalEnv)
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "gsea_emapplot"), emapplot, envir=.GlobalEnv)
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "gsea_cnetplot"), cnetplot, envir=.GlobalEnv)
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "gsea_ridgeplot"), ridgeplot, envir=.GlobalEnv)
  return()
}


#Midgut GSEA Figures
gsea_figures(Cmp24h_VS_Chl24h_midgut__parsed)
gsea_figures(Sit24h_VS_Chl24h_midgut__parsed)
gsea_figures(Sit4h_VS_Chl4h_midgut__parsed)
gsea_figures(Stg24h_VS_Chl24h_midgut__parsed)
gsea_figures(SCh48h_VS_Chl48h_midgut__parsed)
gsea_figures(Sit8h_VS_Chl8h_midgut__parsed)

#Carcass GSEA Figures
gsea_figures(Cmp24h_VS_Chl24h_carcass__parsed)
gsea_figures(Sit24h_VS_Chl24h_carcass__parsed)
gsea_figures(Sit4h_VS_Chl4h_carcass__parsed)
gsea_figures(Stg24h_VS_Chl24h_carcass__parsed)
gsea_figures(SCh48h_VS_Chl48h_carcass__parsed)
gsea_figures(Sit8h_VS_Chl8h_carcass__parsed)

#Now KEGG Pathway analysis
#Function for generating GSEA
kegg_figures <- function(list){
  gene_list <- list[[1]]
  gene_list <- unlist(gene_list)
  kk2 <- gseKEGG(geneList     = gene_list,
                 organism     = kegg_organism,
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr",
                 keyType       = "ncbi-geneid")
  kk2 <- pairwise_termsim(kk2)
  dotplot <- dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
  emapplot <- emapplot(kk2)
  cnetplot <- cnetplot(kk2, foldChange=gene_list, circular = TRUE, colorEdge=TRUE)
  ridgeplot <- ridgeplot(kk2) + labs(x = "enrichment distribution")
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "kegg_dotplot"), dotplot, envir=.GlobalEnv)
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "kegg_emapplot"), emapplot, envir=.GlobalEnv)
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "kegg_cnetplot"), cnetplot, envir=.GlobalEnv)
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "kegg_ridgeplot"), ridgeplot, envir=.GlobalEnv)
  return()
}

#Midgut KEGG Figures
kegg_figures(Cmp24h_VS_Chl24h_midgut__parsed)
kegg_figures(Sit24h_VS_Chl24h_midgut__parsed)
kegg_figures(Sit4h_VS_Chl4h_midgut__parsed)
kegg_figures(Stg24h_VS_Chl24h_midgut__parsed)
kegg_figures(SCh48h_VS_Chl48h_midgut__parsed)
kegg_figures(Sit8h_VS_Chl8h_midgut__parsed)

#carcass KEGG Figures
kegg_figures(Cmp24h_VS_Chl24h_carcass__parsed)
kegg_figures(Sit24h_VS_Chl24h_carcass__parsed)
kegg_figures(Sit4h_VS_Chl4h_carcass__parsed)
kegg_figures(Stg24h_VS_Chl24h_carcass__parsed)
kegg_figures(SCh48h_VS_Chl48h_carcass__parsed)
kegg_figures(Sit8h_VS_Chl8h_carcass__parsed)


########## CLUSTERPROFILER For GO Enrichment ############

#GO functions
enrGO <- function(list){
  universe <- list[[1]]
  genes_up <- list[[2]]
  genes_down <- list[[3]]
  enGO.Up <- enrichGO(gene = genes_up,
           universe = names(universe),
           OrgDb = org.Hzea.eg.db, 
           keyType = "GID",
           readable = F,
           ont = "ALL",
           pvalueCutoff = 0.05, 
           qvalueCutoff = 0.10)
  enGO.Dn <- enrichGO(gene = genes_down,
                      universe = names(universe),
                      OrgDb = org.Hzea.eg.db, 
                      keyType = "GID",
                      readable = F,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
  dtplt1 <- dotplot(enGO.Up)
  dtplt2 <- dotplot(enGO.Dn)
  upsetty1 <- upsetplot(enGO.Up)
  upsetty2 <- upsetplot(enGO.Dn)
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "Up_GOs_dotplot"), dtplt1, envir=.GlobalEnv)
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "Down_GOs_dotplot"), dtplt2, envir=.GlobalEnv)
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "Up_GOs_upsetplot"), upsetty1, envir=.GlobalEnv)
  assign(paste(gsub(print(substitute(list)), pattern="_parsed",replacement=""), sep ="", "Down_GOs_upsetplot"), upsetty2, envir=.GlobalEnv)
}

#Midgut GO Enrichment
enrGO(Cmp24h_VS_Chl24h_midgut__parsed) #NO output
enrGO(Sit24h_VS_Chl24h_midgut__parsed) #NO output
enrGO(Sit4h_VS_Chl4h_midgut__parsed)
enrGO(Stg24h_VS_Chl24h_midgut__parsed) #NO output
enrGO(SCh48h_VS_Chl48h_midgut__parsed) #NO output
enrGO(Sit8h_VS_Chl8h_midgut__parsed)

#carcass GO Enrichment
enrGO(Cmp24h_VS_Chl24h_carcass__parsed) #NO output
enrGO(Sit24h_VS_Chl24h_carcass__parsed)
enrGO(Sit4h_VS_Chl4h_carcass__parsed)
enrGO(Stg24h_VS_Chl24h_carcass__parsed) #NO output
enrGO(SCh48h_VS_Chl48h_carcass__parsed) 
enrGO(Sit8h_VS_Chl8h_carcass__parsed)



########## NORAMLIZATION ########

#Exporting normalized counts for later analyses

#midgut
normCounts_mid <- counts(dds_mid, normalized = T)

#midgut normalized counts written to new csv file
write.csv(as.data.frame(normCounts_mid), file='midgut_DEGs_out/nornmalized_midgut_counts.csv')

#carcass
normCounts_car <- counts(dds_car, normalized = T)

#carcass normalized counts written to new csv file
write.csv(as.data.frame(normCounts_car), file='carcass_DEGs_out/nornmalized_carcass_counts.csv')

########## GENOMICFEATURES for CHROMOSOME-LEVEL ANALYSIS #########3

dir1 <- ("/Users/masonclark/Dropbox/Mac/Desktop/RNA-Seq Dealkylation Assay/Data/RNA_seq/sequence_data/my_analysis/new_r_analysis_071322/")


#Need to first create a data.frame of supplemental information for the chromosome
chr.info <- data.frame(chrom = as.character(c("NC_061452.1",
                                              "NC_061453.1",
                                              "NC_061454.1",
                                              "NC_061455.1",
                                              "NC_061456.1",
                                              "NC_061457.1",
                                              "NC_061458.1",
                                              "NC_061459.1",
                                              "NC_061460.1",
                                              "NC_061461.1",
                                              "NC_061462.1",
                                              "NC_061463.1",
                                              "NC_061464.1",
                                              "NC_061465.1",
                                              "NC_061466.1",
                                              "NC_061467.1",
                                              "NC_061468.1",
                                              "NC_061469.1",
                                              "NC_061470.1",
                                              "NC_061471.1",
                                              "NC_061472.1",
                                              "NC_061473.1",
                                              "NC_061474.1",
                                              "NC_061475.1",
                                              "NC_061476.1",
                                              "NC_061477.1",
                                              "NC_061478.1",
                                              "NC_061479.1",
                                              "NC_061480.1",
                                              "NC_061481.1",
                                              "NC_061482.1",
                                              "NW_025899711.1",
                                              "NW_025899712.1",
                                              "NW_025899713.1",
                                              "NW_025899714.1",
                                              "NW_025899715.1",
                                              "NW_025899716.1",
                                              "NW_025899717.1",
                                              "NW_025899718.1",
                                              "NW_025899719.1",
                                              "NW_025899720.1",
                                              "NW_025899721.1",
                                              "NC_061507.1")),
                            length = as.integer(c(15512169,
                                                  15061584,
                                                  14695159,
                                                  14165716,
                                                  14051263,
                                                  14023251,
                                                  13649923,
                                                  13649271,
                                                  13615264,
                                                  13258385,
                                                  13088840,
                                                  12972780,
                                                  12841304,
                                                  12757947,
                                                  12605921,
                                                  12332656,
                                                  12254882,
                                                  12245473,
                                                  11629846,
                                                  11599754,
                                                  11461517,
                                                  10871041,
                                                  10422137,
                                                  9742508,
                                                  9280487,
                                                  9278413,
                                                  7355033,
                                                  7211800,
                                                  7203478,
                                                  6316813,
                                                  18805280,
                                                  500001,
                                                  250000,
                                                  125021,
                                                  88156,
                                                  51777,
                                                  47286,
                                                  42317,
                                                  34946,
                                                  25017,
                                                  20924,
                                                  20253,
                                                  15351)), 
                                            is_circular = as.vector(c(rep(FALSE,42), TRUE)), stringsAsFactors = FALSE)

#First need to create a txdb object of Hz genome from gtf
txHz <- makeTxDbFromGFF(file = "GCF_022581195.2_ilHelZeax1.1_genomic.gtf", 
                        chrominfo = chr.info,
                        organism = "Helicoverpa zea",
                        taxonomyId = 7113,
                        format = "gtf")

#Ran some console commands just to test out the db and after some struggling, things seem to check out

#First step is to create a genomic ranges object of gene start and stop positions
exon.ranges <- exonsBy(txHz, "gene") %>%
  range() %>% #ranges reduces it down to one start and stop position
  unlist()

exon.ranges <- keepSeqlevels(exon.ranges, value = c("NC_061452.1","NC_061453.1","NC_061454.1","NC_061455.1",
                                    "NC_061456.1","NC_061457.1","NC_061458.1","NC_061459.1",
                                    "NC_061460.1","NC_061461.1","NC_061462.1","NC_061463.1",
                                    "NC_061464.1","NC_061465.1","NC_061466.1","NC_061467.1",
                                    "NC_061468.1","NC_061469.1","NC_061470.1","NC_061471.1",
                                    "NC_061472.1","NC_061473.1","NC_061474.1","NC_061475.1",
                                    "NC_061476.1","NC_061477.1","NC_061478.1","NC_061479.1",
                                    "NC_061480.1","NC_061481.1","NC_061482.1"), pruning.mode="tidy")

#rename chromosomes
chromosome_order <- c("chr1",
                      "chr2",
                      "chr3",
                      "chr4",
                      "chr5",
                      "chr6",
                      "chr7",
                      "chr8",
                      "chr9",
                      "chr10",
                      "chr11",
                      "chr12",
                      "chr13",
                      "chr14",
                      "chr15",
                      "chr16",
                      "chr17",
                      "chr18",
                      "chr19",
                      "chr20",
                      "chr21",
                      "chr22",
                      "chr23",
                      "chr24",
                      "chr25",
                      "chr26",
                      "chr27",
                      "chr28",
                      "chr29",
                      "chr30",
                      "chr31")


exon.ranges <- renameSeqlevels(exon.ranges, chromosome_order)


#Let's test combining metadata
#midgut tissue for sitosterol at 8 hours
dat <- read.csv("midgut_DEGSs_out/csvs/Sit8h_VS_Chl8h_midgut_contrast.csv", header=TRUE)

#Get sig genes only (FDR < 0.05)
sigGenes <- subset(dat, padj<0.05)

#I needed to append "LOC" to the sigGenes data.frame 
#to match the names() in exon.ranges
sigGenes$gene_ID <- paste0("LOC",sigGenes$gene_ID)

#Only keep genes in exon.ranges that matches the sigGene gene_IDs
sigRegions <- exon.ranges[na.omit(match(sigGenes$gene_ID, names(exon.ranges)))]

#Only evaluate expression across the 31 nuclear chromosomes
sigRegions <- keepSeqlevels(sigRegions, value = c("NC_061452.1","NC_061453.1","NC_061454.1","NC_061455.1",
                                                   "NC_061456.1","NC_061457.1","NC_061458.1","NC_061459.1",
                                                   "NC_061460.1","NC_061461.1","NC_061462.1","NC_061463.1",
                                                   "NC_061464.1","NC_061465.1","NC_061466.1","NC_061467.1",
                                                   "NC_061468.1","NC_061469.1","NC_061470.1","NC_061471.1",
                                                   "NC_061472.1","NC_061473.1","NC_061474.1","NC_061475.1",
                                                   "NC_061476.1","NC_061477.1","NC_061478.1","NC_061479.1",
                                                   "NC_061480.1","NC_061481.1","NC_061482.1"), pruning.mode="tidy")

#rename chromosomes
chromosome_order <- c("1",
                      "2",
                      "3",
                      "4",
                      "5",
                      "6",
                      "7",
                      "8",
                      "9",
                      "10",
                      "11",
                      "12",
                      "13",
                      "14",
                      "15",
                      "16",
                      "17",
                      "18",
                      "19",
                      "20",
                      "21",
                      "22",
                      "23",
                      "24",
                      "25",
                      "26",
                      "27",
                      "28",
                      "29",
                      "30",
                      "31")
#rename chromosomes
sigRegions <- renameSeqlevels(sigRegions, chromosome_order)

#Combine the exon-ranges data with the DGE metadata
mcols(sigRegions) <- sigGenes[match(names(sigRegions), sigGenes$gene_ID), ]

#Define Upregulated genes
mcols(sigRegions)$UpRegulated <- mcols(sigRegions)$log2FoldChange > 1.5

#karyogram plotting
autoplot(sigRegions, layout="karyogram", aes(color=UpRegulated, fill=UpRegulated)) +
  scale_fill_discrete(name = "Gene Regulation", labels = c("Upregulated", "Downregulated"))

  
#Onto Gviz for visualizing expression at gene level


#I want to create a granges object that has the genomic coordinates of all of the genes
#Need to fix this first
dat$gene_ID <- paste0("LOC",dat$gene_ID)

#Only keep genes in exon.ranges that matches the gene_IDs
allgenes <- exon.ranges[na.omit(match(dat$gene_ID, names(exon.ranges)))]

#Combine
#Combine the exon-ranges data with the DGE metadata
mcols(allgenes) <- dat[match(names(allgenes), dat$gene_ID), ]

#color code for significance
status <- factor(ifelse(allgenes$padj < 0.05 & abs(allgenes$log2FoldChange) > 1.5 & !is.na(allgenes$padj),
                        "sig", "notsig"))

allgenes$group <- names(allgenes)
 

#Build the annotation track
#atrack shows all DEGs and their orientation
#Each chromosome is displayed separately
#I used chromosome(atrack) <- '25' to visualize just chr25
atrack <- AnnotationTrack(allgenes, name="Genes Range", group=allgenes$gene_ID, ucscChromosomeNames=FALSE, feature=status)

#Now I will add the genomic coordinates
gtrack  <- GenomeAxisTrack()
#Datatrack
dtrack <- DataTrack(allgenes, data="log2FoldChange", name = "log2FC", strand = "*", type='histogram')

#Plotting track surrounding region of HzDHCR24-1
#chr25
#~38mb position
plotTracks(list(dtrack, gtrack, atrack), chromosome='chr25', groupAnnotation = 'group', sig="hotpink", notsig="grey", from=3710000, to=3900000)

#Plotting track surrounding region of delta(7)-sterol 5(6)-desaturase erg32-like
#chr7
#coordinates 474019	482300
#LOC124631934
#adjacent gene is LOC124632161 (uncharacterized)
plotTracks(list(dtrack, gtrack, atrack), chromosome='chr7', groupAnnotation = 'group', sig="hotpink", notsig="grey", from=424019, to=512300)

#chr19
#positon 9939144	9941330
#A cluster of CYP50 6B5-likes
plotTracks(list(dtrack, gtrack, atrack), chromosome='chr19', groupAnnotation = 'group', sig="hotpink", notsig="grey", from=9890000, to=9970000)

#NPC intracellular cholesterol transporter 1 isoform X1 
#chr14
#positions 7723066	7725702
plotTracks(list(dtrack, gtrack, atrack), chromosome='chr14', groupAnnotation = 'group', sig="hotpink", notsig="grey", from=7715000, to=7750000)


#chr19 cluster
#not very interesting
plotTracks(list(dtrack, gtrack, atrack), chromosome='chr19', groupAnnotation = 'group', sig="hotpink", notsig="grey", from=5910000, to=6500000)

#DHCR24-2
#chr17
#positions 520523	523198
plotTracks(list(dtrack, gtrack, atrack), chromosome='chr17', groupAnnotation = 'group', sig="hotpink", notsig="grey", from=490000, to=710000)


#chr8 
#CYP450 4C1
#LOC124632321
#Position 3760574	3775015
plotTracks(list(dtrack, gtrack, atrack), chromosome='chr8', groupAnnotation = 'group', sig="hotpink", notsig="grey", from=3660574, to=3875015)
###NOT VERY INTERESTING


#chr24
#JH epoxide hydrolase
#LOC124642336
#position 1684567	1693568
plotTracks(list(dtrack, gtrack, atrack), chromosome='chr24', groupAnnotation = 'group', sig="hotpink", notsig="grey", from=1670000, to=1720000)


############## LIKLIHOOD RATIO TEST FOR SITOSTEROL SUBSET ###############

#importing raw counts data from HTSeq counts
sampleFiles1 <- list.files(path="~/Dropbox/Mac/Desktop/RNA-Seq Dealkylation Assay/Data/RNA_seq/sequence_data/my_analysis/new_r_analysis_071322/Sit_v_Chl", pattern="*.csv")
sampleNames1 <- gsub(sampleFiles1, pattern="_paired_kraken2_unclassified_out_1.fq.Aligned.sortedByCoord.out.bam.reverse_stranded.csv", replacement="")

#directory
dirsit = c("~/Dropbox/Mac/Desktop/RNA-Seq Dealkylation Assay/Data/RNA_seq/sequence_data/my_analysis/new_r_analysis_071322/Sit_v_Chl")

#Creating sampleData file
treatment1 = c(rep("Chl", 24), rep("Sit", 24))
tissue1 = c(rep(c(rep("carcass", 4), rep("midgut", 4)), 6))
time1=c(rep(c(rep("24h",8), rep("4h",8), rep("8h", 8)),2))
sampleTable1 <- data.frame(sampleName = sampleNames1, fileName = sampleFiles1, treatment=treatment1, tissue=tissue1, time=time1)

#Next we need to create a DESeqDataset Object
#Design is a mathematic representation of the experiment
#Tilde means: "given"
#Here, gene expression is being evaluate given treatment (diet) and time (0-48H)
#The design is not full rank, so I will manually edit the model matrix by removing the zeroes


#Set model 
des1 = formula(~treatment + time + treatment:time)
ddsHTSeq1 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable1, directory = dirsit, design=des1)


##### PREFILTERING #####
#The DESeq2 vignette suggests removing rows that have approximately no information
#about the expression of a locus for a given sample. They suggest removing row where the raw count value is =< 1

#Only keep rows with counts >1
keep1 <-  rowSums(counts(ddsHTSeq1)) > 1
dds1 <- ddsHTSeq1[keep1,]

#relevel
dds1$treatment <- factor(dds1$treatment, levels = c("Chl","Sit"))


#subset by tissue type
dds_mid1 <- dds1[ , c(5:8, 13:16, 21:24, 29:32, 37:40, 45:48)]
dds_car1 <- dds1[ , c(1:4, 9:12, 17:20, 25:28, 33:36, 41:44)]

#Analysis
dds_mid_sitchl <- DESeq(dds_mid1, test="LRT", full = des1, reduced = ~treatment + time)
dds_car_sitchl <- DESeq(dds_car1, test="LRT", full = des1, reduced = ~treatment + time)

#Get results
sitchl_midgut_results <- results(dds_mid_sitchl)
sitchl_carcass_results <- results(dds_car_sitchl)

#Create cut-off column to distinguish DEGs (midgut)
padj_thres <- sitchl_midgut_results$padj < 0.05 
log_fold_thres <- abs(sitchl_midgut_results$log2FoldChange) > 1.5
sitchl_midgut_results$test <- padj_thres & log_fold_thres

#Create cut-off column to distinguish DEGs (carcass)
padj_thres <- sitchl_carcass_results$padj < 0.05 
log_fold_thres <- abs(sitchl_carcass_results$log2FoldChange) > 1.5
sitchl_carcass_results$test <- padj_thres & log_fold_thres

#Export results
write.csv(as.data.frame(sitchl_midgut_results), file='Sit_versus_Chl_midgut_full_rank.csv')
write.csv(as.data.frame(sitchl_carcass_results), file='Sit_versus_Chl_carcass_full_rank.csv')


###### SIT SUBSET GSEA, KEGG, GO ANALYSES ######
#Using all same code from above pairwise analyses


#files
sit.midgut <- c("/Users/masonclark/Dropbox/Mac/Desktop/RNA-Seq Dealkylation Assay/Data/RNA_seq/sequence_data/my_analysis/new_r_analysis_071322/Sit_v_Chl/csvs/Sit_versus_Chl_midgut_full_rank.csv")
sit.carcass <- c("/Users/masonclark/Dropbox/Mac/Desktop/RNA-Seq Dealkylation Assay/Data/RNA_seq/sequence_data/my_analysis/new_r_analysis_071322/Sit_v_Chl/csvs/Sit_versus_Chl_carcass_full_rank.csv")

sit.midgut.parsed <- file_degs_parse(sit.midgut)
sit.carcass.parsed <- file_degs_parse(sit.carcass)

go.test <- gseGO(geneList=sit.midgut.parsed[[1]], 
             ont ="BP", 
             keyType = "GID", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hzea.eg.db, 
             pAdjustMethod = "fdr")


#GSEA
gsea_figures(sit.midgut.parsed)
gsea_figures(sit.carcass.parsed)

#KEGG
kegg_figures(sit.midgut.parsed)
kegg_figures(sit.carcass.parsed)

#GO
enrGO(sit.midgut.parsed)
enrGO(sit.carcass.parsed)

