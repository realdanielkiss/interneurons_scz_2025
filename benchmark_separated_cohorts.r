#import packages
library(GEOquery)
library(dtangle)
#library(hgu133plus2.db)
#library(AnnotationDbi)
library(limma)
library(ggplot2)
library(reshape2)
library(edgeR)
library(matrixStats)
library(dplyr)
library(tidyverse)
library(matrixStats)
library(broom)
library(knitr)
library(ggpubr)
library(biomaRt)
library(ggrepel)
library(patchwork)
library(ggsignif)
library(modelr)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
library(markerGeneProfile) 
library(DESeq2)
library(factoextra)
library(PCAtools)
library(ggplot2)
library(readr)


#set working directory
#setwd("/nethome/kcni/xzhou/cell_deconv/data")

### Dan's new MGP pipeline:
### Load and format matrices
GVEX_matrix = read.csv("/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/GVEX_count_matrix.csv")
names(GVEX_matrix) = gsub("X", '', names(GVEX_matrix))
names(GVEX_matrix)[1] = "gene_symbol"
LIBD_matrix = read.csv("/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/LIBD_count_matrix.csv")
names(LIBD_matrix) = gsub("X", '', names(LIBD_matrix))
names(LIBD_matrix)[1] = "gene_symbol"
CMC_matrix = read.csv("/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/CMC_count_matrix.csv")
names(CMC_matrix) = gsub("X", '', names(CMC_matrix))
names(CMC_matrix)[1] = "gene_symbol"

#Rename EnsemblIDs to be compatible with gene names
GVEX_matrix$gene_symbol = gsub("\\..*", '', GVEX_matrix$gene_symbol)
LIBD_matrix$gene_symbol = gsub("\\..*", '', LIBD_matrix$gene_symbol)
CMC_matrix$gene_symbol = gsub("\\..*", '', CMC_matrix$gene_symbol)

#Make EnsemblIDs row names
row.names(GVEX_matrix) = GVEX_matrix$gene_symbol
GVEX_matrix = GVEX_matrix[,-1]
row.names(LIBD_matrix) = LIBD_matrix$gene_symbol
LIBD_matrix = LIBD_matrix[,-1]
row.names(CMC_matrix) = CMC_matrix$gene_symbol
CMC_matrix = CMC_matrix[,-1]

#Combine into one large matrix for easier analysis
combined_matrix = cbind(GVEX_matrix, LIBD_matrix, CMC_matrix)

###########################

###Load and format metadata
psychencode_metadata = read.csv(("//external/rprshnas01/external_data/psychencode/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv"))

GVEX_metadata = read.delim("/external/rprshnas01/external_data/psychencode/PsychENCODE//BrainGVEX/RNAseq/SYNAPSE_METADATA_MANIFEST.tsv") %>% filter(dataType == "geneExpression")  %>% left_join(psychencode_metadata)
LIBD_metadata = read.delim("/external/rprshnas01/external_data/psychencode/PsychENCODE/LIBD__szControl/RNAseq/SYNAPSE_METADATA_MANIFEST.tsv") %>% filter(dataType == "geneExpression")  %>% left_join(psychencode_metadata)
CMC_metadata = read.csv("/external/rprshnas01/external_data/psychencode/PsychENCODE/CMC/Metadata/SYNAPSE_TABLE_QUERY_123020650.csv") %>% filter(dataType == "geneExpression", fileFormat == "tsv")  %>% left_join(psychencode_metadata)
names(CMC_metadata)[names(CMC_metadata) == "Individual_ID"] = "individualID"

#Change LIBD NAs to LIBD_szControl
LIBD_metadata[which(is.na(LIBD_metadata$individualIdSource)), 'individualIdSource'] = 'LIBD_szControl'

columns = names(GVEX_metadata) %>% intersect(names(LIBD_metadata)) %>% intersect(names(CMC_metadata))
#Combine all metadata, filter based on inclusion criteria, and add newStudy column for cohort-based MGP
combined_metadata = rbind(GVEX_metadata[, columns], LIBD_metadata[, columns], CMC_metadata[, columns])%>% filter(!is.na(individualIdSource), primaryDiagnosis %in% c("control", "Schizophrenia"), ageDeath >= 15)
#Additional filtering to remove duplicates from Pitt and GVEX
combined_metadata = combined_metadata %>% filter(!grepl("_BP_", specimenID), contributingStudy != "[\"BrainGVEX\"]", contributingStudy != "[\"LIBD_szControl\"]")
#Create newStudy column for filtering
combined_metadata = combined_metadata %>%
  mutate(newStudy = case_when(
    contributingStudy == "LIBD_szControl"           ~ "LIBD_szControl",
    contributingStudy == "[\"CMC_HBCC\"]"           ~ "NIMH_HBCC",
    grepl("SMRI", individualIdSource)               ~ "GVEX",
    individualIdSource == "MSSM"                    ~ "MSSM",
    individualIdSource == "Penn"                    ~ "Penn",
    individualIdSource == "Pitt"                    ~ "Pitt",
    TRUE                                            ~ "Not_Used"
  )) %>% filter(newStudy != "Not_Used")
############################################

###Subset matrices based on study cohort
#GVEX and LIBD matrices are named based on sampleID, whereas CMC is named based on individualID
#Create standardID to pull samples from the same column
combined_metadata = combined_metadata %>%
  mutate(standardID = ifelse(newStudy %in% c("MSSM", "Penn", "Pitt"), individualID, specimenID)) %>%
  mutate(standardID = ifelse(newStudy == "GVEX", 
                             str_replace_all(standardID, "-", "."), 
                             standardID))
combined_metadata = combined_metadata %>%
  mutate(dataset = case_when(
    newStudy %in% c("NIMH_HBCC", "MSSM", "Pitt", "Penn") ~ "CMC",
    newStudy == "GVEX"                                   ~ "GVEX",
    newStudy == "LIBD_szControl"                         ~ "LIBD_szControl",
    TRUE                                                 ~ NA_character_
  ))

#Create cohort specific matrices
libd_samples = combined_metadata %>% filter(newStudy == "LIBD_szControl") %>% pull(standardID)
nimh_samples = combined_metadata %>% filter(newStudy == "NIMH_HBCC") %>% pull(standardID)
gvex_samples = combined_metadata %>% filter(newStudy == "GVEX") %>% pull(standardID) 
mssm_samples = combined_metadata %>% filter(newStudy == "MSSM") %>% pull(standardID)
penn_samples = combined_metadata %>% filter(newStudy == "Penn") %>% pull(standardID)
pitt_samples = combined_metadata %>% filter(newStudy == "Pitt") %>% pull(standardID)
#Create matrices
gene_symbol = rownames(combined_matrix)
libd_matrix = combined_matrix[, colnames(combined_matrix) %in% libd_samples] 
nimh_matrix = combined_matrix[, colnames(combined_matrix) %in% nimh_samples] 
gvex_matrix = combined_matrix[, colnames(combined_matrix) %in% gvex_samples] 
mssm_matrix = combined_matrix[, colnames(combined_matrix) %in% mssm_samples] 
penn_matrix = combined_matrix[, colnames(combined_matrix) %in% penn_samples]
pitt_matrix = combined_matrix[, colnames(combined_matrix) %in% pitt_samples] 
############################################

### Cell-Type Proportion Estimation: Updated marker list from Micaela's paper
### GET MARKERS FOR MGP ANALYSIS
# note that this is the list of markers from micaela's paper - you get similar but diff results if you use the markers from the aging paper
sonny_markers = read_csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'))
colnames(sonny_markers) = colnames(sonny_markers) %>% make.names() %>% tolower()

hgnc_mapping = read_tsv('/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/hgnc_complete_set.txt')

# now, this is the list of sonnys markers with entrez ids and ensembl ids where possible
sonny_hgnc_merged_markers = left_join(sonny_markers %>% dplyr::rename(entrez_id = entrez.gene.id), 
                                      hgnc_mapping %>% distinct(entrez_id, .keep_all = T)%>% 
                                        dplyr::select(entrez_id, ensembl_gene_id) %>% 
                                        dplyr::rename(ensembl_id = ensembl_gene_id)) %>% 
  dplyr::select(gene, entrez_id, ensembl_id, -ensembl.gene.id, everything()) %>% 
  group_by(subclass) %>% 
  arrange(subclass, -average.log.fold.change) %>% 
  ungroup()

# get ensembl list of markers
new_markers = sonny_hgnc_merged_markers %>% filter(used.in.mgp == "TRUE")
new_cell_types = new_markers %>% filter(!is.na(subclass)) %>% pull(subclass) %>% unique
new_marker_list  = lapply(new_cell_types, function(cell_type){
  cell_type_marker_list = new_markers %>% 
    filter(subclass == cell_type, ensembl_id %in% rownames(combined_matrix)) %>% 
    pull(ensembl_id)
  return(cell_type_marker_list)
})
names(new_marker_list) = new_cell_types
print(new_cell_types)
#write the new marker list to a csv file
# Find the maximum length of the lists
max_length <- max(sapply(new_marker_list, length))

# Pad each list with NA or empty strings to make them the same length
padded_list <- lapply(new_marker_list, function(x) {
  length(x) <- max_length
  return(x)
})

# Convert the padded list to a data frame
new_marker_df <- as.data.frame(padded_list)

# Write the data frame to a TSV file
#write.table(new_marker_df, file = "/nethome/kcni/xzhou/dan_paper/new_marker_list_symbol.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
################################################################################################################################################################################

# ###############################################################################################################################################################################
# #---------------------------------------------- Revision for Comment 1 from Reviewer 2 --------------------------------------------------
# ### Remove the eponymous marker genes (ie: SST gene for SST cell type)
# ### Can we still get the accurate cell type proportions without using the marker genes?
# #----------------------------------------------------------------------------------------------------------------------------------------
# # remove the SST gene from the marker list for SST cell type
# new_marker_list$SST = new_marker_list$SST[new_marker_list$SST != "ENSG00000157005"]

# # remove the PVALB gene from the marker list for PVALB cell type
# new_marker_list$PVALB = new_marker_list$PVALB[new_marker_list$PVALB != "ENSG00000100362"]

# # remove the VIP gene from the marker list for VIP cell type
# new_marker_list$VIP = new_marker_list$VIP[new_marker_list$VIP != "ENSG00000146469"]

################################################################################################################################################################################
### Update: Replace the marker genes for SST, PVALB, and VIP cell types with canonical markers
### Does the deconvolutionr esult change/improve when using the canonical markers?
#----------------------------------------------------------------------------------------------------------------------------------------
# # remove the SST gene from the marker list for SST cell type
# new_marker_list$SST = c("ENSG00000157005", "ENSG00000154645","ENSG00000006128","ENSG00000172137","ENSG00000122585","ENSG00000189056","ENSG00000181449","ENSG00000089250","ENSG00000136352")
# # remove the PVALB gene from the marker list for PVALB cell type
# new_marker_list$PVALB = c("ENSG00000100362","ENSG00000125740","ENSG00000178568","ENSG00000155657","ENSG00000144285")

# # remove the VIP gene from the marker list for VIP cell type
# new_marker_list$VIP = c("ENSG00000146469","ENSG00000172137","ENSG00000187094","ENSG00000147571","ENSG00000166736","ENSG00000185736","ENSG00000144724","ENSG00000106852")



################################################################################################################################################################################
###Normalize and process count matrices and perform MGP analysis
cell_types = new_marker_list$subclass %>% unique()

#Define vector of file names to perform MGP on and initialize empty variable before looping
matrix_names = c("libd_matrix",  "nimh_matrix", "gvex_matrix", "mssm_matrix", "penn_matrix", "pitt_matrix")
mgp_estimations = setNames(data.frame(matrix(ncol = length(c(colnames(combined_metadata), cell_types)), nrow = 0)), c(colnames(combined_metadata), cell_types))

for(matrix_name in matrix_names) {
  #Preprocessing
  matrix = get(matrix_name)
  cpm = cpm(matrix, log = TRUE, prior.count = 0.1)
  sds = rowSds(cpm, na.rm = TRUE)
  matrix = cpm[sds > 0.1,] %>% as.data.frame() %>% rownames_to_column(var = "gene_symbol") #Consider setting SD cutoff to 0 
  genes_only = matrix %>% subset(gene_symbol != "")
if(length(which(duplicated(genes_only$gene_symbol))) != 0){
  genes_only = genes_only[-which(duplicated(genes_only$gene_symbol)),]
}

estimations = mgpEstimate(
  exprData = genes_only,
  genes = new_marker_list,
  geneColName = 'gene_symbol',
  outlierSampleRemove = F, # should outlier samples removed. This is done using boxplot stats
  geneTransform = NULL, # this is the default option for geneTransform
  groups = NULL, # if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority = F)

#Coerce estimations list into data frame 
estimations_scaled <- estimations$estimates %>% as.data.frame() %>% scale() %>% as.data.frame() %>% tibble::rownames_to_column(var = "standardID")
estimations_metadata <- right_join(combined_metadata, estimations_scaled, by = "standardID")

#Merge cell type proportions with sample metadata
estimations_metadata = right_join(combined_metadata, estimations_scaled )

#Remove '+' from ageDeath for modelling
estimations_metadata$ageDeath = as.numeric(gsub("[+]", "", estimations_metadata$ageDeath))

mgp_estimations = rbind(mgp_estimations, estimations_metadata) 
}

#write.csv(mgp_estimations, "/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/psychencode_mgp_estimations.csv")

################################################################################################################################################################################
#load data from Dan's new MGP pipeline (separated the estimation for each brain bank institution)
#psy_est <- read.csv('/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/psychencode_mgp_estimations.csv') # dan changed the input file here, the org one should be /external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/psychencode_mgp_estimations.csv
psy_est <- mgp_estimations
#subset the data to only include the samples from MSSM, and from PITT
MSSM_est <- psy_est[psy_est$newStudy == 'MSSM', ] # 312 samples  55 features #NOT correct here!!@ Start from here!!!!!
PITT_est <- psy_est[psy_est$newStudy == 'Pitt', ] # 106 samples  55 features col

#load data from snCTP as ground truth
psy_snCTP <- read.csv('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/PsychEncode_label_transferred_snCTP.csv') # 140 samples 19 cell types
MSSM_snCTP <- psy_snCTP[psy_snCTP$Institution =='MtSinai',] # 92 samples; 20 cell types; 23 cols

#cleaning the data so that the sample names and cell types are consistent in both datasets with sames order
##samples
common_samples <- intersect(MSSM_est$specimenID, MSSM_snCTP$subject_id) #92 common samples -> 79 after remove low cell counts samples ie: <500
MSSM_est_common <- MSSM_est[MSSM_est$specimenID %in% common_samples,] #92 55
MSSM_snCTP_common <- MSSM_snCTP[MSSM_snCTP$subject_id %in% common_samples,] #92 24
rownames(MSSM_est_common) <- MSSM_est_common$specimenID
rownames(MSSM_snCTP_common) <- MSSM_snCTP_common$subject_id
MSSM_est_common <- MSSM_est_common[, (ncol(MSSM_est_common)-18):ncol(MSSM_est_common)]
MSSM_snCTP_common <- MSSM_snCTP_common[,3:21]
##cell types  
new_order <- c(1,2,3,4,5,6,7,8,9,10,11,13,12,14,16,15,17,18,19)
MSSM_est_common <- MSSM_est_common[,new_order]
if (all(colnames(MSSM_est_common) == colnames(MSSM_snCTP_common))) {
  colnames(MSSM_est_common) <- colnames(MSSM_snCTP_common)
} else {
  print("The column names of MSSM_est_common and MSSM_snCTP_common are not the same.")
}
# ##cell types
# new_order <- c(1,2,10,11,12,13,15,16,17,18)
# MSSM_est_common <- MSSM_est_common[,new_order]
# new_order <- c(1,2,4,5,6,7,8,9,10,11)
# MSSM_snCTP_common <- MSSM_snCTP_common[,new_order]
# colnames(MSSM_snCTP_common) <- colnames(MSSM_est_common)

#calculate the correlation between the estimated cell proportions and the snCTP cell proportions
predicted_df <- data.frame(Sample = rownames(MSSM_est_common), MSSM_est_common)
actual_df <- data.frame(Sample = rownames(MSSM_snCTP_common), MSSM_snCTP_common)

merged_df <- merge(predicted_df, actual_df, by = "Sample")

plots_list <- list()

for (cell_type in colnames(MSSM_est_common)) {
  # Calculate the prediction accuracy (e.g., correlation coefficient)
  x <- paste(cell_type, ".x", sep = "")
  y <- paste(cell_type, ".y", sep = "")
  accuracy <- cor(merged_df[[x]], merged_df[[y]], use = "pairwise.complete.obs")

  # Create the scatter plot with a trend line
  plot <- ggplot(merged_df, aes(x = .data[[y]], y = .data[[x]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
    labs(x = "Actual Proportions", y = "Predicted Proportions") +
    ggtitle(paste("Prediction Accuracy for", cell_type, ":", round(accuracy, 2)))

  # Store the scatter plot in the list
  plots_list[[cell_type]] <- plot
}

# Arrange and display the scatter plots in a grid
# png("/nethome/kcni/xzhou/label_transfer/benchmark_sep_cohorts_hodge_tax.png", width = 1600, height = 1000)
# grid.arrange(grobs = plots_list, ncol = 5)
# dev.off()

############################################# dtangle ratio for estimating cell type proportions for sep cohorts #######################################
#load ref without empty annotation
load('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/ref_matrix_all_annotated.RData') #ref_matrix
load('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/ref_anno.RData')
ref_anno_filtered <- ref_anno[ref_anno$subclass_label != '',]

#obtain updated mssm count matrix
matrix = get("mssm_matrix")
cpm = cpm(matrix, log = TRUE, prior.count = 0.1)
sds = rowSds(cpm, na.rm = TRUE)
matrix = cpm[sds > 0.1,] %>% as.data.frame() %>% rownames_to_column(var = "gene_symbol") #Consider setting SD cutoff to 0 
genes_only = matrix %>% subset(gene_symbol != "")
if(length(which(duplicated(genes_only$gene_symbol))) != 0){
genes_only = genes_only[-which(duplicated(genes_only$gene_symbol)),] # [1] 46417   313
}
MSSM_meta <- combined_metadata[combined_metadata$newStudy == 'MSSM',] #[1] 312  35

count_cpm <- genes_only
rowname <- genes_only$gene_symbol
rownames(count_cpm) <- rowname

#Converting Ensembl ID to Gene Name
#mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
ensembl = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id","description","gene_biotype","percentage_gene_gc_content"), mart = mart) 
ensembl_to_gene = (data.frame(ensembl$ensembl_gene_id, ensembl$hgnc_symbol))
names(ensembl_to_gene) = c("gene_symbol", "gene_name")
#remove duplicates
ensembl_to_gene = ensembl_to_gene[!duplicated(ensembl_to_gene[,1]),]
count_cpm = merge(x=count_cpm, y=ensembl_to_gene, by = "gene_symbol", all.x = T) # 46417 genes   314 samples for MSSM cohort
# some gene names are duplicated after matching (‘’, ‘DUXAP8’, ‘GOLGA8M’, ‘ITFG2-AS1’, ‘LINC01238’, ‘PINX1’, ‘POLR2J4’, ‘RN7SL274P’, ‘SIGLEC5’, ‘TUBB7P’), in order to set them as unqiue row names, we only keep the first occurence of the duplicated gene names
count_cpm_filtered <- count_cpm[!duplicated(count_cpm$gene_name), ] #[1] 31918 genes  314 samples ##################START HERE######################!!!!!!!!!!!!!!!!3.13
count_cpm_filtered <- count_cpm_filtered[count_cpm_filtered$gene_name != '', ] #remove the empty gene name
count_cpm_filtered <- count_cpm_filtered[!is.na(count_cpm_filtered$gene_name), ]
rownames(count_cpm_filtered) = count_cpm_filtered$gene_name
count_cpm_filtered <- count_cpm_filtered[ , !names(count_cpm_filtered) %in% c("gene_symbol", "gene_name")] #[1] 31916 genes   312 samples

# #find the genes that are common to Count and Reference. 
# NOTE: The following commented steps are IN ANOTHER FILE CALLED /nethome/kcni/xzhou/dan_paper/CMC_dtangle_y.r, since we need to submit a slurm job with larger memory for this.
# commongenes <- intersect(rownames(count_cpm_filtered), rownames(ref_matrix))
# count_final_common <- count_cpm_filtered[pmatch(commongenes, rownames(count_cpm_filtered)), ] #  23349   312
# ref_final_common <- ref_matrix[pmatch(commongenes, rownames(ref_matrix)), ] #23349 47432
# # saveRDS(ref_final_common, file="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/MC_dtangle_ref_final_common.RData")
# # saveRDS(count_final_common, file="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/CMC_dtangle_count_final_common.RData")

# #join the datasets
# y <- cbind(ref_final_common, count_final_common)
# #apply quantile normalization to these datasets in order to ensure that they are indeed comparable.
# y <- normalizeQuantiles(y)
# y <- t(y)

y <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/CMC_dtangle_y.RData")

all_cell_type <- unique(ref_anno_filtered$subclass_label)
pure_samples <- lapply(1:length(all_cell_type), function(i) {
    which(ref_anno_filtered$subclass_label == all_cell_type[i])
})
names(pure_samples) = all_cell_type

marker_ratio = find_markers(y,pure_samples=pure_samples,data_type="rna-seq",marker_method='ratio')

q = .1
quantiles = lapply(marker_ratio$V,function(x)quantile(x,1-q))
K = length(pure_samples)
n_markers = sapply(1:K,function(i){max(which(marker_ratio$V[[i]] > quantiles[[i]]))})

#run the deconvolution
marks = marker_ratio$L
dc_ratio <- dtangle(y, pure_samples=pure_samples, n_markers=n_markers, data_type = 'rna-seq', markers = marks)

final_est <- dc_ratio$estimates[(dim(ref_anno_filtered)[1]+1):dim(y)[1],]
colnames(final_est) <- all_cell_type
head(final_est) #312sampels 19 cell types

# # updates: Reviewer comment 1: for the new marker list, we need to remove the eponymous marker genes (ie: SST gene for SST cell type)
# #run the deconvolution
# marks = marker_ratio$L

# marks$SST <- marks$SST[names(marks$SST) != "SST"]
# marks$PVALB <- marks$PVALB[names(marks$PVALB) != "PVALB"]
# marks$VIP <- marks$VIP[names(marks$VIP) != "VIP"]
# dc_ratio_new <- dtangle(y, pure_samples=pure_samples, n_markers=n_markers, data_type = 'rna-seq', markers = marks)

# final_est_new <- dc_ratio_new$estimates[(dim(ref_anno_filtered)[1]+1):dim(y)[1],]
# colnames(final_est_new) <- all_cell_type
# head(final_est_new) #312sampels 19 cell types
################################################################################################################################################################################

#load the snRNA-seq data as ground truth for pearson correlation
# snRNA <- read.csv('/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/psychencode_snCTPs.csv') #140 samples
# snRNA_MSSM <- snRNA[snRNA$Institution =='MSSM',] 
# use this label transferred ref data: MSSM_snCTP


################################################### Find common cells using specimenID and bulk sample id as identifiers ########################################################

### For MGP
#find the common specimen id between the estimated cell type proportions and the bulk RNA-seq data
cmc_estimations_mgp <- MSSM_est_common
cmc_estimations_mgp_common <- cmc_estimations_mgp
common_specimen_id <- intersect(rownames(cmc_estimations_mgp_common), MSSM_snCTP$subject_id) #92 common samples
cmc_estimations_mgp_common <- cmc_estimations_mgp_common[rownames(cmc_estimations_mgp_common) %in% common_specimen_id,]
rownames(MSSM_snCTP) <- MSSM_snCTP$subject_id
MSSM_snCTP_common <- MSSM_snCTP[MSSM_snCTP$subject_id %in% common_specimen_id,c(3:21)]
colnames(cmc_estimations_mgp_common) == colnames(MSSM_snCTP_common) #92samples 19 cell types



### For dtangle
#Merge cell type proportions with specimenID
cmc_estimations_metadata <- read.csv("/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/cmc_cell_prop.csv")
cmc_estimations_metadata <- cmc_estimations_metadata[,colnames(cmc_estimations_metadata) %in% c("specimenID","individualID")]
final_est <- as.data.frame(final_est) # update here
final_est$individualID <- row.names(final_est)
cmc_estimations_dt <- left_join(final_est, cmc_estimations_metadata, by = c("individualID" = "individualID")) # 312 samples  55 features
nrow(cmc_estimations_dt[cmc_estimations_dt$specimenID %in% MSSM_snCTP$subject_id,]) # 92 samples
cmc_estimations_dt <- cmc_estimations_dt[cmc_estimations_dt$specimenID %in% MSSM_snCTP$subject_id,]
# clean up the cmc_estimations_mgp_common 
rownames(cmc_estimations_dt) <- cmc_estimations_dt$specimenID
cmc_estimations_dt_common <- cmc_estimations_dt[,c(1:(ncol(cmc_estimations_dt)-2))]

# re-order the columns in cmc_estimations_dt_common to match the cell types order in MSSM_snCTP_common
new_order <- c(6,17,3,18,15,7,8,10,13,2,14,5,11,4,16,12,9,1,19)
cmc_estimations_dt_common <- cmc_estimations_dt_common[,new_order]
colnames(cmc_estimations_dt_common) <- colnames(MSSM_snCTP_common) #92samples 19 celltypes

# Create function to calculate the spearman/pearson correlation coefficient and CI for two dataframes.
cor.test.plus <- function(x) {
  se <- sqrt((1 - x$estimate^2) / x$parameter)
  return(se)
}
cor_func <- function(estimations_df, actual_prop) {
  predicted_df <- data.frame(Sample = rownames(estimations_df), estimations_df)
  actual_df <- data.frame(Sample = rownames(actual_prop), actual_prop)

  merged_df <- merge(predicted_df, actual_df, by = "Sample")

  accuracy_list <- c()
  SE_list <- c()
  for (cell_type in colnames(estimations_df)) {
    x <- paste(cell_type, ".x", sep = "")
    y <- paste(cell_type, ".y", sep = "")
    test_result <- cor.test(merged_df[[x]], merged_df[[y]], method = "pearson")
    accuracy_list <- c(accuracy_list, test_result$estimate)
    SE_list <- c(SE_list, cor.test.plus(test_result))
  }
  return(list(accuracy_list, SE_list))
}

mgp_sn_cor <- cor_func(cmc_estimations_mgp_common, MSSM_snCTP_common)
dt_sn_cor <- cor_func(cmc_estimations_dt_common, MSSM_snCTP_common)

##################################################################################################################

# Create dataframe for plotting
df <- data.frame(
  Name = colnames(cmc_estimations_mgp_common),
  Category = "MGP",
  Accuracy = mgp_sn_cor[[1]],
  SE = mgp_sn_cor[[2]]
)

# Create dataframe for plotting
df_dt <- data.frame(
  Name = colnames(cmc_estimations_dt_common),
  Category = "dtangle",
  Accuracy = dt_sn_cor[[1]],
  SE = dt_sn_cor[[2]]
)

# Calculate lower and upper bounds for error bars
df$SE_Lower <- df$Accuracy - df$SE
df$SE_Upper <- df$Accuracy + df$SE

# Calculate lower and upper bounds for error bars
df_dt$SE_Lower <- df_dt$Accuracy - df_dt$SE
df_dt$SE_Upper <- df_dt$Accuracy + df_dt$SE
# ### The CI is too wide for some cell types, so we need to check if data falls into normal distribution, which is the assumption for Pearson correlation
# hist(cmc_estimations_mgp_common$SST, main = "Histogram", xlab = "Data", breaks = "Sturges") # SST is normally distributed
# shapiro.test(cmc_estimations_mgp_common$Oligodendrocyte) #but oligo does not. That's wgy the CI is too wide
# # Thus, We need to use pearson correlation instead of pearson correlation!

# Add a column for labeling the cell types to boarder cell classifications: excitaatory, inhibitory, and non-neuron
for (i in 1:length(df$Name)) {
  if (df$Name[i] %in% c("IT", "L4.IT", "L5.6.IT.Car3", "L5.ET", "L6.CT", "L6b", "L5.6.IT.Car3", "L5.6.NP")) { #Exc_L5.6.IT.Car3 -> L5.6.IT.Car3
    df$Label[i] <- "Excitatory"
  }
  else if (df$Name[i] %in% c("LAMP5", "PAX6", "PVALB", "SST", "VIP")) {
    df$Label[i] <- "Inhibitory"
  }
  else {
    df$Label[i] <- "Non-Neuronal"
}
}

## for dtangle
for (i in 1:length(df_dt$Name)) {
  if (df_dt$Name[i] %in% c("IT", "L4.IT", "L5.6.IT.Car3", "L5.ET", "L6.CT", "L6b", "L5.6.IT.Car3", "L5.6.NP")) {
    df_dt$Label[i] <- "Excitatory"
  }
  else if (df$Name[i] %in% c("LAMP5", "PAX6", "PVALB", "SST", "VIP")) {
    df_dt$Label[i] <- "Inhibitory"
  }
  else {
    df_dt$Label[i] <- "Non-Neuronal"
}
}

#Step 1: Reorder cell types within "Inhibitory" class
inhibitory_order <- c("PVALB", "SST", "VIP") # Specify the desired order for these cell types
df$Name <- factor(df$Name, levels = c(inhibitory_order, setdiff(df$Name, inhibitory_order)))
df_dt$Name <- factor(df_dt$Name, levels = c(inhibitory_order, setdiff(df_dt$Name, inhibitory_order)))
#Step 2: Reorder the classification in the facet grid
# Ensure "Inhibitory" class appears first
df$Label <- factor(df$Label, levels = c("Inhibitory", "Excitatory", "Non-Neuronal"))
df_dt$Label <- factor(df_dt$Label, levels = c("Inhibitory", "Excitatory", "Non-Neuronal"))
# # Use ggplot2 to create the bar plot
# color <- brewer.pal(length(categories), "Set3")
# p <- ggplot(df, aes(x = Name, y = Accuracy, fill = Category)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   scale_fill_manual(values = color) +
#   labs(title = "Pearson Correlation Coefficient for Cell Type Proportion Prediction by MGP",
#        x = "Cell Types",
#        y = "Pearson Correlation") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
#   facet_grid(~Label, drop = T, scale = "free", space = "free_x")

# print(p)

### Dan's code to make barplot
#saveRDS(df, file = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/correlation_barplot_danpaper_MGP.RData")
# Assuming df is your dataframe
theme_set(theme_classic2())
#Colour palette
cbPalette = c("#56B4E9", "#009E73","#E69F00", "#0072B2", "#D55E00", "#CC79A7","#000000","#F0E442")

figure_1b <- df %>% 
  ggplot(aes(x = Name, y = Accuracy)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = SE_Lower, ymax = SE_Upper), width = 0.2) +
  facet_grid(~Label, scales = "free", space = "free_x") +
  #geom_hline(yintercept = 0.2, color = "red", linewidth=1.5) + 
  ylab("Pearson's R") + 
  xlab('Cell Types') + 
  theme(
    panel.background = element_blank(),  # Keeps the panel background transparent
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),  # Removes minor grid lines
    strip.background = element_rect(fill = "white", colour = "black"),
    axis.line = element_line(color = "black"),  # Ensures x-axis and y-axis lines are visible
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 30),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(vjust = 8, size = 32),
    axis.title.y = element_text(size = 32),
    strip.text.x = element_text(size = 30)
  ) +
  guides(fill = FALSE) +
  scale_y_continuous(limits = c(-0.35, 0.7)) 
  #scale_fill_manual(values = "#A9A9A9") # Set the fill color for bars to grey

print(figure_1b)

# # #save this plot p into a RData file for future retrieval
# png("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/figure1_barplot_danpaper_MGP.png", width = 800, height = 600)
# print(figure_1b)
# dev.off()
# #saveRDS(p, file = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/correlation_barplot_danpaper_MGP.RData")
# new_p <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/correlation_barplot_danpaper_MGP.RData")

# print(new_p)


## Supplemenraty figure for dtangle barplot
#combine df and df_dt into a single merged_df and creat a barplot to show the pearson correlation of dt_sn and mgp_sn side by side
cbPalette = c("#56B4E9", "#009E73","#E69F00", "#0072B2", "#D55E00", "#CC79A7","#000000","#F0E442")
merged_df <- rbind(df, df_dt)  
merged_df$Category <- factor(merged_df$Category, levels = c("MGP", "dtangle")) 

supplementary_figure_1 <- merged_df %>% 
  ggplot(aes(x = Name, y = Accuracy, fill=Category)) + 
  geom_hline(yintercept = 0) + 
  #geom_hline(yintercept = 0.2, color = "red") +  # Add this line for the threshold
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7) + 
  geom_errorbar(aes(ymin = SE_Lower, ymax = SE_Upper), width = 0.7, position = position_dodge(width = 0.75)) +
  facet_grid(~Label, scales = "free", space = "free_x") +
  ylab("Pearson's R") + 
  xlab('Cell Types') + 
  theme(
    panel.background = element_blank(),  # Keeps the panel background transparent
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),  # Removes minor grid lines
    strip.background = element_rect(fill = "white", colour = "black"),
    axis.line = element_line(color = "black"),  # Ensures x-axis and y-axis lines are visible
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 20),
    axis.text.y = element_text(size = 28),
    axis.title.x = element_text(vjust = 8, size = 28),
    axis.title.y = element_text(size = 28),
    strip.text.x = element_text(size = 28),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  ) +
  guides(fill = guide_legend(title = "Deconvolution Tools")) +
  scale_fill_manual(values = c("MGP" = "#cc0000", "dtangle" = "#56B4E9")) 
  #scale_y_continuous(limits = c(-0.4, 0.6)) 
 

print(supplementary_figure_1)
#saveRDS(supplementary_figure_1, file="/nethome/kcni/xzhou/dan_paper/supplymentary_figure_1.RData")
supplementary_figure_1 <- readRDS("/nethome/kcni/xzhou/dan_paper/supplymentary_figure_1.RData")
#Save to plots folder
output_filename <- "/nethome/kcni/xzhou/dan_paper/figure_S1.png"

# Use ggsave to save the plot as a png image
# ggsave(output_filename, supplementary_figure_1, width = 15, height = 8, dpi = 300)

figure_S1 = supplementary_figure_1 / figure_1a / figure_1a_dt
figure_S1 = figure_S1 + plot_annotation(tag_levels = 'A') + plot_layout(heights = c(2,1,1)) & theme(plot.tag = element_text(size = 25))
#figure_1

### Update the supplymentary figure 1 with scatterplots of eponymous marker genes removed
#Save to plots folder
output_filename <- "/nethome/kcni/xzhou/dan_paper/figure_S1_add_no_eponymous_marker.png"

# Use ggsave to save the plot as a JPG image
#ggsave(output_filename, figure_S1, width = 22, height = 30, dpi = 300)



########################## Scatter plot ###############################
# plots_list <- list()

# predicted_df <- data.frame(Sample = rownames(cmc_estimations_mgp_common), cmc_estimations_mgp_common[,16:18])
# actual_df <- data.frame(Sample = rownames(MSSM_snCTP_common), MSSM_snCTP_common[,16:18])
# merged_df <- merge(predicted_df, actual_df, by = "Sample")

# ### without case/control info
# for (cell_type in colnames(cmc_estimations_mgp_common[16:18])) {
#   # Calculate the prediction accuracy (e.g., correlation coefficient)
#   x <- paste(cell_type, ".x", sep = "")
#   y <- paste(cell_type, ".y", sep = "")
#   accuracy <- cor(merged_df[[x]], merged_df[[y]], method="pearson", use = "pairwise.complete.obs")

#   # Create the scatter plot with a trend line, swapping x and y aesthetics
#   plot <- ggplot(merged_df, aes(x = .data[[y]], y = .data[[x]])) +  # Add color mapping
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
#     labs(x = "snCTP", y = "rCTP") +  # Swap axis labels
#     ggtitle(paste("Pearson Correlation", cell_type, ":", round(accuracy, 2)))

#   # Store the scatter plot in the list
#   plots_list[[cell_type]] <- plot
# }

### Dan's color and settings for ggolot
#combine the cmc_estimations_mgp_common and MSSM_snCTP_common into single dataframe

### RMEMEBER TO UPDATE THE MARKER GENE LIST!!!!!!! Then you can proceed to the next step:
predicted_df <- data.frame(Sample = rownames(cmc_estimations_mgp_common), cmc_estimations_mgp_common[,16:18])
predicted_df_long <- pivot_longer(predicted_df, cols = c(2:4), names_to = "cell_type", values_to = "rCTP")
actual_df <- data.frame(Sample = rownames(MSSM_snCTP_common), MSSM_snCTP_common[,16:18])
actual_df_long <- pivot_longer(actual_df, cols = c(2:4), names_to = "cell_type", values_to = "snCTP")

#merge the two dataframes into a single one using sampls and cell_type as keys
merged_df <- merge(predicted_df_long, actual_df_long, by = c("Sample", "cell_type")) #276 for 3 cell types, so 92 samples for each cell type

# Step 1: Calculate pearson r for each cell type
correlation_df <- merged_df %>%
  group_by(cell_type) %>%
  summarize(cor_value = cor(snCTP, rCTP, method = "pearson", use="pairwise.complete.obs")) %>%
  ungroup()



# saveRDS(merged_df, file = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/Fig1A_scatterplot_merged_df.RData")

# Your plotting code, now including geom_text for the correlation coefficient
figure_1a <- merged_df %>%
  ggplot(aes(x = snCTP*100, y = rCTP)) + 
  geom_smooth(se = FALSE, method = 'lm', fullrange = TRUE) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~")), 
           color = "black", 
           geom = "label", 
           size = 10, 
           label.y = max(merged_df$rCTP) + 1) + # Adjust label position
  geom_point(alpha = 1, size = 3, color = "black") +
  ylab('rCTP (AU)') + 
  xlab('snCTP (%)') +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(vjust = -3, size = 32),
    plot.margin = margin(t = 2, r = 1, b = 1, l = 1, unit = "cm"), # Match margin
    axis.title.y = element_text(size = 32),
    legend.position = c(0.9, 0.9),
    legend.title = element_blank(),
    legend.text = element_text(size = 32),
    strip.text.x = element_text(size = 32)
  ) +
  facet_grid(~cell_type, scale = "free_x") +
  theme(panel.spacing = unit(1, "cm", data = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) # Match y expansion

print(figure_1a)


#saveRDS(figure_1a, file="/nethome/kcni/xzhou/dan_paper/figure_1a_revised_no_eponymous_marker.RData")
#figure_1a <- readRDS("/nethome/kcni/xzhou/dan_paper/figure_1a_revised_no_eponymous_marker.RData")
#ggsave(figure_1a, file = "/nethome/kcni/xzhou/dan_paper/figure_1a_canonical_marker.png", width = 20, height = 15, dpi = 300)

figure_1 = figure_1a / figure_1b
figure_1 = figure_1 + plot_annotation(tag_levels = 'A') + plot_layout(heights = c(2,1)) & theme(plot.tag = element_text(size = 25))
#figure_1

#Save to plots folder
output_filename <- "/nethome/kcni/xzhou/dan_paper/figure_1_july9.png"

# Use ggsave to save the plot as a JPG image
# ggsave(output_filename, figure_1, width = 20, height = 15, dpi = 300)

### Update: March, 2025
# Create a spacer plot for Figure A
spacer <- plot_spacer() + theme_void()

# Arrange the plots
figure_1 = spacer / figure_1a / figure_1b +
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(1, 2, 1)) #& # Adjust heights to allocate space
  #theme(plot.tag = element_text(size = 25))


figure_1 = figure_1a / figure_1b
figure_1 = figure_1 + plot_annotation(tag_levels = 'A') + plot_layout(heights = c(2,1)) & theme(plot.tag = element_text(size = 25))
# Save the plot
output_filename <- "/nethome/kcni/xzhou/dan_paper/figure_1_Mar_2025.svg"
#ggsave(output_filename, figure_1, width = 20, height = 20, dpi = 300, device = "svg")

print(figure_1)

output_filename <- "/nethome/kcni/xzhou/dan_paper/figure_1_Mar_2025.png"
#ggsave(output_filename, figure_1, width = 20, height = 15, dpi = 300)

# #------------------------------------------------ Plot the RMSE for each cell type------------------------------------------------
# library(ggplot2)
# library(dplyr)
# library(tidyr)

# # Calculate RMSE for each cell type
# rmse_df <- merged_df %>%
#   group_by(cell_type) %>%
#   summarize(rmse_value = sqrt(mean((snCTP - rCTP)^2, na.rm = TRUE))) %>%
#   ungroup()

# # Plot with RMSE
# figure_1a <- merged_df %>%
#   ggplot(aes(x = snCTP * 100, y = rCTP)) +
#   geom_smooth(se = FALSE, method = 'lm', fullrange = TRUE) +
#   geom_point(alpha = 1, size = 3, color = "black") +
#   geom_text(
#     data = rmse_df, 
#     aes(label = sprintf("RMSE = %.2f", rmse_value)),
#     x = 0, y = Inf, # Adjust y position as needed
#     hjust = 0, vjust = 1.5, 
#     size = 7, color = "blue",
#     inherit.aes = FALSE
#   ) +
#   ylab('rCTP (AU)') + 
#   xlab('snCTP (%)') +
#   theme(
#     panel.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     strip.background = element_rect(fill = "white", colour = "black"),
#     axis.line = element_line(color = "black"),
#     axis.text.x = element_text(size = 32),
#     axis.text.y = element_text(size = 32),
#     axis.title.x = element_text(vjust = -3, size = 32),
#     plot.margin = margin(1, 1, 1, 1, "cm"),
#     axis.title.y = element_text(size = 32),
#     legend.position = c(0.9, 0.9),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 32),
#     strip.text.x = element_text(size = 32)
#   ) +
#   facet_grid(~cell_type, scale = "free_x") +
#   theme(panel.spacing = unit(1, "cm"))

# print(figure_1a)
# ggsave(figure_1a, file = "/nethome/kcni/xzhou/dan_paper/figure_1a_rmse.png", width = 20, height = 15, dpi = 300)
# #------------------------------------------------------------------------------------------------------------


##########################################################################################################################################################
#------------------------------------------------------------------------------------------------------------
# 2025-03-15 Update: scatter plot for dtangle estimtaion aftr removing the eponymous marker genes -  REMEBER TO  UPDATE THE INPUT FILE at #update here !!!!!

### REMEBER TO  UPDATE THE INPUT FILE at #update here !!!!! Then we can proceed to the next step:

#combine the cmc_estimations_mgp_common and MSSM_snCTP_common into single dataframe
predicted_df <- data.frame(Sample = rownames(cmc_estimations_dt_common), cmc_estimations_dt_common[,16:18])
predicted_df_long <- pivot_longer(predicted_df, cols = c(2:4), names_to = "cell_type", values_to = "rCTP")
actual_df <- data.frame(Sample = rownames(MSSM_snCTP_common), MSSM_snCTP_common[,16:18])
actual_df_long <- pivot_longer(actual_df, cols = c(2:4), names_to = "cell_type", values_to = "snCTP")

#merge the two dataframes into a single one using sampls and cell_type as keys
merged_df <- merge(predicted_df_long, actual_df_long, by = c("Sample", "cell_type")) #276 for 3 cell types, so 92 samples for each cell type

# Step 1: Calculate pearson r for each cell type
correlation_df <- merged_df %>%
  group_by(cell_type) %>%
  summarize(cor_value = cor(snCTP, rCTP, method = "pearson", use="pairwise.complete.obs")) %>%
  ungroup()



# saveRDS(merged_df, file = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/Fig1A_scatterplot_merged_df.RData")

figure_1a_dt <- merged_df %>%
  ggplot(aes(x = snCTP*100, y = rCTP*100)) +
  geom_smooth(se = FALSE, method = 'lm', fullrange = TRUE) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~")), 
           color = "black", 
           geom = "label", 
           size = 10,
           label.y = max(merged_df$rCTP*100) + 1) + # Adjust this value to move labels up/down
  geom_point(alpha = 1, size = 3, color="black") +
  ylab('dtangle_absCTP (%)') + 
  xlab('snCTP (%)') +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(vjust = -3, size = 32),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"), # Increased top margin
    axis.title.y = element_text(size = 32),
    legend.position = c(0.9, 0.9),
    legend.title = element_blank(),
    legend.text = element_text(size = 32),
    strip.text.x = element_text(size = 32),
  ) +
  facet_grid(~cell_type, scale="free_x") +
  theme(panel.spacing = unit(1, "cm", data = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) # Adds 20% expansion to top of plot

print(figure_1a_dt)
#saveRDS(figure_1a_dt, file="/nethome/kcni/xzhou/dan_paper/figure_1a_dt_no_eponymous_marker.RData")
#figure_1a_dt <- readRDS("/nethome/kcni/xzhou/dan_paper/figure_1a_dt_no_eponymous_marker.RData")
#ggsave(figure_1a_dt, file = "/nethome/kcni/xzhou/dan_paper/figure_1a_dt.png", width = 20, height = 15, dpi = 300)


#------------------------------------------------ Plot the scatter plot with full marker genes ------------------------------------------------

#combine the cmc_estimations_mgp_common and MSSM_snCTP_common into single dataframe
predicted_df <- data.frame(Sample = rownames(cmc_estimations_dt_common), cmc_estimations_dt_common[,16:18])
predicted_df_long <- pivot_longer(predicted_df, cols = c(2:4), names_to = "cell_type", values_to = "rCTP")
actual_df <- data.frame(Sample = rownames(MSSM_snCTP_common), MSSM_snCTP_common[,16:18])
actual_df_long <- pivot_longer(actual_df, cols = c(2:4), names_to = "cell_type", values_to = "snCTP")

#merge the two dataframes into a single one using sampls and cell_type as keys
merged_df <- merge(predicted_df_long, actual_df_long, by = c("Sample", "cell_type")) #276 for 3 cell types, so 92 samples for each cell type

# Step 1: Calculate pearson r for each cell type
correlation_df <- merged_df %>%
  group_by(cell_type) %>%
  summarize(cor_value = cor(snCTP, rCTP, method = "pearson", use="pairwise.complete.obs")) %>%
  ungroup()


figure_1a_dt_full_marker <- merged_df %>%
  ggplot(aes(x = snCTP*100, y = rCTP*100)) +
  geom_smooth(se = FALSE, method = 'lm', fullrange = TRUE) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~")), 
           color = "black", 
           geom = "label", 
           size = 10,
           label.y = max(merged_df$rCTP*100) + 1) + # Adjust this value to move labels up/down
  geom_point(alpha = 1, size = 3, color="black") +
  ylab('absCTP (%)') + 
  xlab('snCTP (%)') +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(vjust = -3, size = 32),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"), # Increased top margin
    axis.title.y = element_text(size = 32),
    legend.position = c(0.9, 0.9),
    legend.title = element_blank(),
    legend.text = element_text(size = 32),
    strip.text.x = element_text(size = 32),
  ) +
  facet_grid(~cell_type, scale="free_x") +
  theme(panel.spacing = unit(1, "cm", data = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) # Adds 20% expansion to top of plot

print(figure_1a_dt_full_marker)
#ggsave(figure_1a_dt_full_marker, file = "/nethome/kcni/xzhou/dan_paper/figure_1a_dt_full_marker.png", width = 20, height = 15, dpi = 300)


##########################################################################################################################################################



##########################################################################################################################################################
########################################################## Test if sex and age variabilities influcing correlation ##########################################################
### get metadata that inlcudes sex and age with sample names
MSSM_est_with_meta <- MSSM_est[MSSM_est$specimenID %in% common_samples,] #92 55
rownames(MSSM_est_with_meta) <- MSSM_est_with_meta$specimenID
MSSM_est_with_meta <- data.frame(
  reportedGender = MSSM_est_with_meta$reportedGender,
  ageDeath = MSSM_est_with_meta$ageDeath,
  sampleID = MSSM_est_with_meta$specimenID
)

#Merge cell type proportions with sample metadata
MSSM_est_metadata <- right_join(MSSM_est_with_meta, predicted_df, by = c("sampleID" = "Sample"))

# create a new column for msex, 0 represents female and 1 represents male
MSSM_est_metadata$msex <- ifelse(MSSM_est_metadata$reportedGender == "female", 0, 1)

#create a new column for age_death_sex
MSSM_est_metadata$age_death_sex <- MSSM_est_metadata$ageDeath * MSSM_est_metadata$msex
#create a new column for age_death2
MSSM_est_metadata$age_death2 <- MSSM_est_metadata$ageDeath*MSSM_est_metadata$ageDeath
#create a new column for age_death2_sex
MSSM_est_metadata$age_death2_sex <- MSSM_est_metadata$ageDeath*MSSM_est_metadata$ageDeath*MSSM_est_metadata$msex

## perform linear regression to regress out the sex and age, and then apply rank-based inverse normal transformation
cell_types = c("PVALB", "SST", "VIP")
pheno_lms = lapply(cell_types, function(cell_type){
  lm = paste0("scale(", cell_type, ")", " ~ msex + scale(ageDeath) + scale(age_death_sex) + scale(age_death2) + scale(age_death2_sex)") 
  results = lm(lm, data = MSSM_est_metadata) %>% 
  residuals()
  results <- results[!is.na(results)]  %>% #the last two columns are always NA???
 # RNOmni::RankNorm() %>%  #RINT here: apply rank-based inverse normal transformation
  as.data.frame()
  results$cell_type = cell_type
  results$FID = MSSM_est_metadata$sampleID
  return(results)
}) %>% bind_rows()


#format pheno_lms into phenotype dataframe
names(pheno_lms)[1] <- "transformed_residuals"

# Pivot the dataframe
pivot_df <- pivot_wider(data = pheno_lms, 
                        names_from = "cell_type", 
                        values_from = "transformed_residuals")


### same code for snRNAseq
#Merge cell type proportions with sample metadata
MSSM_sn_metadata <- right_join(MSSM_est_with_meta, actual_df, by = c("sampleID" = "Sample"))

# create a new column for msex, 0 represents female and 1 represents male
MSSM_sn_metadata$msex <- ifelse(MSSM_sn_metadata$reportedGender == "female", 0, 1)

#create a new column for age_death_sex
MSSM_sn_metadata$age_death_sex <- MSSM_sn_metadata$ageDeath * MSSM_sn_metadata$msex
#create a new column for age_death2
MSSM_sn_metadata$age_death2 <- MSSM_sn_metadata$ageDeath*MSSM_sn_metadata$ageDeath
#create a new column for age_death2_sex
MSSM_sn_metadata$age_death2_sex <- MSSM_sn_metadata$ageDeath*MSSM_sn_metadata$ageDeath*MSSM_sn_metadata$msex

## perform linear regression to regress out the sex and age, and then apply rank-based inverse normal transformation
cell_types = c("PVALB", "SST", "VIP")
pheno_lms_sn = lapply(cell_types, function(cell_type){
  lm = paste0("scale(", cell_type, ")", " ~ msex + scale(ageDeath) + scale(age_death_sex) + scale(age_death2) + scale(age_death2_sex)") 
  results = lm(lm, data = MSSM_sn_metadata) %>% 
  residuals()
  results <- results[!is.na(results)]  %>% #the last two columns are always NA???
 # RNOmni::RankNorm() %>%  #RINT here: apply rank-based inverse normal transformation
  as.data.frame()
  results$cell_type = cell_type
  results$FID = MSSM_sn_metadata$sampleID
  return(results)
}) %>% bind_rows()


#format pheno_lms into phenotype dataframe
names(pheno_lms_sn)[1] <- "transformed_residuals" 

# Pivot the dataframe
pivot_df_sn <- pivot_wider(data = pheno_lms_sn, 
                        names_from = "cell_type", 
                        values_from = "transformed_residuals") 



predicted_df_long <- pivot_longer(pivot_df, cols = c(2:4), names_to = "cell_type", values_to = "rCTP")
actual_df_long <- pivot_longer(pivot_df_sn, cols = c(2:4), names_to = "cell_type", values_to = "snCTP")

#merge the two dataframes into a single one using sampls and cell_type as keys
merged_df <- merge(predicted_df_long, actual_df_long, by = c("FID", "cell_type")) #276 for 3 cell types, so 92 samples for each cell type

# Step 1: Calculate pearson r for each cell type
correlation_df <- merged_df %>%
  group_by(cell_type) %>%
  summarize(cor_value = cor(snCTP, rCTP, method = "pearson", use="pairwise.complete.obs")) %>%
  ungroup()


# Your plotting code, now including geom_text for the correlation coefficient
figure_1a_reg_sex_age <- merged_df %>%
  ggplot(aes(x = snCTP*100, y = rCTP)) +
  geom_smooth(se = FALSE, method = 'lm', fullrange = TRUE) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~")), color = "black", geom = "label", size=10, label.y=3)+
  geom_point(alpha = 1, size = 3, color="black") +
  # geom_text(
  #   data = correlation_df, 
  #   aes(label = sprintf("Pearson's R = %.2f", cor_value)),
  #   x = 0, y = 1, # You'll need to adjust these values based on your actual plot range
  #   hjust = 0, vjust = 1, 
  #   size = 5, color = "blue",
  #   inherit.aes = FALSE # This is to avoid inheriting the other aesthetics like fill or group
  # ) +
  ylab('rCTP (AU)') + 
  xlab('snCTP (%)') +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(vjust = -3,size = 32),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title.y = element_text(size = 32),
    legend.position = c(0.9, 0.9),
    legend.title = element_blank(),
    legend.text = element_text(size = 32),
    strip.text.x = element_text(size = 32)
  ) +
  facet_grid(~cell_type, scale="free_x") +
  theme(panel.spacing = unit(1, "cm", data = NULL))
  # +
  #scale_y_continuous(limits = c(-4, 4))

print(figure_1a_reg_sex_age)
#ggsave(figure_1a_reg_sex_age, file = "/nethome/kcni/xzhou/dan_paper/figure_1a_reg_sex_age.png", width = 20, height = 15, dpi = 300)
##########################################################################################################################################################




##########################################################################################################################################################
### 2025-03-04 Update: Addressing Cross-Specificity
# Correlation Analysis: Perform a correlation analysis between the estimated proportions of all pairs of interneuron types (e.g., SST MGP vs. PVALB snCTPs, SST MGP vs. VIP snCTPs, etc.). 
# This will help determine if there are any unexpected correlations between different cell types that could indicate cross-specificity.
# Assuming merged_df contains columns for each cell type's rCTP and snCTP
predicted_df <- data.frame(Sample = rownames(cmc_estimations_mgp_common), cmc_estimations_mgp_common[,16:18])
predicted_df_long <- pivot_longer(predicted_df, cols = c(2:4), names_to = "cell_type", values_to = "rCTP")
actual_df <- data.frame(Sample = rownames(MSSM_snCTP_common), MSSM_snCTP_common[,16:18])
actual_df_long <- pivot_longer(actual_df, cols = c(2:4), names_to = "cell_type", values_to = "snCTP")

#merge the two dataframes into a single one using sampls and cell_type as keys
merged_df <- merge(predicted_df_long, actual_df_long, by = c("Sample", "cell_type")) #276 for 3 cell types, so 92 samples for each cell type

# Step 1: Calculate pearson r for each cell type
correlation_df <- merged_df %>%
  group_by(cell_type) %>%
  summarize(cor_value = cor(snCTP, rCTP, method = "pearson", use="pairwise.complete.obs")) %>%
  ungroup()
  
# Example cell types
cell_types <- c("SST", "PVALB", "VIP")  # Add more cell types as needed

# Initialize a matrix to store correlation results
correlation_matrix <- matrix(NA, nrow = length(cell_types), ncol = length(cell_types))
rownames(correlation_matrix) <- cell_types
colnames(correlation_matrix) <- cell_types

# Calculate correlations between all pairs of cell types, including self-comparisons
for (i in seq_along(cell_types)) {
  for (j in seq_along(cell_types)) {
    # Filter data for each cell type
    data_i <- merged_df %>% filter(cell_type == cell_types[i])
    data_j <- merged_df %>% filter(cell_type == cell_types[j])
    
    # Ensure both data frames have the same samples
    common_samples <- intersect(data_i$Sample, data_j$Sample)
    data_i <- data_i %>% filter(Sample %in% common_samples)
    data_j <- data_j %>% filter(Sample %in% common_samples)
    
    # Calculate correlation between rCTP of cell type i and snCTP of cell type j
    correlation_matrix[i, j] <- cor(data_j$snCTP, data_i$rCTP, method = "pearson")
  }
}

# Print the correlation matrix
print(correlation_matrix)

# Melt the correlation matrix for ggplot2
correlation_melted <- melt(correlation_matrix, na.rm = TRUE)

# Load the viridis package
library(viridis)

# Assuming correlation_melted is your melted correlation matrix
# Assuming correlation_melted is your melted correlation matrix
heatmap_plot <- ggplot(correlation_melted, aes(x=Var1, y=Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis(
    option = "viridis",
    name = "Pearson\nCorrelation R",
    limits = c(0, 0.6),
    breaks = seq(0, 0.6, by = 0.1),
    oob = scales::squish
  ) +
  theme_minimal(base_size = 35) +
  theme(
    axis.text.x = element_text(size = 35),
    axis.text.y = element_text(size = 35),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    plot.title = element_text(size = 42, hjust = 0.5),
    legend.title = element_text(size = 36),
    legend.text = element_text(size = 32),
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),
    legend.key.width = unit(1, "cm"),  # Increase legend key width
    legend.key.height = unit(1.5, "cm"),  # Increase legend key height
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = "snCTP (%)", y = "rCTP (AU)", title = "Cross-Specificity Correlation Heatmap for Major Interneurons")

ggsave(heatmap_plot, file = "/nethome/kcni/xzhou/dan_paper/figure_1a_comment3.png", width = 20, height = 15, dpi = 300)


### normalize the diagnal to allow the off diagnal to be more prominent
# Assuming correlation_melted is your melted correlation matrix
# Create a mask for the diagonal
correlation_melted$mask <- ifelse(correlation_melted$Var1 == correlation_melted$Var2, NA, correlation_melted$value)

heatmap_plot <- ggplot(correlation_melted, aes(x=Var1, y=Var2, fill = mask)) +
  geom_tile() +
  scale_fill_viridis(
    option = "viridis",
    name = "Pearson\nCorrelation",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    oob = scales::squish
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    plot.title = element_text(size = 28, hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = "Cell Type (snCTP)", y = "Cell Type (rCTP)", title = "Cross-Specificity Correlation Heatmap for Major Interneurons")


### Draw the correlation for all snRNAseq cell types
# Assuming merged_df contains columns for each cell type's snCTP
# Example cell types
cell_types <- c("SST", "PVALB", "VIP")  # Add more cell types as needed

# Initialize a matrix to store correlation results
correlation_matrix <- matrix(NA, nrow = length(cell_types), ncol = length(cell_types))
rownames(correlation_matrix) <- cell_types
colnames(correlation_matrix) <- cell_types

# Calculate correlations between all pairs of snCTP cell types
for (i in seq_along(cell_types)) {
  for (j in seq_along(cell_types)) {
    # Filter data for each cell type
    data_i <- merged_df %>% filter(cell_type == cell_types[i])
    data_j <- merged_df %>% filter(cell_type == cell_types[j])
    
    # Ensure both data frames have the same samples
    common_samples <- intersect(data_i$Sample, data_j$Sample)
    data_i <- data_i %>% filter(Sample %in% common_samples)
    data_j <- data_j %>% filter(Sample %in% common_samples)
    
    # Calculate correlation between snCTP of cell type i and snCTP of cell type j
    correlation_matrix[i, j] <- cor(data_i$snCTP, data_j$snCTP, method = "pearson")
  }
}

# Print the correlation matrix
print(correlation_matrix)

# Melt the correlation matrix for ggplot2
correlation_melted <- melt(correlation_matrix, na.rm = TRUE)

# Plot the heatmap
heatmap_plot <- ggplot(correlation_melted, aes(x=Var1, y=Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("lightcoral", "salmon", "tomato", "red", "darkred"),
    values = scales::rescale(c(0.2, 0.4, 0.6, 0.8, 1.0)),
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal(base_size = 20) +  # Increase base text size
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 24),  # Increase x-axis title size
    axis.title.y = element_text(size = 24),  # Increase y-axis title size
    plot.title = element_text(size = 28, hjust = 0.5),  # Increase plot title size and center it
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background to white
    plot.background = element_rect(fill = "white", color = NA)  # Set plot background to white
  ) +
  labs(x = "Cell Type (snCTP)", y = "Cell Type (snCTP)", title = "Correlation Heatmap for snRNAseq Cell Types")

#ggsave(heatmap_plot, file = "/nethome/kcni/xzhou/dan_paper/snCTP_correlation_heatmap.png", width = 20, height = 15, dpi = 300)



































################################################################################################################################################################################
# Whether tech cov regression and batch correction improve the prediction accuracy?
## Try to regress out the library batch and sequencing batch for quick assessment

#Preprocessing
matrix = get("mssm_matrix")
cpm = cpm(matrix, log = TRUE, prior.count = 0.1)
sds = rowSds(cpm, na.rm = TRUE)
matrix = cpm[sds > 0.1,] %>% as.data.frame() %>% rownames_to_column(var = "gene_symbol") #Consider setting SD cutoff to 0 
genes_only = matrix %>% subset(gene_symbol != "")
if(length(which(duplicated(genes_only$gene_symbol))) != 0){
genes_only = genes_only[-which(duplicated(genes_only$gene_symbol)),] # [1] 46417   313
}
MSSM_meta <- combined_metadata[combined_metadata$newStudy == 'MSSM',] #[1] 312  35

genes_only_new <- genes_only
rowname <- genes_only$gene_symbol
rownames(genes_only_new) <- rowname
#remove the first column gene_symbol
genes_only_new <- genes_only_new[,-1]


################## PHENOTYPE MANIPULATION ##################
#match the row order of the MSSM_meta by the same order of the the colnames of genes_only_new
MSSM_meta <- MSSM_meta[match(colnames(genes_only_new),MSSM_meta$individualID),]
CMC_rnaseq_meta <- read.csv("/external/rprshnas01/external_data/CommonMind/ControlledAccess/Data/RNAseq/Release4/Metadata/CMC_Human_rnaSeq_metadata.csv")

#join two metadata by individual ID
MSSM_meta_combined <- merge(MSSM_meta, CMC_rnaseq_meta, by.x = "individualID", by.y = "Individual.ID") #308 samples 73col



#Prepare metadata for PCATools::pca()
#remove tech cov with NA
MSSM_meta_combined_noNA <- MSSM_meta_combined[complete.cases(t(MSSM_meta_combined))] #[1]  308  65
#remove tech cov with only one value 
# Step 1: Calculate the number of unique values for each column
num_unique_values <- sapply(MSSM_meta_combined_noNA , function(x) length(unique(x)))
# Step 2: Create a logical vector to identify columns with more than one unique value
cols_to_keep <- num_unique_values > 1
# Step 3: Subset the DataFrame to include only those columns
MSSM_meta_combined_noNA_rmVar <- MSSM_meta_combined_noNA[, cols_to_keep] #[1] 308  38
rownames(MSSM_meta_combined_noNA_rmVar) <- MSSM_meta_combined_noNA_rmVar$individualID

tech_covar_list <- c("PMI","pH","RIN.x","Total_RNA_Yield", "Flowcell_.Given_to_Core", 
                    "Total_Reads", "Mapped_Reads", "Genes_Detected", "Transcripts_Detected", 
                    "Percent_Aligned", "Intragenic_Rate", "Intronic_Rate", "Intergenic_Rate", "Expression_Profiling_Efficiency", "rRNA_Rate") #15
#need to make sure the tech cov matrix contains numerical values for every column
batch_list <- c("RNA_Isolation_Batch","RNA_Prep_Date","Library_Batch","Flowcell_Batch","Ribozero_Batch")
metadata_nobatch <- MSSM_meta_combined_noNA_rmVar[,colnames(MSSM_meta_combined_noNA_rmVar) %in% tech_covar_list] # 308  13
# Check if every column contains all numeric values
non_numeric_columns <- names(which(!sapply(metadata_nobatch, function(x) all(is.numeric(x)))))
non_numeric_columns

#make sure all columns numeric now
for (i in seq_len(length(non_numeric_columns))) {
  metadata_nobatch[, non_numeric_columns[i]] <- as.numeric(as.factor(metadata_nobatch[, non_numeric_columns[i]]))
}


# prepare certain phenotypes: facotr or numeric values
metadata <- metadata_nobatch
for (tech_cov in colnames(metadata)){
  if(is.numeric(metadata[[tech_cov]])){
    metadata[[tech_cov]] <- as.numeric(metadata[[tech_cov]])}
  else{
    metadata[[tech_cov]] <- factor(metadata[[tech_cov]])
    print(tech_cov)
  }
}

# double check to remove samples from the pdata that have any NA value
discard <- apply(metadata, 1, function(x) any(is.na(x)))
metadata <- metadata[!discard,] #308  13 nothing changed

# filter the expression data (log tansfromed and var removed) to match the samples in our pdata
zscore_data <- genes_only_new 
zscore_data <- zscore_data[,which(colnames(zscore_data) %in% rownames(metadata))] #[1] 46417   308; 312 samplses -> 308 samples

# check that sample names match exactly between pdata and expression data 
all(colnames(zscore_data) == rownames(metadata)) #TRUE

# try PCAtools 
p <- PCAtools::pca(zscore_data, metadata = metadata)

# the function to calculate the correlation between the PCs and the tech covariates
eigencorvalue <- function(pcaobj, components = getComponents(pcaobj, seq_len(10)), metavars,  corUSE = 'pairwise.complete.obs', corFUN = 'pearson') {
  data <- pcaobj$rotated
  metadata <- pcaobj$metadata
  corvals <- list()

  # Issue warning if any columns to use are not numeric
  for (i in seq_len(length(components))) {
    if (!is.numeric(data[, components[i]])) {
      warning(components[i], ' is not numeric - please check the source data')
    }
  }
  for (i in seq_len(length(metavars))) {
    if (!is.numeric(metadata[, metavars[i]])) {
      warning(metavars[i], ' is not numeric - please check the source data')
    }
  }

  # Convert the data for x and y to data matrix
  xvals <- data.matrix(data[, which(colnames(data) %in% components), drop = FALSE])
  yvals <- metadata[, which(colnames(metadata) %in% metavars), drop = FALSE]

  # Ensure that non-numeric variables are coerced to numeric
  chararcter_columns <- unlist(lapply(yvals, is.numeric))
  chararcter_columns <- !chararcter_columns
  chararcter_columns <- names(which(chararcter_columns))
  for (c in chararcter_columns) {
    yvals[, c] = as.numeric(as.factor(yvals[, c]))
  }

  yvals <- data.matrix(yvals)

  # Create correlation table
  cor_matrix <- cor(xvals, yvals, use = corUSE, method = corFUN)

  # Store the correlation values and their corresponding variable names in a list
  corvals$x_names <- colnames(xvals)
  corvals$y_names <- colnames(yvals)
  corvals$correlation_values <- cor_matrix

  # Return the list containing the correlation values and variable names
  return(corvals)
}

all_corvals <- eigencorvalue(pcaobj=p, metavars=colnames(metadata), corFUN = "pearson")

# get the top 20 tech covariates that are most correlated with the PCs, in this case we only have 13, so we use all of them
# Extract the correlation values
correlation_matrix <- all_corvals$correlation_values

# Calculate the absolute values of correlation coefficients
abs_correlation_values <- abs(correlation_matrix)

# Flatten the matrix into a vector (since we have multiple x and y variables)
cor_vector <- as.vector(abs_correlation_values)

# Find the indices of the top 13 values
top_indices <- order(cor_vector, decreasing = TRUE)

# Get the corresponding variable names
x_names <- all_corvals$x_names
y_names <- all_corvals$y_names

# Extract the top 13 technical covariates and their names
top_covariates <- data.frame(
  x_var = rep(x_names, each = length(y_names)),
  y_var = rep(y_names, times = length(x_names)),
  correlation = cor_vector,
  stringsAsFactors = FALSE
)

# Sort the top 13 covariates by correlation value
top_covariates <- top_covariates[top_indices, ]

# Print or return the top 13 technical covariates
  top_13 <- unique(top_covariates$y_var)
top_13 
###############################################################################################################################################
#visualize how each tech cov correlates with the PCs
eigencor_plot_p13 <- eigencorplot(p, metavars = top_13, 
                                                 cexLabY = 0.5, rotLabY = 0.8, corFUN = "pearson", 
                                                 main = "Correlation of PCs with technical covariates in CMC and significancies", 
                                                 titleX = "PCs", titleY = "technical covariates", 
                                                 corMultipleTestCorrection = "hochberg", 
                                                 signifSymbols = c('***', '**', '*', ''), signifCutpoints = c(0, 0.001, 0.01, 0.05, 1), scale=FALSE)
#which batch info should we use?
## use PCA to check the the variance of the gene expression that is contributed by the batch effect
pca_result <- prcomp(t(zscore_data), scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
# Assuming you have a vector 'batch_info' containing batch information for each sample
pca_df$Batch <- MSSM_meta_combined[[batch_list[5]]]
print(batch_list[5])
# Plot PCA and color dots by sample batch
plot_pca <- ggplot(pca_df, aes(PC1, PC2, color = Batch)) +
  geom_point() +
  labs(title = "PCA Plot with Batch Effect") +
  theme_minimal()

#no obvious clustering due to batch effect, so we can use the library batch as the batch info
zscore_data_removedBatchEff <- removeBatchEffect(
  x = zscore_data, 
  batch = as.vector(MSSM_meta_combined_noNA_rmVar[,"Library_Batch"]), 
  covariates = metadata_nobatch, 
  design = model.matrix(~ MSSM_meta_combined_noNA_rmVar$ageDeath + 
                               MSSM_meta_combined_noNA_rmVar$primaryDiagnosis + 
                               MSSM_meta_combined_noNA_rmVar$reportedGender + 
                               MSSM_meta_combined_noNA_rmVar$race, data=MSSM_meta_combined_noNA_rmVar)
) #other conditions to include in design to preserve the biological variance??

# #######################################################################################################################################
# ### Draw heatmap to show how different tech covariates contribute to the variance of the gene expression
# rownames(MSSM_meta) <- MSSM_meta$individualID
# p_removedBatchEff_cov <- PCAtools::pca(genes_only_new, metadata = MSSM_meta)

# tech_covar_list <- colnames(MSSM_meta)[complete.cases(t(MSSM_meta))][2:31]
# eigencor_plot_removedBatchEff_cov <- eigencorplot(p_removedBatchEff_cov, metavars = tech_covar_list, cexLabY = 0.5, rotLabY = 0.8, corFUN = "pearson", main = "Correlation of PCs with technical covariates in CMC and significancies", titleX = "PCs", titleY = "technical covariates", corMultipleTestCorrection = "hochberg", signifSymbols = c('***', '**', '*', ''), signifCutpoints = c(0, 0.001, 0.01, 0.05, 1), scale=FALSE)
# ################################################################################################################################################################################

# genes_only_corrected <- removeBatchEffect(x=genes_only_new, batch=as.vector(MSSM_meta[,"RIN"])) 
# genes_only_corrected$gene_symbol <- rowname

# Create a new data frame with the row names as a column
zscore_data_removedBatchEff <- data.frame(gene_symbol = rownames(zscore_data_removedBatchEff), zscore_data_removedBatchEff)
#save(zscore_data_removedBatchEff,new_marker_list, file='/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/MSSM_zscore_data_removedBatchEff_and_makers_list.RData')
estimations_corrected = mgpEstimate(
  exprData = zscore_data_removedBatchEff,
  genes = new_marker_list,
  geneColName = 'gene_symbol',
  outlierSampleRemove = F, # should outlier samples removed. This is done using boxplot stats
  geneTransform = NULL, # this is the default option for geneTransform
  groups = NULL, # if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority = F)

#saveRDS(estimations_corrected, '/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/MSSM_estimations_corrected.RData')
#estimations <- readRDS('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/MSSM_estimations_corrected.RData')
#Coerce estimations list into data frame 
estimations_corrected_scaled <- estimations_corrected$estimates %>% as.data.frame() %>% scale() %>% as.data.frame() %>% tibble::rownames_to_column(var = "standardID")
#Merge cell type proportions with sample metadata
estimations_corrected_metadata = right_join(combined_metadata, estimations_corrected_scaled )

#Remove '+' from ageDeath for modelling
estimations_corrected_metadata$ageDeath = as.numeric(gsub("[+]", "", estimations_corrected_metadata$ageDeath)) #308 54

mgp_corrected_estimations = rbind(mgp_estimations, estimations_corrected_metadata) 
#cleaning the data so that the sample names and cell types are consistent in both datasets with sames order
##samples
###map individualID to specimenID
common_samples <- intersect(estimations_corrected_metadata$specimenID, MSSM_snCTP$subject_id) #90 common samples
rownames(estimations_corrected_metadata) <- estimations_corrected_metadata$specimenID

estimations_corrected_metadata_common <- estimations_corrected_metadata[estimations_corrected_metadata$specimenID %in% common_samples,] #90 54
MSSM_snCTP_common_1 <- MSSM_snCTP[MSSM_snCTP$subject_id %in% common_samples,]
rownames(estimations_corrected_metadata) <- estimations_corrected_metadata$specimenID
rownames(MSSM_snCTP_common_1) <- MSSM_snCTP_common_1$subject_id
estimations_corrected_metadata_common <- estimations_corrected_metadata_common[,36:54]
MSSM_snCTP_common_1 <- MSSM_snCTP_common_1[,3:21]
##cell types 
MSSM_est_corrected_common <- estimations_corrected_metadata_common
new_order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,15,17,18,19)
MSSM_est_corrected_common <- MSSM_est_corrected_common[,new_order]
colnames(MSSM_est_corrected_common) <- colnames(MSSM_snCTP_common_1)
# Dan's estimation: see above scatter plot
# predicted_df <- data.frame(Sample = rownames(MSSM_est_common), MSSM_est_common)
# actual_df <- data.frame(Sample = rownames(MSSM_snCTP_common), MSSM_snCTP_common)
# estimations_scaled  <- read.csv("/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/psychencode_mgp_estimations.csv")
################################################################################################################################################################################

## Compare if batch correction improves the prediction accuracy compared with orignal MGP using barplot
### helper function
cor_func <- function(estimations_df, actual_prop) {
  predicted_df <- data.frame(Sample = rownames(estimations_df), estimations_df)
  actual_df <- data.frame(Sample = rownames(actual_prop), actual_prop)

  # Merge the predicted and actual data frames
  merged_df <- merge(predicted_df, actual_df, by = "Sample")

  accuracy_list <- c()
  for (cell_type in colnames(estimations_df)) {
    # Calculate the prediction accuracy (e.g., pearson correlation coefficient)
    x <- paste(cell_type, ".x", sep = "")
    y <- paste(cell_type, ".y", sep = "")
    accuracy <- cor(merged_df[[x]], merged_df[[y]], method = "pearson", use = "pairwise.complete.obs")
    accuracy_list <- c(accuracy_list, accuracy)
  }
  return(accuracy_list)
}

mgp_corrected_sn_cor <- cor_func(MSSM_est_corrected_common, MSSM_snCTP_common_1)

# common_samples <- intersect(MSSM_est$specimenID, MSSM_snCTP$subject_id) #92 common samples -> 79 after remove low cell counts samples ie: <500
# MSSM_est_common <- MSSM_est[MSSM_est$specimenID %in% common_samples,] #92 55
# MSSM_snCTP_common <- MSSM_snCTP[MSSM_snCTP$subject_id %in% common_samples,] #92 24
# rownames(MSSM_est_common) <- MSSM_est_common$specimenID
# rownames(MSSM_snCTP_common) <- MSSM_snCTP_common$subject_id
# MSSM_est_common <- MSSM_est_common[,37:55]
# MSSM_snCTP_common <- MSSM_snCTP_common[,3:21]
# ##cell types  
# new_order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,15,17,18,19)
# MSSM_est_common <- MSSM_est_common[,new_order]
# colnames(MSSM_est_common) <- colnames(MSSM_snCTP_common)

mgp_sn_cor <- cor_func(MSSM_est_common, MSSM_snCTP_common)

#subset for PVALB, SST, VIP
accuracy_list <- c(mgp_corrected_sn_cor, mgp_sn_cor)

#####
categories <- c("mgp_corrected", "mgp_original")
names <- colnames(MSSM_snCTP_common)

category_vector <- rep(categories, each = 19)
name_vector <- rep(names, 2)
df <- data.frame(Name = name_vector, Category = category_vector, Accuracy = accuracy_list)

color <- brewer.pal(length(categories), "Set3")
# Use ggplot2 to create the barplot.
p <- ggplot(df, aes(x = Name, y = Accuracy, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = color <- brewer.pal(5, "Set1")) +
  coord_flip() +
  labs(title = "Pearson Correlaton Coefficient for Cell Type Proportion Prediction by MGP w/withou Batch Correction",
       x = "Cell Type and parameter settings",
       y = "Pearson Correlation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

print(p)

# #save into png file
# png('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/barplots_for_three_diff_tech_cov_reg.png')
# p
# dev.off()










# ########################################################################################
# # We don't have detailed tech covariates for the snCTP data, so we can only do PCA with the current ones.
# ####RNA-seq alignment and Quality control PCA-based outlier identification
# ### PCA for RAW count matrix without normalization and covariate correction
# #make sure data is matrix, not a data frame
# count.matix <- as.matrix() # raw count
# rosmap_meta_DLPFC_cleaned <- rosmap_meta_DLPFC[!is.na(rosmap_meta_DLPFC$RIN), ] #remove the samples with NA RIN
# count.matix_cleaned <- count.matix[, !(colnames(count.matix) %in% rosmap_meta_DLPFC[is.na(rosmap_meta_DLPFC$RIN), ]$synapseID)]
# dds <- DESeqDataSetFromMatrix(countData=count.matix, colData=rosmap_meta_DLPFC, design=~1) #no design?
# #dds <- DESeqDataSetFromMatrix(countData=count.matix_cleaned, colData=rosmap_meta_DLPFC_cleaned, design=~RIN+sequencingBatch+libraryPrep+readLength)
# # Normalize and transform the data in the `DESeqDataSet` object
# # using the `vst()` function from the `DESeq2` R package
# rld <- vst(dds, blind=FALSE)
# #get matrix of transformed data value
# mat <- assay(rld)

# #PCA using gglot
# tech_covar_list <- colnames(rosmap_meta_DLPFC_cleaned)[complete.cases(t(rosmap_meta_DLPFC_cleaned))]
# pcaData <- plotPCA(rld, intgroup=tech_covar_list, returnData=TRUE)
# # intgroup(): This parameter specifies the variables to use for plotting. The PCA will be based on these variables, and the plot will show the data points (samples) in the space defined by these variables. The variables should be categorical variables and tech covariates, and the plot will show the data points colored by these variables.

# ####################################### calculate HERE!!!!!!!
# #compute the variance for PCs
# #prcomp(rld)
# #######################################



# percentVar <- round(100*attr(pcaData, "percentVar"))
# g <- ggplot(pcaData, aes(PC1,PC2))
# # plot the PCA
# g <- g + geom_point(size=3)
# # colored by sequecing batch
# g1 <- g + geom_point(size=3, aes(col=RnaSeqMetrics__CODING_BASES, fill=RnaSeqMetrics__CODING_BASES))
# # plot(g)
# # plot(g1)

# ### Draw a line where PC1 is more than 4 sd from the mean on the plot
# # Calculate mean and standard deviation for PC1
# pc1_mean <- mean(pcaData$PC1)
# pc1_sd <- sd(pcaData$PC1)
# # Define the threshold (e.g., 4 standard deviations from the mean)
# threshold <- pc1_mean - 4 * pc1_sd
# # Create the ggplot object
# g <- ggplot(pcaData, aes(PC1, PC2))
# # Plot the PCA
# g <- g + geom_point(size = 3)
# # Add a horizontal line for the threshold
# g <- g + geom_vline(xintercept = threshold, color = "red", linetype = "dashed")
# # Display the plot
# # print(g)
# # No outlier detected

# ######################################################################################################################\
# ### Covariate correction
# # Remove genes with no variation
# dds_hasvar <- dds[rowVars(assay(dds)) > 0.1,]  #60603 -> 41749 genes

# # Perform TMM normalization and log2 transformation
# dds_hasvar <- DESeq(dds_hasvar)
# # rld <- rlogTransformation(dds_hasvar): log transformation takes too long for more than 50 samples
# rld <- vst(dds_hasvar, blind=FALSE)
# #saveRDS(rld, '/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_rld.RData')
# #sbatch -J ROSMAP_rld --mem=4G -N 2 -t 0-2:0 -o /external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_rld.slurm.log --wrap="Rscript ls/xzhou/GWAS_tut/ROSMAP/PCA_tech_cov.r"

# # Center genes
# centered_data <- scale(assay(rld), center = TRUE, scale = FALSE)

# # Z-score transform RNA-seq counts per sample
# zscore_data <- t(scale(t(centered_data)))
# saveRDS(zscore_data, '/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_zscore_data.RData')
# zscore_data <- readRDS('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_zscore_data.RData')

# # PCA again to find tech cov
# res.pca <- prcomp(zscore_data)

# #visualize tge scree plot
# fviz_eig(res.pca)





# zscore_data_removedBatchEff <- removeBatchEffect(x=zscore_data_removedBatchEff, batch=as.vector(metadata[,"sequencingBatch"]), batch2=as.vector(metadata[,"libraryPrep"]))

#sbatch -J sep_bench --mem=4G -N 2 -t 0-2:0 -o /external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/sep_bench.slurm.log --wrap="Rscript /nethome/kcni/xzhou/cell_deconv/benchmark_separated_cohorts.r"