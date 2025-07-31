#import packages
library(GEOquery)
library(limma)
library(ggplot2)
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
library(DESeq2)
library(factoextra)
library(PCAtools)
library(markerGeneProfile)
library(edgeR)
library(readr)
library(tidyverse)


# read sibille lab cell densities
sibille_lab_cell_fractions = read_csv('/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/sibille_lab_cell_fractions.csv')

# read ACC sample info - DAN you need to download these from synapse and put them in an appropriate place on the SCC. i've just put them in my local folder for now
cmc_clinical_meta = read_csv('/external/rprshnas01/external_data/CommonMind/ControlledAccess/Data/Clinical/CMC_Human_clinical_metadata.csv')
cmc_acc_rnaseq_meta = read_csv('/external/rprshnas01/external_data/psychencode/PsychENCODE/CMC_ACC/CMC_Human_rnaSeq_metadata_release6.csv')

# figure out overlapping subject IDs between sibille lab densities and cmc ACC
cmc_acc_overlapping = cmc_acc_rnaseq_meta %>% filter(Individual_ID %in% sibille_lab_cell_fractions$Individual.ID)

# merge data from overlapping subjects with clinical metadata
cmc_acc_overlapping_w_meta = left_join(cmc_acc_overlapping, cmc_clinical_meta %>% mutate(Individual_ID = `Individual ID`))

# count up subjects with cell densities from etienne and rnaseq data from acc from commonmind
cmc_acc_overlapping_w_meta %>% group_by(Dx) %>% tally()  

#remove duplicated samples
sibille_lab_cell_fractions <- sibille_lab_cell_fractions %>% distinct(Individual.ID, .keep_all = TRUE)
# merge the cell fractions with the cmc_acc_overlapping_w_meta
CMC_ACC_cell_fractions_combined <- left_join(cmc_acc_overlapping_w_meta, sibille_lab_cell_fractions, by = c("Individual_ID" = "Individual.ID")) #[1] 51 samples 58

#########################################################
# read in the MGPs from CMC_ACC from dan's analysis
CMC_ACC_MGPs <- read_csv('/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/cmc_acc_mgp_estimations.csv')
CMC_ACC_MGPs_common <- CMC_ACC_MGPs[,colnames(CMC_ACC_MGPs) %in% c("SampleID","PVALB","VIP","SST")]
##########################################################
CMC_ACC_cell_fractions_common <- CMC_ACC_cell_fractions_combined[,colnames(CMC_ACC_cell_fractions_combined) %in% c("SampleID","PVALB","VIP","SST")]
# make sure these two dataframes have common samples and the same sample order
common_samples <- intersect(CMC_ACC_MGPs_common$SampleID, CMC_ACC_cell_fractions_common$SampleID)
CMC_ACC_MGPs_common <- CMC_ACC_MGPs_common[CMC_ACC_MGPs_common$SampleID %in% common_samples,]
CMC_ACC_cell_fractions_common <- CMC_ACC_cell_fractions_common[CMC_ACC_cell_fractions_common$SampleID %in% common_samples,]
# remove NA in cell density??


##########################################################################################
#calculate the correlation between the estimated cell proportions and the snCTP cell proportions
predicted_df <- CMC_ACC_MGPs_common
actual_df <- CMC_ACC_cell_fractions_common 

merged_df <- merge(predicted_df, actual_df, by = "SampleID")

plots_list <- list()

for (cell_type in colnames(CMC_ACC_cell_fractions_common[,2:4])) {
  # Calculate the prediction accuracy (e.g., correlation coefficient)
  x <- paste(cell_type, ".x", sep = "")
  y <- paste(cell_type, ".y", sep = "")
  accuracy <- cor(merged_df[[x]], merged_df[[y]], use = "pairwise.complete.obs")

  # Create the scatter plot with a trend line
  plot <- ggplot(merged_df, aes(x = .data[[y]], y = .data[[x]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
    labs(x = "Actual Proportions in Cell Density", y = "Predicted Proportions by MGP") +
    ggtitle(paste("Prediction Accuracy for", cell_type, ":", round(accuracy, 2)))

  # Store the scatter plot in the list
  plots_list[[cell_type]] <- plot
}

# Arrange and display the scatter plots in a grid
png("/nethome/kcni/xzhou/label_transfer/benchmark_CMC_ACC_MGPs.png", width = 800, height = 600)
grid.arrange(grobs = plots_list, ncol = 3)
dev.off()
#############

#remove NA only for SST cells
CMC_ACC_MGPs_common_noNA_SST <- CMC_ACC_MGPs_common[!is.na(CMC_ACC_MGPs_common$SST), colnames(CMC_ACC_MGPs_common) %in% c("SampleID","SST")]
CMC_ACC_cell_fractions_common_noNA_SST <- CMC_ACC_cell_fractions_common[!is.na(CMC_ACC_cell_fractions_common$SST), colnames(CMC_ACC_cell_fractions_common) %in% c("SampleID","SST")] #51 -> 46samples

###############################################
# #manually remove the outlier
# CMC_ACC_cell_fractions_common_noNA_SST <- CMC_ACC_cell_fractions_common_noNA_SST[CMC_ACC_cell_fractions_common_noNA_SST$SampleID != "PITT_RNA_ACC_1367",] #45

# #match the sample between the two dataframes
# CMC_ACC_MGPs_common_noNA_SST <- CMC_ACC_MGPs_common_noNA_SST[CMC_ACC_MGPs_common_noNA_SST$SampleID %in% CMC_ACC_cell_fractions_common_noNA_SST$SampleID,]


# predicted_df <- CMC_ACC_MGPs_common_noNA_SST
# actual_df <- CMC_ACC_cell_fractions_common_noNA_SST 

# merged_df <- merge(predicted_df, actual_df, by = "SampleID")

# # Calculate the prediction accuracy (e.g., correlation coefficient)
# x <- paste("SST.x", sep = "")
# y <- paste("SST.y", sep = "")
# accuracy <- cor(merged_df[[x]], merged_df[[y]], use = "pairwise.complete.obs")

# # Create the scatter plot with a trend line
# plot <- ggplot(merged_df, aes(x = .data[[y]], y = .data[[x]])) +
# geom_point() +
# geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
# labs(x = "Actual Proportions in Cell Density", y = "Predicted Proportions by MGP") +
# ggtitle(paste("Prediction Accuracy for SST", ":", round(accuracy, 2)))
# # Arrange and display the scatter plots in a grid
# png("/nethome/kcni/xzhou/label_transfer/benchmark_CMC_ACC_MGPs_SST.png", width = 800, height = 600)
# plot
# dev.off()
############################################################


# Function to remove outliers beyond 3 standard deviations
remove_outliers <- function(x) {
  for (cell_type in colnames(x)[2:length(colnames(x))]) {
    mean_cell_type <- mean(x[[cell_type]], na.rm = TRUE)
    sd_cell_type <- sd(x[[cell_type]], na.rm = TRUE)
    lower_bound <- mean_cell_type - 3 * sd_cell_type
    upper_bound <- mean_cell_type + 3 * sd_cell_type
    x_filtered <- x[x[[cell_type]] >= lower_bound & x[[cell_type]] <= upper_bound, ]
  }
  return(x_filtered)
}

# Remove outliers from CMC_ACC_MGPs_common for SST
#CMC_ACC_MGPs_common_noNA_SST_filtered <- remove_outliers(CMC_ACC_MGPs_common_noNA_SST)

# Remove outliers from CMC_ACC_cell_fractions_common for SST
CMC_ACC_cell_fractions_common_noNA_SST_filtered <- remove_outliers(CMC_ACC_cell_fractions_common_noNA_SST)

# Ensure both dataset have identical samples
order_common_samples <- function(df1, df2) {
  common_samples <- intersect(df1[[colnames(df1)[1]]], df2[[colnames(df2)[1]]])
  print(df1[[colnames(df1)[1]]])
  df1_common <- df1[df1[[colnames(df1)[1]]] %in% common_samples, ]
  df2_common <- df2[df2[[colnames(df2)[1]]] %in% common_samples, ]
  return(list(df1_common, df2_common))  # Return filtered data frames
}

# Call the function and assign filtered common data frames
common_data_frames <- order_common_samples(CMC_ACC_MGPs_common_noNA_SST, CMC_ACC_cell_fractions_common_noNA_SST_filtered)
CMC_ACC_MGPs_noNA_SST_filtered_common <- common_data_frames[[1]]
CMC_ACC_cell_fractions_noNA_SST_filtered_common <- common_data_frames[[2]]

#multiply cell density by 10 to achieve cells/mm^2
CMC_ACC_cell_fractions_noNA_SST_filtered_common$SST <- CMC_ACC_cell_fractions_noNA_SST_filtered_common$SST*10
#plot
predicted_df <- CMC_ACC_MGPs_noNA_SST_filtered_common 
actual_df <- CMC_ACC_cell_fractions_noNA_SST_filtered_common

merged_df <- merge(predicted_df, actual_df, by = "SampleID")

# Calculate the prediction accuracy (e.g., correlation coefficient)
x <- paste("SST.x", sep = "")
y <- paste("SST.y", sep = "")
accuracy <- cor(merged_df[[x]], merged_df[[y]], use = "pairwise.complete.obs")

# Create the scatter plot with a trend line
plot_SST <- ggplot(merged_df, aes(x = .data[[y]], y = .data[[x]])) +
geom_point() +
geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
labs(x = "Actual Proportions in Cell Density", y = "Predicted Proportions by MGP") +
ggtitle(paste("Prediction Accuracy for SST", ":", round(accuracy, 2)))
# Arrange and display the scatter plots in a grid
png("/nethome/kcni/xzhou/label_transfer/benchmark_CMC_ACC_MGPs_SST_3SD.png", width = 800, height = 600)
plot
dev.off()
###################################

#remove ourlier for both VIP and PVALB
# Call the function and assign filtered common data frames
CMC_ACC_MGPs_common_noNA_VIP <- CMC_ACC_MGPs_common[!is.na(CMC_ACC_MGPs_common$VIP), colnames(CMC_ACC_MGPs_common) %in% c("SampleID","VIP")]
CMC_ACC_cell_fractions_common_noNA_VIP <- CMC_ACC_cell_fractions_common[!is.na(CMC_ACC_cell_fractions_common$VIP), colnames(CMC_ACC_cell_fractions_common) %in% c("SampleID","VIP")]
# Remove outliers from CMC_ACC_MGPs_common for VIP
#CMC_ACC_MGPs_common_noNA_VIP_filtered <- remove_outliers(CMC_ACC_MGPs_common_noNA_VIP)
CMC_ACC_MGPs_common_noNA_VIP_filtered <- CMC_ACC_MGPs_common_noNA_VIP
# Remove outliers from CMC_ACC_cell_fractions_common for VIP
CMC_ACC_cell_fractions_common_noNA_VIP_filtered <- remove_outliers(CMC_ACC_cell_fractions_common_noNA_VIP)
common_data_frames <- order_common_samples(CMC_ACC_MGPs_common_noNA_VIP_filtered, CMC_ACC_cell_fractions_common_noNA_VIP_filtered)
CMC_ACC_MGPs_noNA_VIP_filtered_common <- common_data_frames[[1]]
CMC_ACC_cell_fractions_VIP_filtered_common <- common_data_frames[[2]]

#remove ourlier for both VIP and PVALB
# Call the function and assign filtered common data frames
CMC_ACC_MGPs_common_noNA_PVALB <- CMC_ACC_MGPs_common[!is.na(CMC_ACC_MGPs_common$PVALB), colnames(CMC_ACC_MGPs_common) %in% c("SampleID","PVALB")]
CMC_ACC_cell_fractions_common_noNA_PVALB <- CMC_ACC_cell_fractions_common[!is.na(CMC_ACC_cell_fractions_common$PVALB), colnames(CMC_ACC_cell_fractions_common) %in% c("SampleID","PVALB")]
# Remove outliers from CMC_ACC_MGPs_common for PVALB
#CMC_ACC_MGPs_common_noNA_PVALB_filtered <- remove_outliers(CMC_ACC_MGPs_common_noNA_PVALB)
CMC_ACC_MGPs_common_noNA_PVALB_filtered  <- CMC_ACC_MGPs_common_noNA_PVALB
# Remove outliers from CMC_ACC_cell_fractions_common for PVALB
CMC_ACC_cell_fractions_common_noNA_PVALB_filtered <- remove_outliers(CMC_ACC_cell_fractions_common_noNA_PVALB)
common_data_frames <- order_common_samples(CMC_ACC_MGPs_common_noNA_PVALB_filtered, CMC_ACC_cell_fractions_common_noNA_PVALB_filtered)
CMC_ACC_MGPs_noNA_PVALB_filtered_common <- common_data_frames[[1]]
CMC_ACC_cell_fractions_PVALB_filtered_common <- common_data_frames[[2]]

# ensure three df have identical samples:
common_all <- intersect(CMC_ACC_MGPs_noNA_PVALB_filtered_common$SampleID, CMC_ACC_MGPs_noNA_SST_filtered_common$SampleID)
common_all <- intersect(common_all, CMC_ACC_MGPs_noNA_VIP_filtered_common$SampleID)
#common_test <- readRDS(file="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/common_samples.RData")
## filter for SST
CMC_ACC_MGPs_noNA_SST_filtered_common <- CMC_ACC_MGPs_noNA_SST_filtered_common[CMC_ACC_MGPs_noNA_SST_filtered_common$SampleID %in% common_all,]
CMC_ACC_cell_fractions_noNA_SST_filtered_common <- CMC_ACC_cell_fractions_noNA_SST_filtered_common[CMC_ACC_cell_fractions_noNA_SST_filtered_common$SampleID %in% common_all,]
## filter for VIP
CMC_ACC_MGPs_noNA_VIP_filtered_common <- CMC_ACC_MGPs_noNA_VIP_filtered_common[CMC_ACC_MGPs_noNA_VIP_filtered_common$SampleID %in% common_all,]
CMC_ACC_cell_fractions_VIP_filtered_common <- CMC_ACC_cell_fractions_VIP_filtered_common[CMC_ACC_cell_fractions_VIP_filtered_common$SampleID %in% common_all,]
##filter for PVALB
CMC_ACC_MGPs_noNA_PVALB_filtered_common <- CMC_ACC_MGPs_noNA_PVALB_filtered_common[CMC_ACC_MGPs_noNA_PVALB_filtered_common$SampleID %in% common_all,]
CMC_ACC_cell_fractions_PVALB_filtered_common <- CMC_ACC_cell_fractions_PVALB_filtered_common[CMC_ACC_cell_fractions_PVALB_filtered_common$SampleID %in% common_all,]


#multiply cell density by 10 to achieve cells/mm^2
CMC_ACC_cell_fractions_common_noNA_VIP_filtered$VIP <- CMC_ACC_cell_fractions_common_noNA_VIP_filtered$VIP*10
CMC_ACC_cell_fractions_common_noNA_PVALB_filtered$PVALB <- CMC_ACC_cell_fractions_common_noNA_PVALB_filtered$PVALB*10
CMC_ACC_cell_fractions_common_noNA_SST_filtered$SST <- CMC_ACC_cell_fractions_common_noNA_SST_filtered$SST*10
######
#merge CMC_ACC_MGPs_noNA_PVALB_filtered_common, CMC_ACC_MGPs_noNA_SST_filtered_common, and CMC_ACC_MGPs_noNA_VIP_filtered_common into predict_df
predicted_df <- merge(CMC_ACC_MGPs_noNA_PVALB_filtered_common, CMC_ACC_MGPs_noNA_SST_filtered_common, by = "SampleID")
predicted_df <- merge(predicted_df, as.data.frame(CMC_ACC_MGPs_noNA_VIP_filtered_common), by = "SampleID")
predicted_df <- data.frame(Sample = predicted_df$SampleID, predicted_df[,2:4])
predicted_df_long <- pivot_longer(predicted_df, cols = c(2:4), names_to = "cell_type", values_to = "rCTP")

actual_df <- merge(CMC_ACC_cell_fractions_common_noNA_PVALB_filtered, CMC_ACC_cell_fractions_common_noNA_SST_filtered, by = "SampleID")
actual_df <- merge(actual_df, as.data.frame(CMC_ACC_cell_fractions_common_noNA_VIP_filtered), by = "SampleID")
actual_df <- data.frame(Sample = actual_df$SampleID, actual_df[,2:4])
actual_df_long <- pivot_longer(actual_df, cols = c(2:4), names_to = "cell_type", values_to = "snCTP")

#merge the two dataframes into a single one using sampls and cell_type as keys
merged_df <- merge(predicted_df_long, actual_df_long, by = c("Sample", "cell_type")) # for 3 cell types, so 43 samples for each cell type

# Step 1: Calculate pearson r for each cell type
correlation_df <- merged_df %>%
  group_by(cell_type) %>%
  summarize(cor_value = cor(snCTP, rCTP, method = "pearson")) %>%
  ungroup()

# Example of what correlation_df might look like
correlation_df$x_pos <- c(0, 6, 2)
correlation_df <- as.data.frame(correlation_df)

# # Your plotting code, now including geom_text for the correlation coefficient
# figure_S2 <- merged_df %>%
#   ggplot(aes(x = snCTP, y = rCTP)) +
#   geom_smooth(se = FALSE, method = 'lm', fullrange = TRUE) +
#   geom_point(alpha = 1, size = 1.5, color="#56B4E9") +
#   stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~")), color = "black", geom = "label")+
#   # geom_text(
#   #   data = correlation_df, 
#   #   aes(label = sprintf("Pearson's R = %.2f", cor_value), x = x_pos, y = 3), # You'll need to adjust these values based on your actual plot range
#   #   hjust = 0, vjust = 1, 
#   #   size = 5, color = "blue",
#   #   inherit.aes = FALSE # This is to avoid inheriting the other aesthetics like fill or group
#   # ) +
#   ylab('rCTP(AU)') + 
#   xlab('Cell Density(cells/mm^2)') +
#   theme(
#     panel.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     strip.background = element_rect(fill = "white", colour = "black"),
#     axis.line = element_line(color = "black"),
#     axis.text.x = element_text(size = 15),
#     axis.text.y = element_text(size = 15),
#     axis.title.x = element_text(vjust = -3,size = 15),
#     plot.margin = margin(1, 1, 1, 1, "cm"),
#     axis.title.y = element_text(size = 15),
#     legend.position = c(0.9, 0.9),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 12),
#     strip.text.x = element_text(size = 14)
#   ) +
#   facet_grid(~cell_type, scale="free_x") #+
#   #scale_y_continuous(limits = c(-3, 3))

# print(figure_S2)



figure_S2 <- merged_df %>%
  ggplot(aes(x = snCTP, y = rCTP)) +
  geom_smooth(se = FALSE, method = 'lm', fullrange = TRUE) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~")), color = "black", geom = "label", size=10, label.y=4)+
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
  xlab('Cell Density(cells/mm^2)') +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 28),
    axis.text.y = element_text(size = 28),
    axis.title.x = element_text(vjust = -3,size = 28),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title.y = element_text(size = 28),
    legend.position = c(0.9, 0.9),
    legend.title = element_blank(),
    legend.text = element_text(size = 28),
    strip.text.x = element_text(size = 28)
  ) +
  facet_grid(~cell_type, scale="free_x") +
  theme(panel.spacing = unit(1, "cm", data = NULL)) +
  scale_y_continuous(limits = c(-4, 4))

#print(figure_S2)

#Save to plots folder
output_filename <- "/nethome/kcni/xzhou/dan_paper/figure_S2.png"

# Use ggsave to save the plot as a png image
ggsave(output_filename, figure_S2, width = 15, height = 8, dpi = 200)






#########################################
#plot
predicted_df <- CMC_ACC_MGPs_noNA_VIP_filtered_common 
actual_df <- CMC_ACC_cell_fractions_common_noNA_VIP_filtered

merged_df <- merge(predicted_df, actual_df, by = "SampleID")

# Calculate the prediction accuracy (e.g., correlation coefficient)
x <- paste("VIP.x", sep = "")
y <- paste("VIP.y", sep = "")
accuracy <- cor(merged_df[[x]], merged_df[[y]], use = "pairwise.complete.obs")

# Create the scatter plot with a trend line
plot <- ggplot(merged_df, aes(x = .data[[y]], y = .data[[x]])) +
geom_point() +
geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
labs(x = "Actual Proportions in Cell Density", y = "Predicted Proportions by MGP") +
ggtitle(paste("Prediction Accuracy for VIP", ":", round(accuracy, 2)))
# Arrange and display the scatter plots in a grid
png("/nethome/kcni/xzhou/label_transfer/benchmark_CMC_ACC_MGPs_VIP.png", width = 800, height = 600)
plot
dev.off()
#######################

CMC_ACC_cell_fractions_common_noNA_PVALB_filtered$PVALB <- CMC_ACC_cell_fractions_common_noNA_PVALB_filtered$PVALB*10
#plot
predicted_df <- CMC_ACC_MGPs_noNA_PVALB_filtered_common 
actual_df <- CMC_ACC_cell_fractions_common_noNA_PVALB_filtered

merged_df <- merge(predicted_df, actual_df, by = "SampleID")

# Calculate the prediction accuracy (e.g., correlation coefficient)
x <- paste("PVALB.x", sep = "")
y <- paste("PVALB.y", sep = "")
accuracy <- cor(merged_df[[x]], merged_df[[y]], use = "pairwise.complete.obs")

# Create the scatter plot with a trend line
plot <- ggplot(merged_df, aes(x = .data[[y]], y = .data[[x]])) +
geom_point() +
geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
labs(x = "Actual Proportions in Cell Density", y = "Predicted Proportions by MGP") +
ggtitle(paste("Prediction Accuracy for PVALB", ":", round(accuracy, 2)))
# Arrange and display the scatter plots in a grid
png("/nethome/kcni/xzhou/label_transfer/benchmark_CMC_ACC_MGPs_PVALB.png", width = 800, height = 600)
plot
dev.off()