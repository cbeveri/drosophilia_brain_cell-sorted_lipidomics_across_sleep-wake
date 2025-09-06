# Modified R Script for FA, DAG, TAG Analysis


library(plyr)
library(ggridges)
library(tidyverse)
library(readxl)
library(edgeR)
library(tidyverse)
library(factoextra)
library(grid)
library(ggsignif)
library(stringr)
library(reshape2)
library(gtools)
library(tidyr)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(ggplot2)
library(cowplot)
library(ggridges)
library(ggforce)

# Function to plot stacked bars with blank normalization for FA, DAG, TAG only
plot_combined_values_Stacked_with_blank_normalized_filtered <- function(df, Title1, Title2) {
  # Filter lipid classes to only FA, DAG, TAG
  lipid_classes <- c("FA", "DAG", "TAG")
  lipid_colors <- c("#33a02c", "#8dd3c7", "#6a3d9a")
  
  lipid_colors_alpha <- alpha(lipid_colors, 0.5)
  lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
  
  # Filter the dataframe to only include FA, DAG, TAG
  df <- df %>% filter(type %in% lipid_classes)
  
  # Extract column lengths
  len1 <- as.numeric(df$Length1[1])
  len2 <- as.numeric(df$Length2[1])
  
  # Columns after 'type' are the ones of interest
  all_cols <- colnames(df)
  start_idx <- which(all_cols == "type") + 1
  cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
  cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
  
  # Get the blank column (which is the penultimate column)
  blank_col <- all_cols[length(all_cols)]
  
  # Subtract the blank column from the columns of interest and set negative values to 0
  df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
  df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
  
  df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
  df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
  
  # Compute the means for the groups of columns
  df <- df %>%
    mutate(mean1 = rowMeans(select(., one_of(cols1))),
           mean2 = rowMeans(select(., one_of(cols2))))
  
  # Divide the means by the sum of their respective columns
  total_mean1 <- sum(df$mean1)
  total_mean2 <- sum(df$mean2)
  
  df$mean1 <- df$mean1 / total_mean1
  df$mean2 <- df$mean2 / total_mean2
  
  # Sum by type for all values
  df_all <- df %>%
    group_by(type) %>%
    summarise(sum_mean1 = sum(mean1), 
              sum_mean2 = sum(mean2), .groups = "keep")
  
  # Combine the dataframes for plotting
  df_all_long <- df_all %>% 
    select(type, sum_mean1, sum_mean2) %>% 
    gather(key = "group", value = "value", -type)
  
  # Ensure the order is FA, DAG, TAG from bottom to top
  df_all_long$type <- factor(df_all_long$type, levels = c("FA", "DAG", "TAG"))
  
  # Plotting all values
  plot_all <- df_all_long %>%
    ggplot(aes(x=group, y = value, fill = type)) +
    geom_col() +
    theme_classic() +
    ylab("Total lipid content") +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
    scale_x_discrete(labels = c(Title1, Title2)) +
    labs(title = paste(Title1, "vs", Title2, " (FA, DAG, TAG only)"), y = "Sum of Means")
  
  # Create the new directory for filtered lipids
  dir.create("FA_DG_TG_bars", showWarnings = FALSE)
  
  # Save the plot with modified filename
  ggsave(filename = paste0("FA_DG_TG_bars/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_SUM_all_normalized_DG_FA_TG_only.pdf"), plot = plot_all)
}

# Function to plot pie charts normalized for FA, DAG, TAG only
plot_pie_charts_normalized_filtered <- function(df, Title1, Title2, filename, filename_csv) {
  library(gridExtra)
  
  # Filter lipid classes to only FA, DAG, TAG
  lipid_classes <- c("FA", "DAG", "TAG")
  lipid_colors <- c("#33a02c", "#8dd3c7", "#6a3d9a")
  
  lipid_colors_alpha <- alpha(lipid_colors, 0.5)
  lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
  
  # Filter the dataframe to only include FA, DAG, TAG
  df <- df %>% filter(type %in% lipid_classes)
  
  # Extract column lengths
  len1 <- as.numeric(df$Length1[1])
  len2 <- as.numeric(df$Length2[1])
  
  # Columns after 'type' are the ones of interest
  all_cols <- colnames(df)
  start_idx <- which(all_cols == "type") + 1
  cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
  cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
  
  # Get the blank column (which is the penultimate column)
  blank_col <- all_cols[length(all_cols)]
  
  # Subtract the blank column from the columns of interest and set negative values to 0
  df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
  df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
  
  df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
  df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
  
  # Compute the means for the groups of columns
  df <- df %>%
    mutate(mean1 = rowMeans(select(., one_of(cols1))),
           mean2 = rowMeans(select(., one_of(cols2))))
  
  # Sum by type for significant values
  df_sig <- df %>%
    group_by(type) %>%
    summarise(sum_mean1 = mean(mean1), 
              sum_mean2 = mean(mean2), .groups = "keep")
  write.csv(df_sig, file = filename_csv, row.names = FALSE)
  
  # Ensure the order is FA, DAG, TAG
  df_sig$type <- factor(df_sig$type, levels = c("FA", "DAG", "TAG"))
  
  # Make the pie charts
  pie1 <- df_sig %>%
    ggplot(aes(x = "", y = sum_mean1, fill = type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(axis.text =element_blank(),
          axis.line = element_blank()) +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
    labs(title = Title1)
  
  pie2 <- df_sig %>%
    ggplot(aes(x = "", y = sum_mean2, fill = type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(axis.text =element_blank(),
          axis.line = element_blank()) +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
    labs(title = Title2)
  
  # Save the plots to an SVG file side by side
  svg(filename, width=10, height=5)
  grid.arrange(pie1, pie2, ncol=2)
  dev.off()
  
  return(list(pie1 = pie1, pie2 = pie2))
}

file_list = list.files(path="results", pattern=NULL, all.files=FALSE,
                       full.names=FALSE)
# Main processing loop with filtered functions
for (jj in file_list){
  excel_file <- read_csv(paste0("results/",jj,sep=""))
  
  if (!grepl("full", jj, ignore.case = TRUE)) {
    next
  }
  
  Title1 <- gsub("[:| ]", "_", excel_file$Title_1[1])
  Title2 <- gsub("[:| ]", "_", excel_file$Title_2[1])
  Title1 <- gsub("__", "_", Title1)
  Title2 <- gsub("__", "_", Title2) 
  
  title_for_plot <- paste0(Title1,Title2,sep="_")
  title_for_plot <- gsub("_| |-", "", title_for_plot)
  Title1 <- gsub("_| |-", "", Title1)
  Title2 <- gsub("_| |-", "", Title2)
  
  # Create the new filtered plots
  plot_combined_values_Stacked_with_blank_normalized_filtered(excel_file, Title1, Title2)
  
  # Create the new directory for filtered pie charts
  dir.create("Pie_DG_FA_TG_only", showWarnings = FALSE)
  
  plot_pie_charts_normalized_filtered(excel_file, Title1, Title2, 
                                      paste0("Pie_DG_FA_TG_only/",Title1,Title2,"_normalized_by_number_DG_FA_TG_only.svg",sep=""),
                                      paste0("Pie_DG_FA_TG_only/",Title1,Title2,"_normalized_by_number_DG_FA_TG_only.csv",sep=""))
}

