
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
library(ggbeeswarm) #https://github.com/eclarke/ggbeeswarm
library(ggrepel)
library(scales)

library(ggplot2)
#library(ggExtra) #https://cran.r-project.org/web/packages/ggExtra/vignettes/ggExtra.html

library(cowplot)
library(ggridges)

# library(limma)
# library(writexl)


# install.packages("ggforce")
library(ggforce)



getwd()


# plot_histogram_plots_TGs_DE_lipids <- function(df, Title1, Title2) {
#   # Define lipid classes and colors
#   # lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG", "DAG", "TAG | DAG", "DAG | CE", "TAG | DAG | CE")
#   # lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
#   # lipid_colors_alpha <- scales::alpha(lipid_colors, 0.5)
#   # lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
#   # 
#   # # Extract lengths and column indices
#   # len1 <- as.numeric(df$Length1[1])
#   # len2 <- as.numeric(df$Length2[1])
#   # all_cols <- colnames(df)
#   # start_idx <- which(all_cols == "type") + 1
#   # cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   # cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   # blank_col <- all_cols[length(all_cols)]
#   # 
#   # # Adjust and clean data
#   # df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   # df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   # df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   # df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   # 
#   dir_name <- "up_and_down_regulated_FA_tag_DAG_count"
#   # Further processing
#   df_processed <- df %>%
#     # mutate(mean1 = rowMeans(select(., all_of(cols1))),
#     #        mean2 = rowMeans(select(., all_of(cols2)))) %>%
#     filter(type == "TAG") %>%filter(FDR < 0.1) %>%
#     mutate(lipid_sub = str_extract(lipid, "(?<=_FA).+$")) %>%
#     mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub)) %>%
#     arrange(lipid_sub) %>%
#     select(lipid, type, lipid_sub, logFC,FDR)
#   
#   
#   plot_title <- paste( Title1,"_vs_", Title2,  sep = "")
#   
#   df_summary <- df_processed %>%
#     group_by(lipid_sub) %>%
#     summarize(
#       Upregulated = sum(logFC > 0),
#       Downregulated = sum(logFC < 0)
#     ) %>%
#     pivot_longer(cols = c("Upregulated", "Downregulated"), names_to = "Regulation", values_to = "Count")
#   filename <- plot_title
#   
#   directory_name <- "TG_DG_FA_lengths"
#   file_name <- paste0("TG_", filename, ".csv")
#   full_path <- file.path(directory_name, file_name)
#   
#   # Create the directory if it doesn't exist
#   if (!dir.exists(directory_name)) {
#     dir.create(directory_name)
#   }
#   
#   # Save the summarized data frame to a CSV file
#   write.csv(df_summary, full_path, row.names = FALSE)
#   
#   # Create the bar plot
#   p<- ggplot(df_summary, aes(x = lipid_sub, y = Count, fill = Regulation)) +
#     geom_bar(stat = "identity", position = position_dodge()) +
#     theme_minimal() +
#     labs(x = "Lipid Subtype", y = "Count", title = plot_title)
#   
#   # Save the plot
#   filename <- paste(dir_name, "/", Title1, Title2, "_TAG_LogFC.pdf", sep = "")
#   ggsave(filename, plot = p, width = 12, height = 8, units = "in")
#   
# 
#   write.csv(df_summary, paste(dir_name, "/", Title1, Title2, "_TAG_.csv", sep = ""), row.names = FALSE)
# }


# plot_histogram_plots_DGs_DE_lipids <- function(df, Title1, Title2) {
#   # Define lipid classes and colors
#   # lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG", "DAG", "TAG | DAG", "DAG | CE", "TAG | DAG | CE")
#   # lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
#   # lipid_colors_alpha <- scales::alpha(lipid_colors, 0.5)
#   # lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
#   # 
#   # # Extract lengths and column indices
#   # len1 <- as.numeric(df$Length1[1])
#   # len2 <- as.numeric(df$Length2[1])
#   # all_cols <- colnames(df)
#   # start_idx <- which(all_cols == "type") + 1
#   # cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   # cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   # blank_col <- all_cols[length(all_cols)]
#   # 
#   # # Adjust and clean data
#   # df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   # df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   # df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   # df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   # 
#   dir_name <- "up_and_down_regulated_FA_tag_DAG_count"
#   # Further processing
#   df_processed <- df %>%
#     # mutate(mean1 = rowMeans(select(., all_of(cols1))),
#     #        mean2 = rowMeans(select(., all_of(cols2)))) %>%
#     filter(type == "DAG") %>%filter(FDR < 0.1) %>%
#     mutate(lipid_sub = str_extract(lipid, "(?<=_C).+$")) %>%
#     mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub)) %>%
#     arrange(lipid_sub) %>%
#     select(lipid, type, lipid_sub, logFC,FDR)
#   
#   write.csv(df_processed, paste(dir_name, "/", Title1, Title2, "All_lipids_DAG_.csv", sep = ""), row.names = FALSE)
#   
#   
#   plot_title <- paste( Title1,"_vs_", Title2,  sep = "")
#   
#   df_summary <- df_processed %>%
#     group_by(lipid_sub) %>%
#     summarize(
#       Upregulated = sum(logFC > 0),
#       Downregulated = sum(logFC < 0)
#     ) %>%
#     pivot_longer(cols = c("Upregulated", "Downregulated"), names_to = "Regulation", values_to = "Count")
#   
#   filename <- plot_title
#   directory_name <- "TG_DG_FA_lengths"
#   file_name <- paste0("DG_", filename, ".csv")
#   full_path <- file.path(directory_name, file_name)
#   
#   # Create the directory if it doesn't exist
#   if (!dir.exists(directory_name)) {
#     dir.create(directory_name)
#   }
#   
#   # Save the summarized data frame to a CSV file
#   write.csv(df_summary, full_path, row.names = FALSE)
#   
#   
#   
#   # Create the bar plot
#   p<- ggplot(df_summary, aes(x = lipid_sub, y = Count, fill = Regulation)) +
#     geom_bar(stat = "identity", position = position_dodge()) +
#     theme_minimal() +
#     labs(x = "Lipid Subtype", y = "Count", title = plot_title)
#   
#   # Save the plot
#   filename <- paste(dir_name, "/", Title1, Title2, "_DAG_LogFC.pdf", sep = "")
#   ggsave(filename, plot = p, width = 12, height = 8, units = "in")
#   
#   
#   write.csv(df_summary, paste(dir_name, "/", Title1, Title2, "_DAG_.csv", sep = ""), row.names = FALSE)
# }

plot_histogram_plots_TGs_DE_lipids <- function(df, Title1, Title2) {
  dir_name <- "up_and_down_regulated_FA_tag_DAG_count"
  
  df_processed <- df %>%
    filter(type == "TAG") %>% filter(FDR < 0.1) %>%
    mutate(lipid_sub = str_extract(lipid, "(?<=_FA).+$")) %>%
    mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub)) %>%
    arrange(lipid_sub) %>%
    select(lipid, type, lipid_sub, logFC, FDR)
  
  
  write.csv(df_processed, paste(dir_name, "/", Title1, Title2, "all_lipids_TAG_.csv", sep = ""), row.names = FALSE)
  
  
  plot_title <- paste(Title1, "_vs_", Title2, sep = "")
  
  df_summary <- df_processed %>%
    group_by(lipid_sub) %>%
    summarize(
      Upregulated = sum(logFC > 0),
      Downregulated = sum(logFC < 0)
    ) %>%
    pivot_longer(cols = c("Upregulated", "Downregulated"), names_to = "Regulation", values_to = "Count") %>%
    mutate(Count = ifelse(Regulation == "Downregulated", -Count, Count))
  
  directory_name <- "TG_DG_FA_lengths"
  file_name <- paste0("TG_", plot_title, ".csv")
  full_path <- file.path(directory_name, file_name)
  
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
  }
  
  write.csv(df_summary, full_path, row.names = FALSE)
  
  p <- ggplot(df_summary, aes(x = lipid_sub, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
    coord_flip() +
    theme_minimal() +
    labs(x = "Lipid Subtype", y = "Count", title = plot_title) +
    scale_y_continuous(labels = abs)
  
  filename <- paste(dir_name, "/", Title1, Title2, "_TAG_LogFC.pdf", sep = "")
  ggsave(filename, plot = p, width = 12, height = 8, units = "in")
  
  write.csv(df_summary, paste(dir_name, "/", Title1, Title2, "_TAG_.csv", sep = ""), row.names = FALSE)
}

plot_histogram_plots_DGs_DE_lipids <- function(df, Title1, Title2) {
  dir_name <- "up_and_down_regulated_FA_tag_DAG_count"
  
  df_processed <- df %>%
    filter(type == "DAG") %>% filter(FDR < 0.1) %>%
    mutate(lipid_sub = str_extract(lipid, "(?<=_C).+$")) %>%
    mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub)) %>%
    arrange(lipid_sub) %>%
    select(lipid, type, lipid_sub, logFC, FDR)
  
  write.csv(df_processed, paste(dir_name, "/", Title1, Title2, "all_lipids_DAG_.csv", sep = ""), row.names = FALSE)
  
  
  
  plot_title <- paste(Title1, "_vs_", Title2, sep = "")
  
  df_summary <- df_processed %>%
    group_by(lipid_sub) %>%
    summarize(
      Upregulated = sum(logFC > 0),
      Downregulated = sum(logFC < 0)
    ) %>%
    pivot_longer(cols = c("Upregulated", "Downregulated"), names_to = "Regulation", values_to = "Count") %>%
    mutate(Count = ifelse(Regulation == "Downregulated", -Count, Count))
  
  directory_name <- "TG_DG_FA_lengths"
  file_name <- paste0("DG_", plot_title, ".csv")
  full_path <- file.path(directory_name, file_name)
  
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
  }
  
  write.csv(df_summary, full_path, row.names = FALSE)
  
  p <- ggplot(df_summary, aes(x = lipid_sub, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
    coord_flip() +
    theme_minimal() +
    labs(x = "Lipid Subtype", y = "Count", title = plot_title) +
    scale_y_continuous(labels = abs)
  
  filename <- paste(dir_name, "/", Title1, Title2, "_DAG_LogFC.pdf", sep = "")
  ggsave(filename, plot = p, width = 12, height = 8, units = "in")
  
  write.csv(df_summary, paste(dir_name, "/", Title1, Title2, "_DAG_.csv", sep = ""), row.names = FALSE)
}



file_list = list.files(path="results", pattern=NULL, all.files=FALSE,
                       full.names=FALSE)
# 
excel_file <- df
df <- read_csv(paste0("results/",file_list[1],sep=""))
jj <-file_list[1]
file_list[2]
for (jj in file_list){
  excel_file <- read_csv(paste0("results/",jj,sep=""))
  dir.create("plots", F)
  jj
  if (!grepl("full", jj, ignore.case = TRUE)) {
    next
  }
  # title_for_plot <- gsub(":", "_",excel_file$Title[1])
  Title1 <- gsub("[:| ]", "_", excel_file$Title_1[1])
  Title2 <- gsub("[:| ]", "_", excel_file$Title_2[1])
  Title1 <- gsub("__", "_", Title1)
  Title2 <- gsub("__", "_", Title2) 
  # Title1 <-
  if (!dir.exists("Stacked_bars")) {
    dir.create("Stacked_bars")
  }
  
  # plot_ridge_plots_TGs_DGs_FA_length(excel_file, Title1, Title2)
  plot_histogram_plots_TGs_DE_lipids(excel_file, Title1, Title2)
  plot_histogram_plots_DGs_DE_lipids(excel_file, Title1, Title2)
  
  
}

