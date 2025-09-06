# Load necessary libraries
library(plyr)
library(ggridges)
library(tidyverse)
library(readxl)
library(edgeR)
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
library(ggforce)
library(gridExtra)
library(patchwork)

# Create necessary directories
if (!dir.exists("FA_lengths_stacked_bar")) {
  dir.create("FA_lengths_stacked_bar")
}
if (!dir.exists("PIE_CHART_FA_lengths")) {
  dir.create("PIE_CHART_FA_lengths")
}
if (!dir.exists("processed_results_2")) {
  dir.create("processed_results_2")
}
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Stacked Bar + Pie Chart Function for FA Chain Lengths
FA_plot_lengths <- function(df, Title1, Title2) {
  # Define FA chain length categories
  # chain_length_classes <- c("Very Short", "Short", "Medium", "Long", "Very Long")
  # chain_length_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99")
  # 
  # categorize_chain_length <- function(fa_type) {
  #   chain_length <- as.numeric(sub(".*\\((\\d+):\\d+\\).*", "\\1", fa_type))
  #   if (is.na(chain_length)) return(NA)
  #   if (chain_length < 6) return("Very Short")
  #   else if (chain_length < 12) return("Short")
  #   else if (chain_length < 22) return("Medium")
  #   else if (chain_length < 28) return("Long")
  #   else return("Very Long")
  # }
  # 
  
  # Define FA chain length categories
  chain_length_classes <- c("Short-Chain", "Medium-Chain", "Long-Chain", "Very Long-Chain", "Ultra Long-Chain")
  # Define color palette for chain length categories
  chain_length_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")
  
  # Function to categorize fatty acid chain lengths based on the number before the colon (e.g., 16 in FA(16:0))
  categorize_chain_length <- function(fa_type) {
    chain_length <- as.numeric(sub(".*\\((\\d+):\\d+\\).*", "\\1", fa_type))
    if (is.na(chain_length)) return(NA)
    if (chain_length < 6) {
      return("Short-Chain")
    } else if (chain_length < 13) {
      return("Medium-Chain")
    } else if (chain_length < 22) {
      return("Long-Chain")
    } else if (chain_length < 30) {
      return("Very Long-Chain")
    } else {
      return("Ultra Long-Chain")
    }
  }
  
  
  # Filter FA lipids
  df <- df %>% filter(type == "FA")
  
  # Extract group columns
  len1 <- as.numeric(df$Length1[1])
  len2 <- as.numeric(df$Length2[1])
  all_cols <- colnames(df)
  start_idx <- which(all_cols == "type") + 1
  cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
  cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
  blank_col <- all_cols[length(all_cols)]
  
  df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
  df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
  
  df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
  df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
  
  df <- df %>%
    mutate(mean1 = rowMeans(select(., one_of(cols1))),
           mean2 = rowMeans(select(., one_of(cols2))),
           chain_length = factor(sapply(lipid, categorize_chain_length),
                                 levels = chain_length_classes))
  
  df_all <- df %>%
    group_by(chain_length) %>%
    summarise(sum_mean1 = sum(mean1),
              sum_mean2 = sum(mean2), .groups = "keep")
  
  df_all_long <- df_all %>%
    select(chain_length, sum_mean1, sum_mean2) %>%
    pivot_longer(cols = starts_with("sum_mean"), names_to = "group", values_to = "value")
  
  # Stacked bar plot
  plot_fa <- ggplot(df_all_long, aes(x = group, y = value, fill = chain_length)) +
    geom_col() +
    theme_classic() +
    ylab("Total lipid content") +
    scale_fill_manual(values = setNames(chain_length_colors, chain_length_classes), name = "Chain Length") +
    scale_x_discrete(labels = c("sum_mean1" = Title1, "sum_mean2" = Title2)) +
    labs(title = paste(Title1, "vs", Title2, "(FA Chain Length)"))
  
  ggsave(filename = paste0("FA_lengths_stacked_bar/", Title1, "_vs_", Title2, "_STACKED_FA_ChainLength.pdf"),
         plot = plot_fa)
  
  # Pie chart function
  pie_plot <- function(data, group_label) {
    ggplot(data, aes(x = "", y = value, fill = chain_length)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0) +
      theme_void() +
      scale_fill_manual(values = setNames(chain_length_colors, chain_length_classes)) +
      labs(title = paste(group_label, "(FA Chain Length)"))
  }
  
  pie1 <- pie_plot(df_all_long %>% filter(group == "sum_mean1"), Title1)
  pie2 <- pie_plot(df_all_long %>% filter(group == "sum_mean2"), Title2)
  
  ggsave(filename = paste0("PIE_CHART_FA_lengths/", Title1, "_PIE_FA_ChainLength.pdf"), plot = pie1)
  ggsave(filename = paste0("PIE_CHART_FA_lengths/", Title2, "_PIE_FA_ChainLength.pdf"), plot = pie2)
}

# Main loop over result files
file_list <- list.files(path = "results", pattern = NULL, full.names = FALSE)

for (jj in file_list) {
  if (!grepl("full", jj, ignore.case = TRUE)) next
  
  excel_file <- read_csv(paste0("results/", jj))
  Title1 <- gsub("[:| ]", "_", excel_file$Title_1[1])
  Title2 <- gsub("[:| ]", "_", excel_file$Title_2[1])
  Title1 <- gsub("__", "_", Title1)
  Title2 <- gsub("__", "_", Title2)
  
  FA_plot_lengths(excel_file, Title1, Title2)
}

