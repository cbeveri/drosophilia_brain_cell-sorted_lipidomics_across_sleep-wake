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
if (!dir.exists("FA_lengths_stacked_bar_select_FA")) {
  dir.create("FA_lengths_stacked_bar_select_FA")
}
if (!dir.exists("PIE_CHART_FA_lengths_select_FA")) {
  dir.create("PIE_CHART_FA_lengths_select_FA")
}
if (!dir.exists("processed_results_2")) {
  dir.create("processed_results_2")
}
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Stacked Bar + Pie Chart Function for FA Chain Lengths
FA_plot_lengths <- function(df, Title1, Title2) {
  # Define FA chain length categories for selected lengths
  chain_length_classes <- c("FA(2:0)", "FA(14:0)", "FA(15:0)", "FA(16:0)")
  # Define color palette for chain length categories
  chain_length_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
  
  # Function to categorize fatty acid chain lengths based on the number before the colon
  categorize_chain_length <- function(fa_type) {
    chain_length <- as.numeric(sub(".*\\((\\d+):\\d+\\).*", "\\1", fa_type))
    if (is.na(chain_length)) return(NA)
    if (chain_length == 2) {
      return("FA(2:0)")
    } else if (chain_length == 14) {
      return("FA(14:0)")
    } else if (chain_length == 15) {
      return("FA(15:0)")
    } else if (chain_length == 16) {
      return("FA(16:0)")
    } else {
      return(NA)
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
           chain_length = sapply(lipid, categorize_chain_length))
  
  # Filter out NA values (chain lengths not in our selected categories)
  df <- df %>% filter(!is.na(chain_length))
  
  # Convert to factor with specified order
  df$chain_length <- factor(df$chain_length, levels = chain_length_classes)
  
  df_all <- df %>%
    group_by(chain_length) %>%
    summarise(sum_mean1 = sum(mean1),
              sum_mean2 = sum(mean2), .groups = "keep")
  
  df_all_long <- df_all %>%
    select(chain_length, sum_mean1, sum_mean2) %>%
    pivot_longer(cols = starts_with("sum_mean"), names_to = "group", values_to = "value")
  
  # Stacked bar plot with adjusted size and spacing
  plot_fa <- ggplot(df_all_long, aes(x = group, y = value, fill = chain_length)) +
    geom_col() +
    theme_classic(base_size = 14) +  # Increased base font size
    ylab("Total lipid content") +
    scale_fill_manual(values = setNames(chain_length_colors, chain_length_classes), name = "Chain Length") +
    scale_x_discrete(labels = c("sum_mean1" = Title1, "sum_mean2" = Title2)) +
    labs(title = paste(Title1, "vs", Title2, "(FA Chain Length)")) +
    theme(
      plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),  # Larger title with more bottom margin
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.2, "cm"),  # Larger legend keys
      plot.margin = margin(20, 20, 20, 20)  # More margin around the plot
    )
  
  ggsave(filename = paste0("FA_lengths_stacked_bar_select_FA/", Title1, "_vs_", Title2, "_STACKED_FA_ChainLength_select_FA.pdf"),
         plot = plot_fa, width = 10, height = 8, units = "in")  # Increased plot dimensions
  
  # Pie chart function with adjusted size
  pie_plot <- function(data, group_label) {
    ggplot(data, aes(x = "", y = value, fill = chain_length)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0) +
      theme_void(base_size = 14) +  # Increased base font size
      scale_fill_manual(values = setNames(chain_length_colors, chain_length_classes)) +
      labs(title = paste(group_label, "(FA Chain Length)")) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1.2, "cm"),
        plot.margin = margin(20, 20, 20, 20)
      )
  }
  
  pie1 <- pie_plot(df_all_long %>% filter(group == "sum_mean1"), Title1)
  pie2 <- pie_plot(df_all_long %>% filter(group == "sum_mean2"), Title2)
  
  ggsave(filename = paste0("PIE_CHART_FA_lengths_select_FA/", Title1, "_PIE_FA_ChainLength_select_FA.pdf"), 
         plot = pie1, width = 10, height = 8, units = "in")  # Increased plot dimensions
  ggsave(filename = paste0("PIE_CHART_FA_lengths_select_FA/", Title2, "_PIE_FA_ChainLength_select_FA.pdf"), 
         plot = pie2, width = 10, height = 8, units = "in")  # Increased plot dimensions
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