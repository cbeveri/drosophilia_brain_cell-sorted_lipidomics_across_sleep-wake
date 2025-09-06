
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

# read the variable from the text file
# cwd <- readLines("Variable_Storage/folder_path.txt")[1]
# cwd
# setwd(cwd)



# setwd(cwd)




file_list = list.files(path="results", pattern=NULL, all.files=FALSE,
                       full.names=FALSE)
# 


df <- excel_file

plot_combined_values_Stacked_with_blank <- function(df, Title1, Title2) {
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  
  lipid_colors_alpha <- alpha(lipid_colors, 0.5)
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
  # df -> excel_file
  # df -> excel_file
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
    filter(FDR < 0.1) %>%
    group_by(type) %>%
    summarise(sum_mean1 = sum(mean1), 
              sum_mean2 = sum(mean2), .groups = "keep")
  
  df_all <- df %>%
    group_by(type) %>%
    summarise(sum_mean1 = sum(mean1), 
              sum_mean2 = sum(mean2), .groups = "keep")
  
  # Combine the two dataframes to have a 'group' column and a 'value' column for plotting
  df_sig_long <- df_sig %>% 
    select(type, sum_mean1, sum_mean2) %>% 
    gather(key = "group", value = "value", -type)
  
  df_all_long <- df_all %>% 
    select(type, sum_mean1, sum_mean2) %>% 
    gather(key = "group", value = "value", -type)
  
  # Plotting significant values
  # plot_sig <- df_sig_long %>%
  #   group_by(type, group) %>%
  #   # mutate(total_intensity = sum(summed_intensity)) %>%
  #   # filter(FDR_40DP_vs_40GFP < 0.10 | FDR_40DP_vs_40DN < 0.10) %>% #filtering for lipids that are FDR<0.10 in DP vs GFP AND/OR DP vs DN 
  #   ungroup() %>%
  #   select(group, type, value) %>%
  #   unique() %>%
  #   ggplot(aes(x=group, y = value, fill = type)) +
  #   geom_col() +
  #   theme_classic() +
  #   ylab("Total lipid content ") +
  #   scale_fill_manual(values = lipid_class_colors, name = "Lipid class", guide = 'none') +
  #   scale_x_discrete(labels = c(Title1, Title2)) + # This line has been added
  #   labs(title = paste(Title1, "vs", Title2, " (All values)"), y = "Sum of Means") 
  
  # Plotting all values
  plot_all <-df_all_long %>%
    group_by(type, group) %>%
    # mutate(total_intensity = sum(summed_intensity)) %>%
    # filter(FDR_40DP_vs_40GFP < 0.10 | FDR_40DP_vs_40DN < 0.10) %>% #filtering for lipids that are FDR<0.10 in DP vs GFP AND/OR DP vs DN 
    ungroup() %>%
    select(group, type, value) %>%
    unique() %>%
    ggplot(aes(x=group, y = value, fill = type)) +
    geom_col() +
    theme_classic() +
    ylab("Total lipid content ") +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class", guide = 'none') +
    scale_x_discrete(labels = c(Title1, Title2)) + # This line has been added
    labs(title = paste(Title1, "vs", Title2, " (All values)"), y = "Sum of Means") 
  # Saving the plots
  # ggsave(filename = paste0("Stacked_bars/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_SUM_sig.svg"), plot = plot_sig)
  ggsave(filename = paste0("Stacked_bars/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_SUM_all.pdf"), plot = plot_all)
}

library(tidyverse)
library(tidyverse)

plot_combined_values_Stacked_with_blank_individual <- function(df, Title1, Title2) {
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  lipid_colors_alpha <- alpha(lipid_colors, 0.5)
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
  
  all_cols <- colnames(df)
  start_idx <- which(all_cols == "type") + 1
  cols1 <- all_cols[start_idx:(start_idx + as.numeric(df$Length1[1]) - 1)]
  cols2 <- all_cols[(start_idx + as.numeric(df$Length1[1])):(start_idx + as.numeric(df$Length1[1]) + as.numeric(df$Length2[1]) - 1)]
  
  blank_col <- all_cols[length(all_cols)]
  
  df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
  df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
  
  all_cols_to_plot <- c(cols1, cols2)
  
  df_long <- df %>%
    pivot_longer(cols = all_of(all_cols_to_plot), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = all_cols_to_plot)) %>%  # Ensuring order
    group_by(type, variable) %>%
    summarise(sum_value = sum(value, na.rm = TRUE), .groups = "keep")
  
  df_long_sig <- df %>%
    filter(FDR < 0.1) %>%
    pivot_longer(cols = all_of(all_cols_to_plot), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = all_cols_to_plot)) %>%  # Ensuring order
    group_by(type, variable) %>%
    summarise(sum_value = sum(value, na.rm = TRUE), .groups = "keep")
  
  plot_all <- df_long %>%
    ggplot(aes(x=variable, y=sum_value, fill=type)) +
    geom_col() +
    theme_classic() +
    ylab("Total lipid content") +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
    labs(title = paste(Title1, "vs", Title2, " (All values)"), y = "Sum of Values") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  plot_sig <- df_long_sig %>%
    ggplot(aes(x=variable, y=sum_value, fill=type)) +
    geom_col() +
    theme_classic() +
    ylab("Total lipid content") +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
    labs(title = paste(Title1, "vs", Title2, " (Significant values)"), y = "Sum of Values") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # Saving the plots
  ggsave(filename = paste0("Stacked_bars_individual/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_all.pdf"), plot = plot_all)
  # ggsave(filename = paste0("Stacked_bars_individual/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_sig.svg"), plot = plot_sig)
}

dir_name <- "processed_results"
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}


plot_pie_charts <- function(df, Title1, Title2, filename, filename_csv) {
  library(gridExtra)
  dir_name <- "processed_results"
  full_path <- file.path(dir_name, filename_csv)
  
  
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  lipid_colors_alpha <- alpha(lipid_colors, 0.5)
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
  
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
  
  
  # write.csv(df, full_path, row.names = FALSE)
  # Compute the means for the groups of columns
  df <- df %>%
    mutate(mean1 = rowMeans(select(., one_of(cols1))),
           mean2 = rowMeans(select(., one_of(cols2))))
  
  # Sum by type for significant values
  df <- df %>%
    #filter(FDR < 0.1) %>%
    group_by(type) %>%
    summarise(sum_mean1 = sum(mean1),
              sum_mean2 = sum(mean2), .groups = "keep")
  
  # write.csv(df_sig, file = filename_csv, row.names = FALSE)
  
  # Make the pie charts
  # pie1 <- df_sig %>%
  #   ggplot(aes(x = "", y = sum_mean1, fill = type)) +
  #   geom_bar(width = 1, stat = "identity", alpha = 0.5, color = "black") +
  #   coord_polar("y", start=0) +
  #   theme_void() + 
  #   theme(axis.text =element_blank(),
  #         axis.line = element_blank()) +
  #   scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
  #   labs(title = Title1)
  # 
  # pie2 <- df_sig %>%
  #   ggplot(aes(x = "", y = sum_mean2, fill = type)) +
  #   geom_bar(width = 1, stat = "identity", alpha = 0.5, color = "black") +
  #   coord_polar("y", start=0) +
  #   theme_void() + 
  #   theme(axis.text =element_blank(),
  #         axis.line = element_blank()) +
  #   scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
  #   labs(title = Title2)
  # Modify pie1
  pie1 <- df %>%
    ggplot(aes(x = "", y = sum_mean1, fill = type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(axis.text =element_blank(),
          axis.line = element_blank()) +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
    labs(title = Title1)
  
  pie2 <- df %>%
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






# 
# 
# plot_ridge_plots_TGs_DGs_FA_length <- function(df, Title1, Title2) {
#   lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG", "DAG", "TAG | DAG", "DAG | CE", "TAG | DAG | CE")
#   lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
#   lipid_colors_alpha <- scales::alpha(lipid_colors, 0.5)
#   lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
#   
#   len1 <- as.numeric(df$Length1[1])
#   len2 <- as.numeric(df$Length2[1])
#   
#   all_cols <- colnames(df)
#   start_idx <- which(all_cols == "type") + 1
#   cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   blank_col <- all_cols[length(all_cols)]
#   
#   df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   
#   df_processed <- df %>%
#     mutate(mean1 = rowMeans(select(., one_of(cols1))),
#            mean2 = rowMeans(select(., one_of(cols2)))) %>%
#     filter(type %in% c("TAG", "DAG")) %>%
#     mutate(lipid_sub = case_when(
#       type == "TAG" ~ str_extract(lipid, "(?<=_FA).+$"),
#       type == "DAG" ~ str_extract(lipid, "(?<=_C).+$"),
#       TRUE ~ NA_character_
#     )) %>%
#     mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub)) %>%
#     arrange(lipid_sub) %>%
#     select(lipid, type, lipid_sub, mean1, mean2)
#   
#   df_long <- df_processed %>%
#     gather(key = "mean_type", value = "intensity", mean1, mean2)
#   
#   # Plotting
#   p <- ggplot(df_long, aes(x = intensity, y = lipid_sub)) +
#     geom_density_ridges() +
#     facet_grid(mean_type ~ type, scales = "free", labeller = labeller(
#       mean_type = c(mean1 = Title1, mean2 = Title2),
#       type = c(TAG = "TAG", DAG = "DAG")
#     )) +
#     labs(x = "Intensity", y = "") +
#     theme_minimal()
#   
#   filename <- paste("FA_distribution/",Title1, Title2, "_TAG_DAG_RidgePlots.pdf", sep = "")
#   ggsave(filename, plot = p)
#   
#   # Save the dataframe
#   # df_to_save <- df %>%
#   #   select(type, lipid, lipid_sub, mean1, mean2) %>%
#   #   spread(key = "mean_type", value = "intensity")
#   
#   write.csv(df_processed, paste("FA_distribution/",Title1, Title2, "_TAG_DAG_RidgePlots.csv", sep = ""), row.names = FALSE)
# }
library(tidyverse)
library(ggridges)

library(tidyverse)
library(ggridges)
# 

library(tidyverse)
library(ggridges)


# Define the function
# plot_histogram_plots_TGs_DGs_FA_length <- function(df, Title1, Title2) {
#   # Define lipid classes and colors
#   lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG", "DAG", "TAG | DAG", "DAG | CE", "TAG | DAG | CE")
#   lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
#   lipid_colors_alpha <- scales::alpha(lipid_colors, 0.5)
#   lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
#   
#   # Extract lengths and column indices
#   len1 <- as.numeric(df$Length1[1])
#   len2 <- as.numeric(df$Length2[1])
#   all_cols <- colnames(df)
#   start_idx <- which(all_cols == "type") + 1
#   cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   blank_col <- all_cols[length(all_cols)]
#   
#   # Adjust and clean data
#   df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   
#   # Further processing
#   df_processed <- df %>%
#     mutate(mean1 = rowMeans(select(., one_of(cols1))),
#            mean2 = rowMeans(select(., one_of(cols2)))) %>%
#     filter(type %in% c("TAG", "DAG")) %>%
#     mutate(lipid_sub = case_when(
#       type == "TAG" ~ str_extract(lipid, "(?<=_FA).+$"),
#       type == "DAG" ~ str_extract(lipid, "(?<=_C).+$"),
#       TRUE ~ NA_character_
#     )) %>%
#     mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub)) %>%
#     arrange(lipid_sub) %>%
#     select(lipid, type, lipid_sub, mean1, mean2)
#   
#   # Convert to long format
#   df_long <- df_processed %>%
#     gather(key = "mean_type", value = "intensity", mean1, mean2)
#   
#   # Extract top 10% intensity values and count the number of occurrences in each lipid_sub
#   # Extract top 10% intensity values and count the number of occurrences in each lipid_sub
#   df_top_10 <- df_long %>%
#     group_by(lipid_sub, mean_type, type) %>%
#     arrange(desc(intensity)) %>%
#     mutate(rank = row_number()) %>%
#     filter(rank <= ceiling(0.1 * n())) %>%
#     summarise(count = n())
#   # Create separate plots
#   plots <- list()
#   for (mean_t in c("mean1", "mean2")) {
#     for (typ in c("TAG", "DAG")) {
#       p <- ggplot(df_top_10 %>% filter(mean_type == mean_t, type == typ), 
#                   aes(x = lipid_sub, y = count, fill = lipid_sub)) +
#         geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
#         labs(x = "Lipid Sub", y = "Count", 
#              title = ifelse(mean_t == "mean1", Title1, Title2),
#              subtitle = typ) +
#         theme_minimal() +
#         scale_fill_manual(values = lipid_class_colors)
#       
#       plots[[paste(mean_t, typ)]] <- p
#     }
#   }
#   
#   # Combine plots using patchwork
#   combined_plot <- plots[["mean1_TAG"]] + plots[["mean2_TAG"]] + plots[["mean1_DAG"]] + plots[["mean2_DAG"]] +
#     plot_layout(ncol = 2)
#   
#   # Ensure the directory exists
#   dir_name <- "FA_distribution_histogram"
#   if (!dir.exists(dir_name)) {
#     dir.create(dir_name)
#   }
#   
#   # Save the combined plot
#   filename <- paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms.pdf", sep = "")
#   ggsave(filename, plot = combined_plot)
#   
#   # Save the dataframe
#   write.csv(df_top_10, paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms.csv", sep = ""), row.names = FALSE)
# }
library(tidyverse)
library(ggridges)
# 
# plot_histogram_plots_TGs_DGs_FA_length_2 <- function(df, Title1, Title2) {
#   # Define lipid classes and colors
#   lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG", "DAG", "TAG | DAG", "DAG | CE", "TAG | DAG | CE")
#   lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
#   lipid_colors_alpha <- scales::alpha(lipid_colors, 0.5)
#   lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
# 
#   # Extract lengths and column indices
#   len1 <- as.numeric(df$Length1[1])
#   len2 <- as.numeric(df$Length2[1])
#   all_cols <- colnames(df)
#   start_idx <- which(all_cols == "type") + 1
#   cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   blank_col <- all_cols[length(all_cols)]
# 
#   # Adjust and clean data
#   df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
# 
#   # Further processing
#   df_processed <- df %>%
#     mutate(mean1 = rowMeans(select(., one_of(cols1))),
#            mean2 = rowMeans(select(., one_of(cols2)))) %>%
#     filter(type == "TAG") %>%
#     mutate(lipid_sub = str_extract(lipid, "(?<=_FA).+$")) %>%
#     mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub)) %>%
#     arrange(lipid_sub) %>%
#     select(lipid, type, lipid_sub, mean1, mean2)
# 
#   # Convert to long format
#   df_long <- df_processed %>%
#     gather(key = "mean_type", value = "intensity", mean1, mean2)
# 
#   # Extract top 10% intensity values and count the number of occurrences in each lipid_sub
#   df_top_10 <- df_long %>%
#     group_by(lipid_sub, mean_type, type) %>%
#     mutate(top_10_threshold = quantile(intensity, 0.9)) %>%
#     filter(intensity >= top_10_threshold) %>%
#     summarise(count = n())
# 
#   # Plotting: Single plot with 4 facets
#   p <- ggplot(df_top_10, aes(x = lipid_sub, y = count, fill = lipid_sub)) +
#     geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
#     facet_grid(mean_type ~ type, scales = "free", labeller = labeller(
#       mean_type = c(mean1 = Title1, mean2 = Title2),
#       type = c(TAG = "TAG")
#     )) +
#     labs(x = "Lipid Sub", y = "Count") +
#     theme_minimal() +
#     scale_fill_manual(values = lipid_class_colors) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
#   # Ensure the directory exists
#   dir_name <- "FA_distribution_histogram"
#   if (!dir.exists(dir_name)) {
#     dir.create(dir_name)
#   }
# 
#   # Save the plot
#   filename <- paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms.pdf", sep = "")
#   ggsave(filename, plot = p)
# 
#   # Save the dataframe
#   write.csv(df_top_10, paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms.csv", sep = ""), row.names = FALSE)
# }


##Works but does not label.
# plot_histogram_plots_TGs_DGs_FA_length_3 <- function(df, Title1, Title2) {
#   # Define lipid classes and colors
#   lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG", "DAG", "TAG | DAG", "DAG | CE", "TAG | DAG | CE")
#   lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
#   lipid_colors_alpha <- scales::alpha(lipid_colors, 0.5)
#   lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
#   
#   # Extract lengths and column indices
#   len1 <- as.numeric(df$Length1[1])
#   len2 <- as.numeric(df$Length2[1])
#   all_cols <- colnames(df)
#   start_idx <- which(all_cols == "type") + 1
#   cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   blank_col <- all_cols[length(all_cols)]
#   
#   # Adjust and clean data
#   df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   
#   # Further processing
#   df_processed <- df %>%
#     mutate(mean1 = rowMeans(select(., all_of(cols1))),
#            mean2 = rowMeans(select(., all_of(cols2)))) %>%
#     filter(type == "TAG") %>%
#     mutate(lipid_sub = str_extract(lipid, "(?<=_FA).+$")) %>%
#     mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub)) %>%
#     arrange(lipid_sub) %>%
#     select(lipid, type, lipid_sub, mean1, mean2)
#   
#   # Helper function to calculate top 10% for a given mean
#   top_10_percentiles <- function(data, column) {
#     threshold <- quantile(data[[column]], 0.9)
#     data %>% filter(data[[column]] >= threshold)
#   }
#   
#   # Calculate top 10% separately for mean1 and mean2
#   df_mean1_top10 <- top_10_percentiles(df_processed, "mean1") %>% 
#     count(lipid_sub, type) %>% 
#     mutate(mean_type = "mean1")
#   
#   df_mean2_top10 <- top_10_percentiles(df_processed, "mean2") %>% 
#     count(lipid_sub, type) %>% 
#     mutate(mean_type = "mean2")
#   
#   # Combine for plotting
#   df_plot <- rbind(df_mean1_top10, df_mean2_top10)
#   
#   # Plotting
#   p <- ggplot(df_plot, aes(x = lipid_sub, y = n, fill = mean_type)) +
#     geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.7) +
#     facet_wrap(~ type, scales = "free") +
#     labs(x = "Lipid Sub", y = "Count", fill = "Mean Type") +
#     scale_fill_manual(values = c("mean1" = "#ff9999", "mean2" = "#9999ff")) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # Ensure the directory exists
#   dir_name <- "FA_distribution_histogram"
#   if (!dir.exists(dir_name)) {
#     dir.create(dir_name)
#   }
#   
#   # Save the plot
#   filename <- paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms.pdf", sep = "")
#   ggsave(filename, plot = p, width = 12, height = 8, units = "in")
#   
#   # Save the dataframes
#   write.csv(df_mean1_top10, paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms_mean1.csv", sep = ""), row.names = FALSE)
#   write.csv(df_mean2_top10, paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms_mean2.csv", sep = ""), row.names = FALSE)
# }



# 
# 
# 
# ##Does not plot legends but does work
# plot_histogram_plots_TGs_DGs_FA_length_3 <- function(df, Title1, Title2) {
#   # Define lipid classes and colors
#   lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG", "DAG", "TAG | DAG", "DAG | CE", "TAG | DAG | CE")
#   lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
#   lipid_colors_alpha <- scales::alpha(lipid_colors, 0.5)
#   lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
#   
#   # Extract lengths and column indices
#   len1 <- as.numeric(df$Length1[1])
#   len2 <- as.numeric(df$Length2[1])
#   all_cols <- colnames(df)
#   start_idx <- which(all_cols == "type") + 1
#   cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   blank_col <- all_cols[length(all_cols)]
#   
#   # Adjust and clean data
#   df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   
#   # Further processing
#   df_processed <- df %>%
#     mutate(mean1 = rowMeans(select(., all_of(cols1))),
#            mean2 = rowMeans(select(., all_of(cols2)))) %>%
#     filter(type == "TAG") %>%
#     mutate(lipid_sub = str_extract(lipid, "(?<=_FA).+$")) %>%
#     mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub)) %>%
#     arrange(lipid_sub) %>%
#     select(lipid, type, lipid_sub, mean1, mean2)
#   
#   # Helper function to calculate top 10% for a given mean
#   top_10_percentiles <- function(data, column) {
#     threshold <- quantile(data[[column]], 0.9)
#     data %>% filter(data[[column]] >= threshold)
#   }
#   
#   # Calculate top 10% separately for mean1 and mean2
#   df_mean1_top10 <- top_10_percentiles(df_processed, "mean1") %>% 
#     count(lipid_sub, type) %>% 
#     mutate(mean_type = Title1)
#   
#   df_mean2_top10 <- top_10_percentiles(df_processed, "mean2") %>% 
#     count(lipid_sub, type) %>% 
#     mutate(mean_type = Title2)
#   
#   # Combine for plotting
#   df_plot <- rbind(df_mean1_top10, df_mean2_top10)
#   
#   # Plotting
#   p <- ggplot(df_plot, aes(x = lipid_sub, y = n, fill = mean_type)) +
#     geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.7) +
#     facet_wrap(~ type, scales = "free") +
#     labs(x = "Lipid Sub", y = "Count", fill = "Mean Type") +
#     scale_fill_manual(values = c(Title1 = "#ff9999", Title2 = "#9999ff")) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # Ensure the directory exists
#   dir_name <- "FA_distribution_histogram"
#   if (!dir.exists(dir_name)) {
#     dir.create(dir_name)
#   }
#   
#   # Save the plot
#   filename <- paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms.pdf", sep = "")
#   ggsave(filename, plot = p, width = 12, height = 8, units = "in")
#   
#   # Save the combined dataframes as a single CSV
#   combined_df <- rbind(df_mean1_top10, df_mean2_top10)
#   write.csv(combined_df, paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms_combined.csv", sep = ""), row.names = FALSE)
# }
plot_histogram_plots_TGs_DGs_FA_length_3 <- function(df, Title1, Title2) {
  # Define lipid classes and colors
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG", "DAG", "TAG | DAG", "DAG | CE", "TAG | DAG | CE")
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
  lipid_colors_alpha <- scales::alpha(lipid_colors, 0.5)
  lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
  
  # Extract lengths and column indices
  len1 <- as.numeric(df$Length1[1])
  len2 <- as.numeric(df$Length2[1])
  all_cols <- colnames(df)
  start_idx <- which(all_cols == "type") + 1
  cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
  cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
  blank_col <- all_cols[length(all_cols)]
  
  # Adjust and clean data
  df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
  df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
  df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
  df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
  
  # Further processing
  df_processed <- df %>%
    mutate(mean1 = rowMeans(select(., all_of(cols1))),
           mean2 = rowMeans(select(., all_of(cols2)))) %>%
    filter(type == "TAG") %>%
    mutate(lipid_sub = str_extract(lipid, "(?<=_FA).+$")) %>%
    mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub)) %>%
    arrange(lipid_sub) %>%
    select(lipid, type, lipid_sub, mean1, mean2)
  
  # Helper function to calculate top 10% for a given mean
  top_10_percentiles <- function(data, column, mean_label) {
    data %>%
      filter(data[[column]] >= quantile(data[[column]], 0.9)) %>%
      count(lipid_sub, type) %>%
      mutate(mean_type = mean_label)
  }
  
  # Calculate top 10% separately for mean1 and mean2
  df_mean1_top10 <- top_10_percentiles(df_processed, "mean1", Title1)
  df_mean2_top10 <- top_10_percentiles(df_processed, "mean2", Title2)
  
  # Combine for plotting
  df_plot <- rbind(df_mean1_top10, df_mean2_top10)
  
  # Plotting
  p <- ggplot(df_plot, aes(x = lipid_sub, y = n, fill = mean_type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.7) +
    facet_wrap(~ type, scales = "free") +
    labs(x = "Lipid Sub", y = "Count", fill = "Mean Type") +
    scale_fill_manual(values = setNames(c("#ff9999", "#9999ff"), c(Title1, Title2))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Ensure the directory exists
  dir_name <- "FA_distribution_histogram"
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  # Save the plot
  filename <- paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms.pdf", sep = "")
  ggsave(filename, plot = p, width = 12, height = 8, units = "in")
  
  # Save the combined dataframes as a single CSV
  combined_df <- rbind(df_mean1_top10, df_mean2_top10)
  write.csv(combined_df, paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms_combined.csv", sep = ""), row.names = FALSE)
}
dir_name <- "processed_results_2"
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}


proccess_results <- function(df, Title1, Title2, filename, filename_csv) {
  library(gridExtra)
  dir_name <- "processed_results_2"
  full_path <- file.path(dir_name, filename_csv)
  
  
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  lipid_colors_alpha <- alpha(lipid_colors, 0.5)
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
  
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
  
  
  # write.csv(df, full_path, row.names = FALSE)
  # Compute the means for the groups of columns
  df <- df %>%#filter(FDR < 0.1) %>%
    mutate(mean1 = rowMeans(select(., one_of(cols1))),
           mean2 = rowMeans(select(., one_of(cols2))))
  write.csv(df, full_path, row.names = FALSE)
  
}
df
library(patchwork)

jj
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
  # plot_histogram_plots_TGs_DGs_FA_length_3(excel_file, Title1, Title2)
  title_for_plot <- paste0(Title1,Title2,sep="_")
  # plot_combined_values_Stacked(excel_file, Title1, Title2)
  plot_combined_values_Stacked_with_blank(excel_file, Title1, Title2)
  plot_combined_values_Stacked_with_blank_individual(excel_file, Title1, Title2)
  plot_pie_charts(excel_file, Title1, Title2, paste0("Pie_charts/",Title1,Title2,".svg",sep=""), paste0("Pie_charts/",Title1,Title2,".csv",sep=""))
  proccess_results(excel_file, Title1, Title2, paste0("Pie_charts/",Title1,Title2,".svg",sep=""), paste0(Title1,Title2,".csv",sep=""))
  
  # Plotting and saving
  # plot_object <- plot_significant_lipids(excel_file, title_for_plot)
  # process_and_plot(excel_file, Title1, Title2)
  # Ensuring the plots directory exists

  

}
  

  
   
  # 
  # make_heatmap <- function(tp, design_mat, gr, contrasts, title, file) {
  #   
  #   DElist <-
  #     tp %>%
  #     get_DE_lipids(design_mat, gr, contrasts)
  #   
  #   if(length(DElist) == 0) {
  #     cat("No significant lipids for ", title)
  #     return()
  #   }
  #   
  #   if(length(DElist) == 1) {
  #     cat("Single significant lipids for ", title, " is ", DElist[1])
  #     return()
  #   }
  #   
  #   Blank <- log2(tp[[blank_name]])
  #   
  #   tp %>%
  #     mutate(lipid = make.unique(lipid)) %>%
  #     filter(lipid %in% DElist) %>%
  #     select( -type) %>%
  #     mutate_if(is.numeric, log2) %>%
  #     mutate_if(is.numeric, list(~ . - Blank)) %>%
  #     select(- blank_name) %>% 
  #     column_to_rownames("lipid") %>%
  #     as.matrix() %>%
  #     pheatmap::pheatmap(main = title,cluster_cols = none,
  #                        cluster_rows = none, filename = file,
  #                        cellheight = 10)
  # }
  # 
  # cells_lipid_expr %>%
  #   make_heatmap(design_expr, gr_expr, contrasts_expr, title_for_plot, 
  #                file = paste("plots/Heatmap_",title_for_plot,".pdf",sep=''))
# }




