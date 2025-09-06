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
  lipid_classes <- c("DAG", "TAG")
  lipid_colors <- c( "#8dd3c7", "#6a3d9a")
  
  lipid_colors_alpha <- alpha(lipid_colors, 0.5)
  lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
  
  # Filter the dataframe to only include FA, DAG, TAG
  df <- df %>% filter(type %in% lipid_classes)
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
  
  
  
  df_all <- df %>%
    group_by(type) %>%
    summarise(sum_mean1 = sum(mean1), 
              sum_mean2 = sum(mean2), .groups = "keep")
  
  
  df_all_long <- df_all %>% 
    select(type, sum_mean1, sum_mean2) %>% 
    gather(key = "group", value = "value", -type) %>%
    group_by(group) %>%
    mutate(percent = 100 * value / sum(value)) %>%
    ungroup()
  
  # Plotting all values
  plot_all <-df_all_long %>%
    group_by(type, group) %>%
    # mutate(total_intensity = sum(summed_intensity)) %>%
    # filter(FDR_40DP_vs_40GFP < 0.10 | FDR_40DP_vs_40DN < 0.10) %>% #filtering for lipids that are FDR<0.10 in DP vs GFP AND/OR DP vs DN 
    ungroup() %>%
    select(group, type, value, percent) %>%
    unique() %>%
    ggplot(aes(x=group, y = value, fill = type)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", percent)),
              position = position_stack(vjust = 0.5),
              size = 3, colour = "black") +
    theme_classic() +
    ylab("Total lipid content ") +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class", guide = 'none') +
    scale_x_discrete(labels = c(Title1, Title2)) + # This line has been added
    labs(title = paste(Title1, "vs", Title2, " (All values)"), y = "Sum of Means") 
  # Saving the plots
  # ggsave(filename = paste0("Stacked_bars_DG_TG_only/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_SUM_sig.svg"), plot = plot_sig)
  ggsave(filename = paste0("Stacked_bars_DG_TG_only/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_SUM_all.pdf"), plot = plot_all)
}

library(tidyverse)
library(tidyverse)

plot_combined_values_Stacked_with_blank_individual <- function(df, Title1, Title2) {
  lipid_classes <- c("DAG", "TAG")
  lipid_colors <- c( "#8dd3c7", "#6a3d9a")
  
  lipid_colors_alpha <- alpha(lipid_colors, 0.5)
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
  ggsave(filename = paste0("Stacked_bars_DG_TG_only_individual/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_all.pdf"), plot = plot_all)
  # ggsave(filename = paste0("Stacked_bars_DG_TG_only_individual/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_sig.svg"), plot = plot_sig)
}

dir_name <- "processed_results"
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}



library(tidyverse)
library(ggridges)

library(tidyverse)
library(ggridges)
# 

library(tidyverse)
library(ggridges)



library(tidyverse)
library(ggridges)
# 


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
  if (!dir.exists("Stacked_bars_DG_TG_only")) {
    dir.create("Stacked_bars_DG_TG_only")
  }
  
  title_for_plot <- paste0(Title1,Title2,sep="_")
  plot_combined_values_Stacked_with_blank(excel_file, Title1, Title2)
  plot_combined_values_Stacked_with_blank_individual(excel_file, Title1, Title2)
  
  
  
  
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



# library(plyr)
# library(ggridges)
# library(tidyverse)
# library(readxl)
# library(edgeR)
# library(tidyverse)
# library(factoextra)
# library(grid)
# library(ggsignif)
# library(stringr)
# library(reshape2)
# library(gtools)
# library(tidyr)
# library(ggbeeswarm) #https://github.com/eclarke/ggbeeswarm
# library(ggrepel)
# library(scales)
# 
# library(ggplot2)
# #library(ggExtra) #https://cran.r-project.org/web/packages/ggExtra/vignettes/ggExtra.html
# 
# library(cowplot)
# library(ggridges)
# 
# # library(limma)
# # library(writexl)
# 
# 
# # install.packages("ggforce")
# library(ggforce)
# getwd()
# 
# # read the variable from the text file
# # cwd <- readLines("Variable_Storage/folder_path.txt")[1]
# # cwd
# # setwd(cwd)
# 
# 
# 
# # setwd(cwd)
# 
# 
# 
# 
# file_list = list.files(path="results", pattern=NULL, all.files=FALSE,
#                        full.names=FALSE)
# # 
# 
# 
# df <- excel_file
# 
# plot_combined_values_Stacked_with_blank <- function(df, Title1, Title2) {
#   lipid_classes <- c("DAG", "TAG")
#   lipid_colors <- c( "#8dd3c7", "#6a3d9a")
#   
#   lipid_colors_alpha <- alpha(lipid_colors, 0.5)
#   lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
#   
#   # Filter the dataframe to only include FA, DAG, TAG
#   df <- df %>% filter(type %in% lipid_classes)
#   # df -> excel_file
#   # Extract column lengths
#   len1 <- as.numeric(df$Length1[1])
#   len2 <- as.numeric(df$Length2[1])
#   
#   # Columns after 'type' are the ones of interest
#   all_cols <- colnames(df)
#   start_idx <- which(all_cols == "type") + 1
#   cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   
#   
#   # Get the blank column (which is the penultimate column)
#   blank_col <- all_cols[length(all_cols)]
#   
#   # Subtract the blank column from the columns of interest and set negative values to 0
#   df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
# 
#   df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   
#   
#   # Compute the means for the groups of columns
#   df <- df %>%
#     mutate(mean1 = rowMeans(select(., one_of(cols1))),
#            mean2 = rowMeans(select(., one_of(cols2))))
#   
# 
#   
#   df_all <- df %>%
#     group_by(type) %>%
#     summarise(sum_mean1 = sum(mean1), 
#               sum_mean2 = sum(mean2), .groups = "keep")
#   
# 
#   df_all_long <- df_all %>% 
#     select(type, sum_mean1, sum_mean2) %>% 
#     gather(key = "group", value = "value", -type)
# 
#   # Plotting all values
#   plot_all <-df_all_long %>%
#     group_by(type, group) %>%
#     # mutate(total_intensity = sum(summed_intensity)) %>%
#     # filter(FDR_40DP_vs_40GFP < 0.10 | FDR_40DP_vs_40DN < 0.10) %>% #filtering for lipids that are FDR<0.10 in DP vs GFP AND/OR DP vs DN 
#     ungroup() %>%
#     select(group, type, value) %>%
#     unique() %>%
#     ggplot(aes(x=group, y = value, fill = type)) +
#     geom_col() +
#     theme_classic() +
#     ylab("Total lipid content ") +
#     scale_fill_manual(values = lipid_class_colors, name = "Lipid class", guide = 'none') +
#     scale_x_discrete(labels = c(Title1, Title2)) + # This line has been added
#     labs(title = paste(Title1, "vs", Title2, " (All values)"), y = "Sum of Means") 
#   # Saving the plots
#   # ggsave(filename = paste0("Stacked_bars_DG_TG_only/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_SUM_sig.svg"), plot = plot_sig)
#   ggsave(filename = paste0("Stacked_bars_DG_TG_only/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_SUM_all.pdf"), plot = plot_all)
# }
# 
# library(tidyverse)
# library(tidyverse)
# 
# plot_combined_values_Stacked_with_blank_individual <- function(df, Title1, Title2) {
#   lipid_classes <- c("DAG", "TAG")
#   lipid_colors <- c( "#8dd3c7", "#6a3d9a")
#   
#   lipid_colors_alpha <- alpha(lipid_colors, 0.5)
#   lipid_class_colors <- setNames(lipid_colors_alpha, lipid_classes)
#   
#   all_cols <- colnames(df)
#   start_idx <- which(all_cols == "type") + 1
#   cols1 <- all_cols[start_idx:(start_idx + as.numeric(df$Length1[1]) - 1)]
#   cols2 <- all_cols[(start_idx + as.numeric(df$Length1[1])):(start_idx + as.numeric(df$Length1[1]) + as.numeric(df$Length2[1]) - 1)]
#   
#   blank_col <- all_cols[length(all_cols)]
#   
#   df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   
#   all_cols_to_plot <- c(cols1, cols2)
#   
#   df_long <- df %>%
#     pivot_longer(cols = all_of(all_cols_to_plot), names_to = "variable", values_to = "value") %>%
#     mutate(variable = factor(variable, levels = all_cols_to_plot)) %>%  # Ensuring order
#     group_by(type, variable) %>%
#     summarise(sum_value = sum(value, na.rm = TRUE), .groups = "keep")
#   
#   df_long_sig <- df %>%
#     filter(FDR < 0.1) %>%
#     pivot_longer(cols = all_of(all_cols_to_plot), names_to = "variable", values_to = "value") %>%
#     mutate(variable = factor(variable, levels = all_cols_to_plot)) %>%  # Ensuring order
#     group_by(type, variable) %>%
#     summarise(sum_value = sum(value, na.rm = TRUE), .groups = "keep")
#   
#   plot_all <- df_long %>%
#     ggplot(aes(x=variable, y=sum_value, fill=type)) +
#     geom_col() +
#     theme_classic() +
#     ylab("Total lipid content") +
#     scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
#     labs(title = paste(Title1, "vs", Title2, " (All values)"), y = "Sum of Values") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#   
#   plot_sig <- df_long_sig %>%
#     ggplot(aes(x=variable, y=sum_value, fill=type)) +
#     geom_col() +
#     theme_classic() +
#     ylab("Total lipid content") +
#     scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
#     labs(title = paste(Title1, "vs", Title2, " (Significant values)"), y = "Sum of Values") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#   
#   # Saving the plots
#   ggsave(filename = paste0("Stacked_bars_DG_TG_only_individual/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_all.pdf"), plot = plot_all)
#   # ggsave(filename = paste0("Stacked_bars_DG_TG_only_individual/", Title1, "_vs_", Title2, "_STACKED__with_BLANK_sig.svg"), plot = plot_sig)
# }
# 
# dir_name <- "processed_results"
# if (!dir.exists(dir_name)) {
#   dir.create(dir_name)
# }
# 
# 
# 
# library(tidyverse)
# library(ggridges)
# 
# library(tidyverse)
# library(ggridges)
# # 
# 
# library(tidyverse)
# library(ggridges)
# 
# 
# 
# library(tidyverse)
# library(ggridges)
# # 
# 
# 
# library(patchwork)
# 
# jj
# df <- read_csv(paste0("results/",file_list[1],sep=""))
# jj <-file_list[1]
# file_list[2]
# for (jj in file_list){
#   excel_file <- read_csv(paste0("results/",jj,sep=""))
#   dir.create("plots", F)
#   jj
#   if (!grepl("full", jj, ignore.case = TRUE)) {
#     next
#   }
#   # title_for_plot <- gsub(":", "_",excel_file$Title[1])
#   Title1 <- gsub("[:| ]", "_", excel_file$Title_1[1])
#   Title2 <- gsub("[:| ]", "_", excel_file$Title_2[1])
#   Title1 <- gsub("__", "_", Title1)
#   Title2 <- gsub("__", "_", Title2) 
#   # Title1 <-
#   if (!dir.exists("Stacked_bars_DG_TG_only")) {
#     dir.create("Stacked_bars_DG_TG_only")
#   }
# 
#   title_for_plot <- paste0(Title1,Title2,sep="_")
#   plot_combined_values_Stacked_with_blank(excel_file, Title1, Title2)
#   plot_combined_values_Stacked_with_blank_individual(excel_file, Title1, Title2)
# 
# 
#   
# 
# }
#   
# 
#   
#    
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




