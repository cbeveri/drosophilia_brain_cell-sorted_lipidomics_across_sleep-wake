
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



plot_scatter_blank_subtraction <- function(df, Title1, Title2) {
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  
  
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors, lipid_classes)
  
  
  # Extract column lengths
  len1 <- as.numeric(df$Length1[1])
  len2 <- as.numeric(df$Length2[1])
  
  # Columns after 'type' are the ones of interest
  all_cols <- colnames(df)
  # Exclude the last column
  all_cols <- colnames(df)[-length(df)]
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
  
  intensity_columns <- c(cols1,cols2)
  # Remove rows where all values in intensity_columns are zeros
  df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
  
  # Compute the means for the groups of columns
  df <- df %>%
    mutate(mean1 = rowMeans(select(., one_of(cols1))),
           mean2 = rowMeans(select(., one_of(cols2))))
  df_sig <- df %>%
    filter(FDR < 0.1) 


  plot_all <- df %>%
    ggplot(aes(x=log10(mean2), y=log10(mean1))) +
    # First, plot all points in grey
    geom_point(aes(color=NULL), color="grey", size=3) +
    
    # Then, plot significant lipids on top, colored by logFC
    geom_point(data=df_sig, 
               aes(color=logFC),
               size=3) +
    
    # Add labels for the first 10 significant lipids
    ggrepel::geom_text_repel(data=head(df_sig, 20),
                             aes(label=lipid),  # Assuming 'lipid' is the column with lipid names
                             size=3, 
                             nudge_y = 0.2,    # Adjust as needed
                             nudge_x = 0.2) +
    
    # Color scale for logFC with squishing out-of-bounds values
    scale_color_gradient2(midpoint=0, low="Blue", mid="#DC9313", 
                          high="Red", space="Lab", 
                          limits=c(-1, 1), breaks=c(-1, -0.5, 0, 0.5, 1), 
                          oob=squish) +
    
    # Themes and title
    theme_bw() +
    labs(
      x = paste("log10(", Title2, ")", sep=""),
      y = paste("log10(", Title1, ")", sep=""),
      title = paste(Title1, " vs ", Title2)
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_all
  
  
  ggsave(filename = paste0("scatter_plots/", Title1, "_vs_", Title2, "_scatter_plot_blank_subtraction.pdf"), plot = plot_all)
  
  ggsave(filename = paste0("scatter_plots/", Title1, "_vs_", Title2, "_scatter_plot_blank_subtraction.svg"), plot = plot_all)
}


plot_scatter_by_type_blank_subtraction <- function(df, Title1, Title2) {
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors, lipid_classes)
  
  # Extract column lengths
  len1 <- as.numeric(df$Length1[1])
  len2 <- as.numeric(df$Length2[1])
  
  # Columns after 'type' are the ones of interest
  all_cols <- colnames(df)
  # Exclude the last column
  all_cols <- colnames(df)[-length(df)]
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
  
  intensity_columns <- c(cols1,cols2)
  # Remove rows where all values in intensity_columns are zeros
  df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
  
  # Compute the means for the groups of columns
  df <- df %>%
    mutate(mean1 = rowMeans(select(., one_of(cols1))),
           mean2 = rowMeans(select(., one_of(cols2))))
  
  # Add transformed variables to the dataframe
  df$x <- ifelse(df$mean2 == 0, 0, log10(df$mean2))
  df$y <- ifelse(df$mean1 == 0, 0, log10(df$mean1))
  df_sig <- df %>%
    filter(FDR < 0.1) 
  
  # Calculate the axis limits
  # x_range <- range(log10(df$mean2), na.rm=TRUE)
  # y_range <- range(log10(df$mean1), na.rm=TRUE)
  # 
  # extended_x_limits <- c(x_range[1] - 0.5, x_range[2] + 0.5)
  # extended_y_limits <- c(y_range[1] - 0.5, y_range[2] + 0.5)
  
  plot_all <- df %>%
    ggplot(aes(x=x, y=y)) +
    # First, plot all points in grey
    # xlim(extended_x_limits) +
    # ylim(extended_y_limits)+
    geom_point(aes(color=NULL), color="grey", size=3) +
    
    # Then, plot significant lipids on top, colored by lipid type
    geom_point(data=df_sig, 
               aes(color=type),
               size=3, alpha=0.5) +
    
    # Add labels for the first 10 significant lipids
    ggrepel::geom_text_repel(data=head(df_sig, 20),
                             aes(label=lipid),  # Assuming 'lipid' is the column with lipid names
                             size=3, 
                             nudge_y = 0.2,    # Adjust as needed
                             nudge_x = 0.2) +
    
    # Color scale for lipid types using lipid_class_colors
    scale_color_manual(values=lipid_class_colors) +
    
    # Themes and title
    theme_bw() +
    labs(
      x = paste("log10(", Title2, ")", sep=""),
      y = paste("log10(", Title1, ")", sep=""),
      title = paste(Title1, " vs ", Title2)
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_all
  
  
  ggsave(filename = paste0("scatter_plots/", Title1, "_vs_", Title2, "_scatter_plot_by_type_blank_subtraction.pdf"), plot = plot_all)
  
  ggsave(filename = paste0("scatter_plots/", Title1, "_vs_", Title2, "_scatter_plot_by_type_blank_subtraction.png"), plot = plot_all)
}



plot_scatter_by_type_blank_subtraction_no_Log <- function(df, Title1, Title2) {
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors, lipid_classes)
  
  # Extract column lengths
  len1 <- as.numeric(df$Length1[1])
  len2 <- as.numeric(df$Length2[1])
  
  # Columns after 'type' are the ones of interest
  all_cols <- colnames(df)
  # Exclude the last column
  all_cols <- colnames(df)[-length(df)]
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
  
  intensity_columns <- c(cols1,cols2)
  # Remove rows where all values in intensity_columns are zeros
  df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
  
  # Compute the means for the groups of columns
  df <- df %>%
    mutate(mean1 = rowMeans(select(., one_of(cols1))),
           mean2 = rowMeans(select(., one_of(cols2))))
  
  
  ##Removed this line
  # Add transformed variables to the dataframe
  # df$x <- ifelse(df$mean2 == 0, 0, log10(df$mean2))
  # df$y <- ifelse(df$mean1 == 0, 0, log10(df$mean1))
  
  df$x <- df$mean2
  df$y <- df$mean1
  
  df_sig <- df %>%
    filter(FDR < 0.1) 
  
  # Calculate the axis limits
  # x_range <- range(log10(df$mean2), na.rm=TRUE)
  # y_range <- range(log10(df$mean1), na.rm=TRUE)
  # 
  # extended_x_limits <- c(x_range[1] - 0.5, x_range[2] + 0.5)
  # extended_y_limits <- c(y_range[1] - 0.5, y_range[2] + 0.5)
  
  plot_all <- df %>%
    ggplot(aes(x=x, y=y)) +
    # First, plot all points in grey
    # xlim(extended_x_limits) +
    # ylim(extended_y_limits)+
    geom_point(aes(color=NULL), color="grey", size=3) +
    
    # Then, plot significant lipids on top, colored by lipid type
    geom_point(data=df_sig, 
               aes(color=type),
               size=3, alpha=0.5) +
    
    # Add labels for the first 10 significant lipids
    ggrepel::geom_text_repel(data=head(df_sig, 20),
                             aes(label=lipid),  # Assuming 'lipid' is the column with lipid names
                             size=3, 
                             nudge_y = 0.2,    # Adjust as needed
                             nudge_x = 0.2) +
    
    # Color scale for lipid types using lipid_class_colors
    scale_color_manual(values=lipid_class_colors) +
    
    # Themes and title
    theme_bw() +
    labs(
      x = paste("(", Title2, ")", sep=""),
      y = paste("(", Title1, ")", sep=""),
      title = paste(Title1, " vs ", Title2)
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_all
  
  
  ggsave(filename = paste0("scatter_plots_no_log_10/", Title1, "_vs_", Title2, "_scatter_plot_by_type_blank_subtraction.pdf"), plot = plot_all)
  
  ggsave(filename = paste0("scatter_plots_no_log_10/", Title1, "_vs_", Title2, "_scatter_plot_by_type_blank_subtraction.png"), plot = plot_all)
}





# 
# plot_scatter <- function(df, Title1, Title2) {
#   lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
#   lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
#   
#   
#   
#   # Create a named vector to map lipid classes to their colors
#   lipid_class_colors <- setNames(lipid_colors, lipid_classes)
#   
#   
#   # Extract column lengths
#   len1 <- as.numeric(df$Length1[1])
#   len2 <- as.numeric(df$Length2[1])
#   
#   # Columns after 'type' are the ones of interest
#   all_cols <- colnames(df)
#   # Exclude the last column
#   all_cols <- colnames(df)[-length(df)]
#   start_idx <- which(all_cols == "type") + 1
#   cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   
#   
#   # # Get the blank column (which is the penultimate column)
#   # blank_col <- all_cols[length(all_cols)]
#   # 
#   # # Subtract the blank column from the columns of interest and set negative values to 0
#   # df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   # df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   # 
#   # df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   # df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   # 
#   
#   # Compute the means for the groups of columns
#   df <- df %>%
#     mutate(mean1 = rowMeans(select(., one_of(cols1))),
#            mean2 = rowMeans(select(., one_of(cols2))))
#   df_sig <- df %>%
#     filter(FDR < 0.1) 
#   
#   
#   plot_all <- df %>%
#     ggplot(aes(x=log10(mean2), y=log10(mean1))) +
#     # First, plot all points in grey
#     geom_point(aes(color=NULL), color="grey", size=3) +
#     
#     # Then, plot significant lipids on top, colored by logFC
#     geom_point(data=df_sig, 
#                aes(color=logFC),
#                size=3) +
#     
#     # Add labels for the first 10 significant lipids
#     ggrepel::geom_text_repel(data=head(df_sig, 20),
#                              aes(label=lipid),  # Assuming 'lipid' is the column with lipid names
#                              size=3, 
#                              nudge_y = 0.2,    # Adjust as needed
#                              nudge_x = 0.2) +
#     
#     # Color scale for logFC with squishing out-of-bounds values
#     scale_color_gradient2(midpoint=0, low="Blue", mid="#DC9313", 
#                           high="Red", space="Lab", 
#                           limits=c(-1, 1), breaks=c(-1, -0.5, 0, 0.5, 1), 
#                           oob=squish) +
#     
#     # Themes and title
#     theme_bw() +
#     labs(
#       x = paste("log10(", Title2, ")", sep=""),
#       y = paste("log10(", Title1, ")", sep=""),
#       title = paste(Title1, " vs ", Title2)
#     ) +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   plot_all
#   
#   
#   ggsave(filename = paste0("scatter_plots/", Title1, "_vs_", Title2, "_scatter_plot.pdf"), plot = plot_all)
#   
#   ggsave(filename = paste0("scatter_plots/", Title1, "_vs_", Title2, "_scatter_plot.svg"), plot = plot_all)
# }

# 
# plot_scatter_by_type <- function(df, Title1, Title2) {
#   lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
#   lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
#   
#   # Create a named vector to map lipid classes to their colors
#   lipid_class_colors <- setNames(lipid_colors, lipid_classes)
#   
#   # Extract column lengths
#   len1 <- as.numeric(df$Length1[1])
#   len2 <- as.numeric(df$Length2[1])
#   
#   # Columns after 'type' are the ones of interest
#   all_cols <- colnames(df)
#   # Exclude the last column
#   all_cols <- colnames(df)[-length(df)]
#   start_idx <- which(all_cols == "type") + 1
#   cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   
#   
#   # # Get the blank column (which is the penultimate column)
#   # blank_col <- all_cols[length(all_cols)]
#   # 
#   # # Subtract the blank column from the columns of interest and set negative values to 0
#   # df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   # df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   # 
#   # df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   # df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   # 
#   # 
#   # Compute the means for the groups of columns
#   df <- df %>%
#     mutate(mean1 = rowMeans(select(., one_of(cols1))),
#            mean2 = rowMeans(select(., one_of(cols2))))
#   df_sig <- df %>%
#     filter(FDR < 0.1) 
#   
#   # Calculate the axis limits
#   # x_range <- range(log10(df$mean2), na.rm=TRUE)
#   # y_range <- range(log10(df$mean1), na.rm=TRUE)
#   # 
#   # extended_x_limits <- c(x_range[1] - 0.5, x_range[2] + 0.5)
#   # extended_y_limits <- c(y_range[1] - 0.5, y_range[2] + 0.5)
#   
#   plot_all <- df %>%
#     ggplot(aes(x=log10(mean2), y=log10(mean1))) +
#     # First, plot all points in grey
#     # xlim(extended_x_limits) +
#     # ylim(extended_y_limits)+
#     geom_point(aes(color=NULL), color="grey", size=3) +
#     
#     # Then, plot significant lipids on top, colored by lipid type
#     geom_point(data=df_sig, 
#                aes(color=type),
#                size=3, alpha=0.5) +
#     
#     # Add labels for the first 10 significant lipids
#     ggrepel::geom_text_repel(data=head(df_sig, 20),
#                              aes(label=lipid),  # Assuming 'lipid' is the column with lipid names
#                              size=3, 
#                              nudge_y = 0.2,    # Adjust as needed
#                              nudge_x = 0.2) +
#     
#     # Color scale for lipid types using lipid_class_colors
#     scale_color_manual(values=lipid_class_colors) +
#     
#     # Themes and title
#     theme_bw() +
#     labs(
#       x = paste("log10(", Title2, ")", sep=""),
#       y = paste("log10(", Title1, ")", sep=""),
#       title = paste(Title1, " vs ", Title2)
#     ) +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   plot_all
#   
#   
#   ggsave(filename = paste0("scatter_plots/", Title1, "_vs_", Title2, "_scatter_plot_by_type.pdf"), plot = plot_all)
#   
#   ggsave(filename = paste0("scatter_plots/", Title1, "_vs_", Title2, "_scatter_plot_by_type.png"), plot = plot_all)
# }
# 

df


jj
df <- read_csv(paste0("results/",jj,sep=""))
jj <-file_list[1]
file_list[2]
for (jj in file_list){
  excel_file <- read_csv(paste0("results/",jj,sep=""))
  dir.create("plots", F)
  jj
  if (!grepl("full", jj, ignore.case = TRUE)) {
    next
  }
  if (!dir.exists("scatter_plots")) {
    dir.create("scatter_plots")
  }
  # title_for_plot <- gsub(":", "_",excel_file$Title[1])
  Title1 <- gsub("[:| ]", "_", excel_file$Title_1[1])
  Title2 <- gsub("[:| ]", "_", excel_file$Title_2[1])
  Title1 <- gsub("__", "_", Title1)
  Title2 <- gsub("__", "_", Title2) 
  # Title1 <-
  
  
  title_for_plot <- paste0(Title1,Title2,sep="_")
  # plot_scatter(excel_file, Title1, Title2)
  # plot_scatter_by_type(excel_file, Title1, Title2)

  plot_scatter_blank_subtraction(excel_file, Title1, Title2)
  # plot_scatter_by_type_blank_subtraction_no_Log(excel_file, Title1, Title2)
  
  
  

}
  

  
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
  #   Blank <- log10(tp[[blank_name]])
  #   
  #   tp %>%
  #     mutate(lipid = make.unique(lipid)) %>%
  #     filter(lipid %in% DElist) %>%
  #     select( -type) %>%
  #     mutate_if(is.numeric, log10) %>%
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




