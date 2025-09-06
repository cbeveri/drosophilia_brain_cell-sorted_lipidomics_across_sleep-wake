# Load necessary libraries
# library(dplyr)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)


library(dplyr)
library(RColorBrewer)
# File path
# Read the csv file


df <- read.csv("5xFAD_within_lipid_data.csv")
min_FDR <- 0.1
df$logFC_FAD_cer_v_hipo <- ifelse(df$FDR_FAD_cer_v_hipo > min_FDR, 0, df$logFC_FAD_cer_v_hipo)
df$logFC_FAD_cor_v_di <- ifelse(df$FDR_FAD_cor_v_di > min_FDR, 0, df$logFC_FAD_cor_v_di)
df$logFC_FAD_cor_v_hippo <- ifelse(df$FDR_FAD_cor_v_hippo > min_FDR, 0, df$logFC_FAD_cor_v_hippo)
df$logFC_FAD_die_vs_hippo <- ifelse(df$FDR_FAD_die_vs_hippo > min_FDR, 0, df$logFC_FAD_die_vs_hippo)
df$logFC_FAD_cer_v_cor <- ifelse(df$FDR_FAD_cer_v_cor > min_FDR, 0, df$logFC_FAD_cer_v_cor)
df$logFC_FAD_cer_v_di <- ifelse(df$FDR_FAD_cer_v_di > min_FDR, 0, df$logFC_FAD_cer_v_di)


###all_5xFAD
names(df)
all_5xFAD <- df[
  (df$FDR_FAD_cer_v_hipo < 0.1) | 
    (df$FDR_FAD_cer_v_cor < 0.1) | (df$FDR_FAD_cor_v_di < 0.1) | (df$FDR_FAD_cor_v_hippo < 0.1) | (df$FDR_FAD_die_vs_hippo < 0.1) | 
    (df$FDR_FAD_cer_v_di < 0.1), 
]

Title1 <-"hippocampus"
Title2 <-"diencephalon"
Title3 <-"cortex"


title <- "5xFAD_all_5xfad"
intensity_columns <- c("logFC_FAD_cer_v_hipo", "logFC_FAD_cer_v_di","logFC_FAD_cer_v_cor","logFC_FAD_cor_v_di","logFC_FAD_cor_v_hippo","logFC_FAD_die_vs_hippo")
# col_titles <- c(Title1, Title2,Title3)
col_titles<- intensity_columns
create_combined_heatmap_logFC_all(all_5xFAD, col_titles, intensity_columns, title, 1, -1)




df <- read.csv("5xFAD_within_lipid_data.csv")
min_FDR <- 0.1
df$logFC_FAD_cer_v_hipo <- ifelse(df$FDR_FAD_cer_v_hipo > min_FDR, 0, df$logFC_FAD_cer_v_hipo)
df$logFC_FAD_cor_v_di <- ifelse(df$FDR_FAD_cor_v_di > min_FDR, 0, df$logFC_FAD_cor_v_di)
df$logFC_FAD_cor_v_hippo <- ifelse(df$FDR_FAD_cor_v_hippo > min_FDR, 0, df$logFC_FAD_cor_v_hippo)
df$logFC_FAD_die_vs_hippo <- ifelse(df$FDR_FAD_die_vs_hippo > min_FDR, 0, df$logFC_FAD_die_vs_hippo)
df$logFC_FAD_cer_v_cor <- ifelse(df$FDR_FAD_cer_v_cor > min_FDR, 0, df$logFC_FAD_cer_v_cor)
df$logFC_FAD_cer_v_di <- ifelse(df$FDR_FAD_cer_v_di > min_FDR, 0, df$logFC_FAD_cer_v_di)


###Cerebellum
names(df)
FAD_cerebellum_vs_all <- df[
  (df$FDR_FAD_cer_v_hipo < 0.1) | 
    (df$FDR_FAD_cer_v_cor < 0.1) | 
    (df$FDR_FAD_cer_v_di < 0.1), 
]

Title1 <-"hippocampus"
Title2 <-"diencephalon"
Title3 <-"cortex"


title <- "5xFAD_cerebellum_vs_other"
intensity_columns <- c("logFC_FAD_cer_v_hipo", "logFC_FAD_cer_v_di","logFC_FAD_cer_v_cor")
col_titles <- c(Title1, Title2,Title3)
create_combined_heatmap_logFC_all(FAD_cerebellum_vs_all, col_titles, intensity_columns, title, 1, -1)



###Cortex
df <- read.csv("5xFAD_within_lipid_data.csv")
min_FDR <- 0.1
df$logFC_FAD_cer_v_hipo <- ifelse(df$FDR_FAD_cer_v_hipo > min_FDR, 0, df$logFC_FAD_cer_v_hipo)
df$logFC_FAD_cor_v_di <- ifelse(df$FDR_FAD_cor_v_di > min_FDR, 0, df$logFC_FAD_cor_v_di)
df$logFC_FAD_cor_v_hippo <- ifelse(df$FDR_FAD_cor_v_hippo > min_FDR, 0, df$logFC_FAD_cor_v_hippo)
df$logFC_FAD_die_vs_hippo <- ifelse(df$FDR_FAD_die_vs_hippo > min_FDR, 0, df$logFC_FAD_die_vs_hippo)
df$logFC_FAD_cer_v_cor <- ifelse(df$FDR_FAD_cer_v_cor > min_FDR, 0, df$logFC_FAD_cer_v_cor)
df$logFC_FAD_cer_v_di <- ifelse(df$FDR_FAD_cer_v_di > min_FDR, 0, df$logFC_FAD_cer_v_di)

names(df)

FAD_cortex_vs_all <- df[
  (df$FDR_FAD_cor_v_di < 0.1) | 
    (df$FDR_FAD_cor_v_hippo < 0.1) | 
    (df$FDR_FAD_cer_v_cor < 0.1), 
]

FAD_cortex_vs_all$logFC_FAD_cer_v_cor <- -1 * FAD_cortex_vs_all$logFC_FAD_cer_v_cor

Title1 <-"hippocampus"
Title2 <-"diencephalon"
Title3 <-"cerebellum"


title <- "5xFAD_cortex_vs_other"
intensity_columns <- c("logFC_FAD_cor_v_hippo", "logFC_FAD_cor_v_di","logFC_FAD_cer_v_cor")
col_titles <- c(Title1, Title2,Title3)
create_combined_heatmap_logFC_all(FAD_cortex_vs_all, col_titles, intensity_columns, title, 1, -1)



###die
df <- read.csv("5xFAD_within_lipid_data.csv")
min_FDR <- 0.1
df$logFC_FAD_cer_v_hipo <- ifelse(df$FDR_FAD_cer_v_hipo > min_FDR, 0, df$logFC_FAD_cer_v_hipo)
df$logFC_FAD_cor_v_di <- ifelse(df$FDR_FAD_cor_v_di > min_FDR, 0, df$logFC_FAD_cor_v_di)
df$logFC_FAD_cor_v_hippo <- ifelse(df$FDR_FAD_cor_v_hippo > min_FDR, 0, df$logFC_FAD_cor_v_hippo)
df$logFC_FAD_die_vs_hippo <- ifelse(df$FDR_FAD_die_vs_hippo > min_FDR, 0, df$logFC_FAD_die_vs_hippo)
df$logFC_FAD_cer_v_cor <- ifelse(df$FDR_FAD_cer_v_cor > min_FDR, 0, df$logFC_FAD_cer_v_cor)
df$logFC_FAD_cer_v_di <- ifelse(df$FDR_FAD_cer_v_di > min_FDR, 0, df$logFC_FAD_cer_v_di)

names(df)

FAD_die_vs_all <- df[
  (df$FDR_FAD_die_vs_hippo < 0.1) | 
    (df$FDR_FAD_cor_v_di < 0.1) | 
    (df$FDR_FAD_cer_v_di < 0.1), 
]

FAD_die_vs_all$logFC_FAD_cor_v_di <- -1 * FAD_die_vs_all$logFC_FAD_cor_v_di
FAD_die_vs_all$logFC_FAD_cer_v_di <- -1 * FAD_die_vs_all$logFC_FAD_cer_v_di

Title1 <-"hippocampus"
Title2 <-"cortex"
Title3 <-"cerebellum"


title <- "5xFAD_dien_vs_other"
intensity_columns <- c("logFC_FAD_die_vs_hippo", "logFC_FAD_cor_v_di","logFC_FAD_cer_v_di")
col_titles <- c(Title1, Title2,Title3)
create_combined_heatmap_logFC_all(FAD_die_vs_all, col_titles, intensity_columns, title, 1, -1)





###hippo
df <- read.csv("5xFAD_within_lipid_data.csv")
min_FDR <- 0.1
df$logFC_FAD_cer_v_hipo <- ifelse(df$FDR_FAD_cer_v_hipo > min_FDR, 0, df$logFC_FAD_cer_v_hipo)
df$logFC_FAD_cor_v_di <- ifelse(df$FDR_FAD_cor_v_di > min_FDR, 0, df$logFC_FAD_cor_v_di)
df$logFC_FAD_cor_v_hippo <- ifelse(df$FDR_FAD_cor_v_hippo > min_FDR, 0, df$logFC_FAD_cor_v_hippo)
df$logFC_FAD_die_vs_hippo <- ifelse(df$FDR_FAD_die_vs_hippo > min_FDR, 0, df$logFC_FAD_die_vs_hippo)
df$logFC_FAD_cer_v_cor <- ifelse(df$FDR_FAD_cer_v_cor > min_FDR, 0, df$logFC_FAD_cer_v_cor)
df$logFC_FAD_cer_v_di <- ifelse(df$FDR_FAD_cer_v_di > min_FDR, 0, df$logFC_FAD_cer_v_di)

names(df)

FAD_hippo_vs_all <- df[
  (df$FDR_FAD_die_vs_hippo < 0.1) | 
    (df$FDR_FAD_cor_v_hippo < 0.1) | 
    (df$FDR_FAD_cer_v_hipo < 0.1), 
]

FAD_hippo_vs_all$logFC_FAD_die_vs_hippo <- -1 * FAD_hippo_vs_all$logFC_FAD_die_vs_hippo
FAD_hippo_vs_all$logFC_FAD_cor_v_hippo <- -1 * FAD_hippo_vs_all$logFC_FAD_cor_v_hippo
FAD_hippo_vs_all$logFC_FAD_cer_v_hipo <- -1 * FAD_hippo_vs_all$logFC_FAD_cer_v_hipo

Title1 <-"diencephalon"
Title2 <-"cortex"
Title3 <-"cerebellum"


title <- "5xFAD_hippo_vs_other"
intensity_columns <- c("logFC_FAD_die_vs_hippo", "logFC_FAD_cor_v_hippo","logFC_FAD_cer_v_hipo")
col_titles <- c(Title1, Title2,Title3)
create_combined_heatmap_logFC_all(FAD_hippo_vs_all, col_titles, intensity_columns, title, 1, -1)






create_combined_heatmap_logFC_all <- function(df, col_titles, intensity_columns, title, limit, limit_lower, min_FDR =0.1) {
  
  if (length(col_titles) <3 || length(intensity_columns) <3) {
    stop("The length of col_titles and intensity_columns must be 3.")
  }
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width=10, height=15) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  tryCatch({
    
    # Update LogFC values based on PValue conditions
    
    # Update LogFC values based on PValue conditions
    # Assuming df is your dataframe and min_FDR is already defined
    # df$logFC_FAD_cer_v_hipo <- ifelse(df$FDR_FAD_cer_v_hipo > min_FDR, 0, df$logFC_FAD_cer_v_hipo)
    # df$logFC_FAD_cor_v_di <- ifelse(df$FDR_FAD_cor_v_di > min_FDR, 0, df$logFC_FAD_cor_v_di)
    # df$logFC_FAD_cor_v_hippo <- ifelse(df$FDR_FAD_cor_v_hippo > min_FDR, 0, df$logFC_FAD_cor_v_hippo)
    # df$logFC_FAD_die_vs_hippo <- ifelse(df$FDR_FAD_die_vs_hippo > min_FDR, 0, df$logFC_FAD_die_vs_hippo)
    # df$logFC_FAD_cer_v_cor <- ifelse(df$FDR_FAD_cer_v_cor > min_FDR, 0, df$logFC_FAD_cer_v_cor)
    # df$logFC_FAD_cer_v_di <- ifelse(df$FDR_FAD_cer_v_di > min_FDR, 0, df$logFC_FAD_cer_v_di)
    
    
    
    # Setup conditions for PValue filtering
    conditions <- list(
      list(filter = TRUE, suffix = "_"),
      list(filter = FALSE, suffix = "")
    )
    
    # Generate the whole heatmap
    for (condition in conditions) {
      if (condition$filter) {
        df_filtered <- df
      } else {
        df_filtered <- df
      }
      

      df_filtered <- df_filtered %>%
        mutate(type = factor(type, levels = c("SM","Cer","FA","CAR","PC","PE","PI","PS","CE", "DAG", "TAG"))) %>%

        arrange(type, df_filtered[[intensity_columns[1]]], df_filtered[[intensity_columns[2]]], df_filtered[[intensity_columns[3]]])

      title_suffix <- condition$suffix
      heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
      # 
      annotation_col_df <- data.frame(Labels = col_titles)
      rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])

      
      first_occurrence <- !duplicated(df_filtered$type)
      
      df_filtered$Row_Label <- ifelse(first_occurrence, as.character(df_filtered$type), "")
      
      
      # annotation_col_df <- data.frame(Labels = annotation_labels)
      # rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
      
      heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                        main = paste0(title, title_suffix),
                                        cluster_rows = FALSE, 
                                        cluster_cols = FALSE, 
                                        show_colnames = TRUE, 
                                        show_rownames = TRUE,
                                        breaks = heatmap_breaks,
                                        border_color = "black", 
                                        fontsize = 10,
                                        fontsize_legend = 25,
                                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                        labels_row = df_filtered$Row_Label,
                                        annotation_col = annotation_col_df)
      # Save the heatmap
      
      save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_new/", title, title_suffix,"Limit_",limit, "_LogFC_All_Comparisons.pdf"), width = 10, height = 15)
    }
    
    # Generate heatmaps for each unique 'type'
    unique_types <- unique(df$type)
    for (single_type in unique_types) {
      
      df_type <- df %>% filter(type == single_type)
      
      for (condition in conditions) {
        if (condition$filter) {
          df_filtered <- df_type 
        } else {
          df_filtered <- df_type
        }
        
        df_filtered <- df_filtered %>% arrange(type, df_filtered[[intensity_columns[1]]], df_filtered[[intensity_columns[2]]], df_filtered[[intensity_columns[3]]])
        
        title_suffix <- condition$suffix
        heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
        
        annotation_col_df <- data.frame(Labels = col_titles)
        rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
        
        heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                          main = paste0(title, title_suffix, " - ", single_type),
                                          cluster_rows = FALSE, 
                                          cluster_cols = FALSE, 
                                          show_colnames = TRUE, 
                                          show_rownames = TRUE,
                                          breaks = heatmap_breaks,
                                          border_color = "black", 
                                          fontsize = 10,
                                          fontsize_legend = 25,
                                          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                          labels_row = df_filtered$lipid,
                                          annotation_col = annotation_col_df)
        
        # Save the heatmap
        save_pheatmap_pdf(heatmap_obj, paste0("by_class_heatmap/", title, title_suffix, "_", single_type,"Limit_",limit, "_LogFC_All_Comparisons.pdf"), width = 10, height = 25)
      }
    }
    
  }, error = function(e) {
    print(paste("Error processing file: Error message:", e$message))
  })
}




