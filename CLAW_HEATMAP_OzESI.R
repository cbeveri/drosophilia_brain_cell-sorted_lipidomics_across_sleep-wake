library(tidyverse)
library(pheatmap)
library(RColorBrewer)


library(dplyr)
library(RColorBrewer)
# File path


create_heatmap_Z_score <- function(file_path) {
  
  # Read and preprocess the data
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Arrange dataframe
  df <- df %>%
    arrange(lipid)
  
  # Identify intensity columns
  intensity_columns <- setdiff(names(df), c("lipid", "logFC", "logCPM", "LR", "PValue", "FDR", "Length1", "Length2", "Title_1", "Title_2", "type", tail(names(df), 1)))
  
  # Get the blank column (assuming it's the last column in the dataframe)
  blank_col <- tail(names(df), 1)
  
  # Subtract the blank column from the intensity columns and set negative values to 0
  df[intensity_columns] <- pmax(df[intensity_columns] - df[[blank_col]], 0)
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.infinite(.), NA)))
  
  
  
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.na(.), 0)))
  new_folder <- "csv_blank_subtracted"
  
  # Create the new folder path
  new_folder_path <- file.path(getwd(), new_folder)
  
  # Check if the folder exists, and create it if it doesn't
  if (!dir.exists(new_folder_path)) {
    dir.create(new_folder_path)
  }
  
  # Define the new file path
  new_file_name <- basename(file_path)
  new_file_path <- file.path(new_folder_path, new_file_name)
  
  # Save the dataframe to the new file path
  write.csv(df, new_file_path, row.names = FALSE)
  
  # Output the new file path
  # new_file_path
  
  # Remove rows where all values in intensity_columns are zeros
  df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
  
  ## Inserted code to divide each value in intensity_columns by the sum of that column
  df <- df %>% 
    mutate(across(all_of(intensity_columns), 
                  ~ . / sum(. , na.rm = TRUE)))
  
  # Compute Z-scores
  df <- df %>%
    rowwise() %>%
    mutate(across(all_of(intensity_columns),
                  ~ ( . - mean(c_across(all_of(intensity_columns)), na.rm = TRUE)) /
                    sd(c_across(all_of(intensity_columns)), na.rm = TRUE))) %>%
    ungroup()

  # df <- df %>%
    # filter(type %in% c("FA", "TAG", "DAG"))
  
  title <- sub("_full.csv", "", basename(file_path))
  labels_title_1 <- unique(df$Title_1)
  labels_title_2 <- unique(df$Title_2)
  annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  # Setup conditions for FDR and PValue filtering
  conditions <- list(
    list(filter = "FDR", threshold = 0.1, suffix = "_FDR_01_"),
    list(filter = "PValue", threshold = 0.05, suffix = "_PVALUE_05_"),
    list(filter = "PValue", threshold = 0.01, suffix = "_PVALUE_01_"),
    list(filter = "NONE", threshold = NULL, suffix = "")
  )
  
  # Loop through conditions to generate and save heatmaps
  for (condition in conditions) {
    
    # Filter data based on conditions
    if (condition$filter == "FDR") {
      df_filtered <- df %>% filter(FDR < condition$threshold)
    } else if (condition$filter == "PValue") {
      df_filtered <- df %>% filter(PValue < condition$threshold)
    } else {
      df_filtered <- df
    }
    
    if(nrow(df_filtered) < 1) {
      next
    }
    
    df_filtered <- df_filtered %>% arrange(type,logFC)
    
    # df_filtered <- df_filtered %>% mutate(type = factor(type, levels = c("FA", "DAG", "TAG"))) %>% arrange(type)
    # df_filtered <- df_filtered %>% mutate(type = factor(type, levels = c("SM","Cer","FA","CAR","PC","PE","PI","PS","CE", "DAG", "TAG"))) %>% arrange(type)
    # 
    # title_suffix <- condition$suffix
    # heatmap_breaks <- seq(-2, 2, length.out = 100)
    # annotation_data <- data.frame(type = df_filtered$type)
    # first_occurrence <- !duplicated(df_filtered$type)
    # df_filtered$Row_Label <- ifelse(first_occurrence, df_filtered$type, "")
    df_filtered <- df_filtered %>%
      mutate(type = factor(type, levels = c("SM","Cer","FA","CAR","PC","PE","PI","PS","CE", "DAG", "TAG"))) %>%
      arrange(type, logFC)  # Sorting by type first, then by logFC
    
    title_suffix <- condition$suffix
    heatmap_breaks <- seq(-2, 2, length.out = 100)
    annotation_data <- data.frame(type = df_filtered$type)
    first_occurrence <- !duplicated(df_filtered$type)
    
    df_filtered$Row_Label <- ifelse(first_occurrence, as.character(df_filtered$type), "")
    
    
    annotation_col_df <- data.frame(Labels = annotation_labels)
    rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
    
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
                                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "BrBG")))(100),
                                      labels_row = df_filtered$Row_Label,
                                      annotation_col = annotation_col_df)
    
    # Save the heatmap
    save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large/", title, title_suffix, "Zscore_2_limit_All_only.pdf"), width = 11, height = 10)
  }
}


create_heatmap_by_class <- function(file_path) {
  
  # Read and preprocess the data
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Arrange dataframe by lipid
  df <- df %>%
    arrange(lipid)
  # # Read the names to be replaced from the Excel file
  # names_to_fix <- readxl::read_xlsx("Names_to_fix.xlsx")
  # 
  # names_to_fix <- names_to_fix[!is.na(names_to_fix$New_name) & names_to_fix$New_name != "",]
  # name_changes <- setNames(names_to_fix$New_name, names_to_fix$lipid)
  # 
  # # Replace lipid names in df using the named vector
  # df$lipid <- ifelse(df$lipid %in% names(name_changes), name_changes[df$lipid], df$lipid)
  # 
  
  # Arrange dataframe by lipid
  df <- df %>%
    arrange(lipid)
  
  # Identify intensity columns
  intensity_columns <- setdiff(names(df), c("lipid", "logFC", "logCPM", "LR", "PValue", "FDR", "Length1", "Length2", "Title_1", "Title_2", "type", tail(names(df), 1)))
  
  # Get the blank column (assuming it's the last column in the dataframe)
  blank_col <- tail(names(df), 1)
  
  # Subtract the blank column from the intensity columns and set negative values to 0
  df[intensity_columns] <- pmax(df[intensity_columns] - df[[blank_col]], 0)
  # Replace infinite values with NA for data frames
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.infinite(.), NA)))
  
  # Replace NAs with 0 in intensity columns
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.na(.), 0)))
  
  # Remove rows where all values in intensity_columns are zeros
  df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
  
  ## Inserted code to divide each value in intensity_columns by the sum of that column
  df <- df %>% 
    mutate(across(all_of(intensity_columns), 
                  ~ . / sum(. , na.rm = TRUE)))
  
  # Compute Z-scores
  df <- df %>%
    rowwise() %>%
    mutate(across(all_of(intensity_columns),
                  ~ ( . - mean(c_across(all_of(intensity_columns)), na.rm = TRUE)) /
                    sd(c_across(all_of(intensity_columns)), na.rm = TRUE))) %>%
    ungroup()
  
  
  title <- sub("_full.csv", "", basename(file_path))
  labels_title_1 <- unique(df$Title_1)
  labels_title_2 <- unique(df$Title_2)
  annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  # Setup conditions for FDR filtering
  conditions <- list(
    list(filter = TRUE, suffix = "_FDR_by_class"),
    list(filter = FALSE, suffix = "_by_class")
  )
  
  # Loop through each unique 'type'
  for (lipid_type in unique(df$type)) {
    df_filtered_by_type <- df %>% filter(type == lipid_type)
    
    # Loop through conditions to generate and save heatmaps
    for (condition in conditions) {
      if (condition$filter) {
        df_filtered <- df_filtered_by_type %>% filter(FDR < 0.1)
      } else {
        df_filtered <- df_filtered_by_type
      }
      if(nrow(df_filtered) < 1) {
        next
      }
      print(colnames(df_filtered))
      # Arrange dataframe by lipid
      df_filtered <- df_filtered %>%
        arrange(logFC)
      
      title_suffix <- condition$suffix
      heatmap_breaks <- seq(-2, 2, length.out = 100)
      annotation_data <- data.frame(type = df_filtered$type)
      
      
      annotation_col_df <- data.frame(Labels = annotation_labels)
      rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
      # Using lipid as row label
      # rownames(df_filtered) <- df_filtered$lipid
      # print(rownames(df_filtered))
      heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                        main = paste0(title, "_", lipid_type, title_suffix),
                                        cluster_rows = FALSE, 
                                        cluster_cols = FALSE, 
                                        show_colnames = TRUE, 
                                        show_rownames = TRUE,
                                        breaks = heatmap_breaks,
                                        border_color = "black", 
                                        fontsize = 10,
                                        fontsize_legend = 25,
                                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "BrBG")))(100),
                                        labels_row = df_filtered$lipid, # Use lipid for row labels
                                        annotation_col = annotation_col_df)
      
      # Save the heatmap
      save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_by_class/", title, "_", lipid_type, title_suffix, ".pdf"), width = 11, height = 10)
    }
  }
}



# Check if heatmaps folder exists, if not, create it
if (!dir.exists("heatmaps_large")) {
  dir.create("heatmaps_large")
}
if (!dir.exists("heatmaps_large_by_class")) {
  dir.create("heatmaps_large_by_class")
}

# Get the list of files with "_full.csv"
files <- list.files(path = "results", pattern = "_full.csv$", full.names = TRUE)
getwd()
# Create heatmap for each file
# lapply(files, create_heatmap_with_blank)

# Create heatmap for each file
lapply(files, create_heatmap_Z_score)
lapply(files, create_heatmap_by_class)









create_heatmap_Z_score_certain_lipids <- function(file_path) {
  
  # Read and preprocess the data
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Arrange dataframe
  df <- df %>%
    arrange(lipid)
  
  # Identify intensity columns
  intensity_columns <- setdiff(names(df), c("lipid", "logFC", "logCPM", "LR", "PValue", "FDR", "Length1", "Length2", "Title_1", "Title_2", "type", tail(names(df), 1)))
  
  # Get the blank column (assuming it's the last column in the dataframe)
  blank_col <- tail(names(df), 1)
  
  # Subtract the blank column from the intensity columns and set negative values to 0
  df[intensity_columns] <- pmax(df[intensity_columns] - df[[blank_col]], 0)
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.infinite(.), NA)))
  
  
  
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.na(.), 0)))
  new_folder <- "csv_blank_subtracted"
  
  # Create the new folder path
  new_folder_path <- file.path(getwd(), new_folder)
  
  # Check if the folder exists, and create it if it doesn't
  if (!dir.exists(new_folder_path)) {
    dir.create(new_folder_path)
  }
  
  # Define the new file path
  new_file_name <- basename(file_path)
  new_file_path <- file.path(new_folder_path, new_file_name)
  
  # Save the dataframe to the new file path
  write.csv(df, new_file_path, row.names = FALSE)
  
  # Output the new file path
  # new_file_path
  
  # Remove rows where all values in intensity_columns are zeros
  df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
  
  ## Inserted code to divide each value in intensity_columns by the sum of that column
  df <- df %>% 
    mutate(across(all_of(intensity_columns), 
                  ~ . / sum(. , na.rm = TRUE)))
  
  # Compute Z-scores
  df <- df %>%
    rowwise() %>%
    mutate(across(all_of(intensity_columns),
                  ~ ( . - mean(c_across(all_of(intensity_columns)), na.rm = TRUE)) /
                    sd(c_across(all_of(intensity_columns)), na.rm = TRUE))) %>%
    ungroup()
  
  df <- df %>%
  filter(type %in% c("FA", "PS", "PI","PE"))
  
  title <- sub("_full.csv", "", basename(file_path))
  labels_title_1 <- unique(df$Title_1)
  labels_title_2 <- unique(df$Title_2)
  annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  # Setup conditions for FDR and PValue filtering
  conditions <- list(
    list(filter = "FDR", threshold = 0.1, suffix = "_FDR_01_"),
    list(filter = "PValue", threshold = 0.05, suffix = "_PVALUE_05_"),
    list(filter = "PValue", threshold = 0.01, suffix = "_PVALUE_01_"),
    list(filter = "NONE", threshold = NULL, suffix = "")
  )
  
  # Loop through conditions to generate and save heatmaps
  for (condition in conditions) {
    
    # Filter data based on conditions
    if (condition$filter == "FDR") {
      df_filtered <- df %>% filter(FDR < condition$threshold)
    } else if (condition$filter == "PValue") {
      df_filtered <- df %>% filter(PValue < condition$threshold)
    } else {
      df_filtered <- df
    }
    
    if(nrow(df_filtered) < 1) {
      next
    }
    
    df_filtered <- df_filtered %>% arrange(type,logFC)
    
    # df_filtered <- df_filtered %>% mutate(type = factor(type, levels = c("FA", "DAG", "TAG"))) %>% arrange(type)
    # df_filtered <- df_filtered %>% mutate(type = factor(type, levels = c("SM","Cer","FA","CAR","PC","PE","PI","PS","CE", "DAG", "TAG"))) %>% arrange(type)
    # 
    # title_suffix <- condition$suffix
    # heatmap_breaks <- seq(-2, 2, length.out = 100)
    # annotation_data <- data.frame(type = df_filtered$type)
    # first_occurrence <- !duplicated(df_filtered$type)
    # df_filtered$Row_Label <- ifelse(first_occurrence, df_filtered$type, "")
    df_filtered <- df_filtered %>%
      mutate(type = factor(type, levels = c("FA","PS","PI","PE"))) %>%
      arrange(type, logFC)  # Sorting by type first, then by logFC
    
    title_suffix <- condition$suffix
    heatmap_breaks <- seq(-2, 2, length.out = 100)
    annotation_data <- data.frame(type = df_filtered$type)
    first_occurrence <- !duplicated(df_filtered$type)
    
    df_filtered$Row_Label <- ifelse(first_occurrence, as.character(df_filtered$type), "")
    
    
    annotation_col_df <- data.frame(Labels = annotation_labels)
    rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
    
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
                                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "BrBG")))(100),
                                      labels_row = df_filtered$Row_Label,
                                      annotation_col = annotation_col_df)
    
    # Save the heatmap
    save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_all_togehter/", title, title_suffix, "Zscore_2_limit_All_only_FA_PS_PI_PE.pdf"), width = 11, height = 10)
  }
}




create_heatmap_Z_score_TG <- function(file_path) {
  
  # Read and preprocess the data
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Arrange dataframe
  df <- df %>%
    arrange(lipid)
  
  # Identify intensity columns
  intensity_columns <- setdiff(names(df), c("lipid", "logFC", "logCPM", "LR", "PValue", "FDR", "Length1", "Length2", "Title_1", "Title_2", "type", tail(names(df), 1)))
  
  # Get the blank column (assuming it's the last column in the dataframe)
  blank_col <- tail(names(df), 1)
  
  # Subtract the blank column from the intensity columns and set negative values to 0
  df[intensity_columns] <- pmax(df[intensity_columns] - df[[blank_col]], 0)
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.infinite(.), NA)))
  
  
  
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.na(.), 0)))
  
  
  # Remove rows where all values in intensity_columns are zeros
  df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
  
  ## Inserted code to divide each value in intensity_columns by the sum of that column
  df <- df %>% 
    mutate(across(all_of(intensity_columns), 
                  ~ . / sum(. , na.rm = TRUE)))
  
  # Compute Z-scores
  df <- df %>%
    rowwise() %>%
    mutate(across(all_of(intensity_columns),
                  ~ ( . - mean(c_across(all_of(intensity_columns)), na.rm = TRUE)) /
                    sd(c_across(all_of(intensity_columns)), na.rm = TRUE))) %>%
    ungroup()
  
  df <- df %>% filter(type == "TAG")# %>% filter(FDR < 0.1)
  
  # Extract the first two numbers (two digits) after [TG(
  df$TG_length <- as.numeric(str_extract(df$lipid, "(?<=\\[TG\\()\\d{2}"))
  
  
  
  
  
  # df <- df %>%
  # filter(type %in% c("FA", "TAG", "DAG"))
  
  title <- sub("_full.csv", "", basename(file_path))
  labels_title_1 <- unique(df$Title_1)
  labels_title_2 <- unique(df$Title_2)
  annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  # Setup conditions for FDR and PValue filtering
  conditions <- list(
    list(filter = "FDR", threshold = 0.1, suffix = "_FDR_01_"),
    list(filter = "PValue", threshold = 0.05, suffix = "_PVALUE_05_"),
    list(filter = "PValue", threshold = 0.01, suffix = "_PVALUE_01_"),
    list(filter = "NONE", threshold = NULL, suffix = "")
  )
  
  # Loop through conditions to generate and save heatmaps
  for (condition in conditions) {
    
    # Filter data based on conditions
    if (condition$filter == "FDR") {
      df_filtered <- df %>% filter(FDR < condition$threshold)
    } else if (condition$filter == "PValue") {
      df_filtered <- df %>% filter(PValue < condition$threshold)
    } else {
      df_filtered <- df
    }
    
    if(nrow(df_filtered) < 1) {
      next
    }
    
    df_filtered <- df_filtered %>% arrange(TG_length)
    
    # df_filtered <- df_filtered %>% mutate(type = factor(type, levels = c("FA", "DAG", "TAG"))) %>% arrange(type)
    
    
    title_suffix <- condition$suffix
    heatmap_breaks <- seq(-2, 2, length.out = 100)
    annotation_data <- data.frame(TG_length = df_filtered$TG_length)
    first_occurrence <- !duplicated(df_filtered$TG_length)
    df_filtered$Row_Label <- ifelse(first_occurrence, df_filtered$TG_length, "")
    
    annotation_col_df <- data.frame(Labels = annotation_labels)
    rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
    
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
                                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "BrBG")))(100),
                                      labels_row = df_filtered$Row_Label,
                                      annotation_col = annotation_col_df)
    
    # Save the heatmap
    save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large/", title, title_suffix, "Zscore_2_limit_All_only.pdf"), width = 11, height = 10)
  }
}


create_heatmap_Z_score_TG_FA_length <- function(file_path) {
  
  # Read and preprocess the data
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Arrange dataframe
  df <- df %>%
    arrange(lipid)
  
  # Identify intensity columns
  intensity_columns <- setdiff(names(df), c("lipid", "logFC", "logCPM", "LR", "PValue", "FDR", "Length1", "Length2", "Title_1", "Title_2", "type", tail(names(df), 1)))
  
  # Get the blank column (assuming it's the last column in the dataframe)
  blank_col <- tail(names(df), 1)
  
  # Subtract the blank column from the intensity columns and set negative values to 0
  df[intensity_columns] <- pmax(df[intensity_columns] - df[[blank_col]], 0)
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.infinite(.), NA)))
  
  
  
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.na(.), 0)))
  
  
  # Remove rows where all values in intensity_columns are zeros
  df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
  
  ## Inserted code to divide each value in intensity_columns by the sum of that column
  df <- df %>% 
    mutate(across(all_of(intensity_columns), 
                  ~ . / sum(. , na.rm = TRUE)))
  
  # Compute Z-scores
  df <- df %>%
    rowwise() %>%
    mutate(across(all_of(intensity_columns),
                  ~ ( . - mean(c_across(all_of(intensity_columns)), na.rm = TRUE)) /
                    sd(c_across(all_of(intensity_columns)), na.rm = TRUE))) %>%
    ungroup()
  
  df <- df %>% filter(type == "TAG")# %>% filter(FDR < 0.1)
  
  # Extract the first two numbers (two digits) after [TG(
  df$TG_length <- as.numeric(str_extract(df$lipid, "(?<=\\[TG\\()\\d{2}"))
  
  df <- df %>%
    # mutate(mean1 = rowMeans(select(., all_of(cols1))),
    #        mean2 = rowMeans(select(., all_of(cols2)))) %>%
    # filter(type == "TAG") %>%filter(FDR < 0.1) %>%
    mutate(lipid_sub = str_extract(lipid, "(?<=_FA).+$")) %>%
    mutate(lipid_sub = ifelse(lipid_sub %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", lipid_sub))
  
  
  
  # df <- df %>%
  # filter(type %in% c("FA", "TAG", "DAG"))
  
  title <- sub("_full.csv", "", basename(file_path))
  labels_title_1 <- unique(df$Title_1)
  labels_title_2 <- unique(df$Title_2)
  annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  # Setup conditions for FDR and PValue filtering
  conditions <- list(
    list(filter = "FDR", threshold = 0.1, suffix = "_FDR_01_"),
    list(filter = "PValue", threshold = 0.05, suffix = "_PVALUE_05_"),
    list(filter = "PValue", threshold = 0.01, suffix = "_PVALUE_01_"),
    list(filter = "NONE", threshold = NULL, suffix = "")
  )
  
  # Loop through conditions to generate and save heatmaps
  for (condition in conditions) {
    
    # Filter data based on conditions
    if (condition$filter == "FDR") {
      df_filtered <- df %>% filter(FDR < condition$threshold)
    } else if (condition$filter == "PValue") {
      df_filtered <- df %>% filter(PValue < condition$threshold)
    } else {
      df_filtered <- df
    }
    
    if(nrow(df_filtered) < 1) {
      next
    }
    
    df_filtered <- df_filtered %>% arrange(lipid_sub)
    
    # df_filtered <- df_filtered %>% mutate(type = factor(type, levels = c("FA", "DAG", "TAG"))) %>% arrange(type)
    
    
    title_suffix <- condition$suffix
    heatmap_breaks <- seq(-2, 2, length.out = 100)
    annotation_data <- data.frame(lipid_sub = df_filtered$lipid_sub)
    first_occurrence <- !duplicated(df_filtered$lipid_sub)
    df_filtered$Row_Label <- ifelse(first_occurrence, df_filtered$lipid_sub, "")
    
    annotation_col_df <- data.frame(Labels = annotation_labels)
    rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
    
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
                                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "BrBG")))(100),
                                      labels_row = df_filtered$Row_Label,
                                      annotation_col = annotation_col_df)
    
    # Save the heatmap
    save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large/", title, title_suffix, "Zscore_2_limit_All_only_FA_TG_length.pdf"), width = 11, height = 10)
  }
}

