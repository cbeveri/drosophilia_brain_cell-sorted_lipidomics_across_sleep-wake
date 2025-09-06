
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

FA_plot_lengths <- function(df, Title1, Title2) {
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
  
  # Step 1: Filter to only include FA types
  df <- df %>%
    filter(type == "FA")
  
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
  
  # Set NA values to 0
  df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
  df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
  
  # Step 2: Compute the means for the groups of columns
  df <- df %>%
    mutate(mean1 = rowMeans(select(., one_of(cols1))),
           mean2 = rowMeans(select(., one_of(cols2))))
  
  # Step 3: Categorize FA lipids by chain length based on the 'type' column
  df <- df %>%
    mutate(chain_length = sapply(lipid, categorize_chain_length))
  
  # Group by chain length and calculate the sum of means for each group
  df_all <- df %>%
    group_by(chain_length) %>%
    summarise(sum_mean1 = sum(mean1), 
              sum_mean2 = sum(mean2), .groups = "keep")
  
  # Reshape for plotting
  df_all_long <- df_all %>%
    select(chain_length, sum_mean1, sum_mean2) %>%
    gather(key = "group", value = "value", -chain_length)
  
  # Plotting stacked bar for FA chain lengths
  plot_fa <- df_all_long %>%
    group_by(chain_length, group) %>%
    ungroup() %>%
    ggplot(aes(x = group, y = value, fill = chain_length)) +
    geom_col() +
    theme_classic() +
    ylab("Total lipid content ") +
    scale_fill_manual(values = setNames(chain_length_colors, chain_length_classes), name = "Chain Length") +
    scale_x_discrete(labels = c(Title1, Title2)) +
    labs(title = paste(Title1, "vs", Title2, " (FA Chain Length)"), y = "Sum of Means")
  plot_fa
  # Save the plot
  ggsave(filename = paste0("FA_lengths_stacked_bar/", Title1, "_vs_", Title2, "_STACKED_FA_ChainLength.pdf"), plot = plot_fa)
}

library(tidyverse)
library(tidyverse)




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
  FA_plot_lengths(excel_file, Title1, Title2)

  # Plotting and saving
  # plot_object <- plot_significant_lipids(excel_file, title_for_plot)
  # process_and_plot(excel_file, Title1, Title2)
  # Ensuring the plots directory exists

  

}
  
