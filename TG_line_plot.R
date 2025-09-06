








library(ggplot2)
# library(dplyr)
library(stringr)






plot_FA_distribution_individual <- function(df, Title1, Title2,filename) {
  # Filter the dataframe
  df_filtered <- df %>% filter(type == "FA")# %>% filter(FDR < 0.1)
  
  # Extract the first two numbers (two digits) after [TG(
  df_filtered$TG_length <- as.numeric(str_extract(df_filtered$lipid, "(?<=FA\\()\\d{1,2}(?=:)"))
  
  # Aggregate the data
  df_summarized <- df_filtered %>% 
    group_by(TG_length) %>% 
    summarize(total_mean1 = sum(mean1), total_mean2 = sum(mean2))
  
  # Calculate percentages
  df_summarized$percent_mean1 <- df_summarized$total_mean1 / sum(df_summarized$total_mean1)
  df_summarized$percent_mean2 <- df_summarized$total_mean2 / sum(df_summarized$total_mean2)
  
  directory_name <- "FA_distribution"
  file_name <- paste0("FA_", filename, ".csv")
  full_path <- file.path(directory_name, file_name)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
  }
  
  # Save the summarized data frame to a CSV file
  write.csv(df_summarized, full_path, row.names = FALSE)
  
  # Create the line plot
  p<- ggplot(df_summarized, aes(x = TG_length)) +
    geom_line(aes(y = percent_mean1), color = "blue") +
    geom_line(aes(y = percent_mean2), color = "red") +
    labs(title = Title1, subtitle = Title2, x = "FA Length", y = "Percentage") +
    theme_minimal()
  
  ggsave(filename, plot = p, width = 12, height = 8, units = "in")
  
  
}



plot_FA_distribution_with_saturation <- function(df, Title1, Title2, filename) {
  # Filter the dataframe for 'FA' type
  df_filtered <- df %>% filter(type == "FA")
  
  # Extract FA length and saturation level from the 'lipid' column
  df_filtered <- df_filtered %>%
    mutate(
      FA_length = as.numeric(str_extract(lipid, "(?<=FA\\()\\d{1,2}(?=:)")),
      saturation_level = as.numeric(str_extract(lipid, "(?<=:)\\d{1,2}(?=\\))")),
      FA_label = paste0(FA_length, ":", saturation_level)
    )
  
  # Remove rows with missing FA_length or saturation_level
  df_filtered <- df_filtered %>% filter(!is.na(FA_length) & !is.na(saturation_level))
  
  # Check if there are any data left
  if (nrow(df_filtered) == 0) {
    stop("No valid FA data after extraction. Please check the lipid column format.")
  }
  
  # Aggregate the data
  df_summarized <- df_filtered %>%
    group_by(FA_length, saturation_level, FA_label) %>%
    summarize(
      total_mean1 = sum(mean1, na.rm = TRUE),
      total_mean2 = sum(mean2, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Remove FA labels where both total_mean1 and total_mean2 are zero
  # df_summarized <- df_summarized %>% filter(total_mean1 > 0 | total_mean2 > 0)
  
  # Calculate percentages
  total_mean1_sum <- sum(df_summarized$total_mean1)
  total_mean2_sum <- sum(df_summarized$total_mean2)
  df_summarized <- df_summarized %>%
    mutate(
      percent_mean1 = total_mean1 / total_mean1_sum,
      percent_mean2 = total_mean2 / total_mean2_sum
    )
  
  # Order FA labels numerically based on FA_length and saturation_level
  df_summarized <- df_summarized %>%
    arrange(FA_length, saturation_level)
  
  # Set FA_label as a factor with levels in the correct order
  df_summarized$FA_label <- factor(df_summarized$FA_label, levels = unique(df_summarized$FA_label))
  
  # Create new directory for outputs
  directory_name <- "FA_l_s_p"
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
  }
  
  # Prepare filenames
  csv_filename <- paste0("FA_saturation_", filename, ".csv")
  plot_filename <- paste0("FA_saturation_", filename, ".pdf")
  csv_full_path <- file.path(directory_name, csv_filename)
  plot_full_path <- file.path(directory_name, plot_filename)
  
  # Save the summarized data frame to a CSV file
  write.csv(df_summarized, csv_full_path, row.names = FALSE)
  
  # Prepare data for plotting
  df_melted <- df_summarized %>%
    select(FA_label, percent_mean1, percent_mean2) %>%
    pivot_longer(cols = c(percent_mean1, percent_mean2), names_to = "Group", values_to = "Percentage")
  
  # Replace 'percent_mean1' and 'percent_mean2' with actual group names
  df_melted$Group <- recode(df_melted$Group,
                            'percent_mean1' = Title1,
                            'percent_mean2' = Title2)
  
  # Plotting
  p <- ggplot(df_melted, aes(x = FA_label, y = Percentage, group = Group, color = Group)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
      title = paste0("FA Distribution with Saturation - ", Title1, " vs ", Title2),
      x = "FA Length: Saturation Level",
      y = "Percentage",
      color = "Group"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Display the plot
  print(p)
  
  # Save the plot
  ggsave(plot_full_path, plot = p, width = 12, height = 8, units = "in")
}



plot_TG_distribution_individual <- function(df, Title1, Title2,filename) {
  # Filter the dataframe
  df_filtered <- df %>% filter(type == "TAG")# %>% filter(FDR < 0.1)
  
  # Extract the first two numbers (two digits) after [TG(
  df_filtered$TG_length <- as.numeric(str_extract(df_filtered$lipid, "(?<=\\[TG\\()\\d{2}"))
  
  # Aggregate the data
  df_summarized <- df_filtered %>% 
    group_by(TG_length) %>% 
    summarize(total_mean1 = sum(mean1), total_mean2 = sum(mean2))
  
  # Calculate percentages
  df_summarized$percent_mean1 <- df_summarized$total_mean1 / sum(df_summarized$total_mean1)
  df_summarized$percent_mean2 <- df_summarized$total_mean2 / sum(df_summarized$total_mean2)
  
  directory_name <- "TG_DG_lengths"
  file_name <- paste0("TG_", filename, ".csv")
  full_path <- file.path(directory_name, file_name)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
  }
  
  # Save the summarized data frame to a CSV file
  write.csv(df_summarized, full_path, row.names = FALSE)
  
  # Create the line plot
  p<- ggplot(df_summarized, aes(x = TG_length)) +
    geom_line(aes(y = percent_mean1), color = "blue") +
    geom_line(aes(y = percent_mean2), color = "red") +
    labs(title = Title1, subtitle = Title2, x = "TG Length", y = "Percentage") +
    theme_minimal()
  
  ggsave(filename, plot = p, width = 12, height = 8, units = "in")
  
  
}




plot_DG_distribution_individual <- function(df, Title1, Title2,filename) {
  # Filter the dataframe
  df_filtered <- df %>% filter(type == "DAG") #%>% filter(FDR < 0.1)
  
  # Extract the first two numbers (two digits) after DG()
  # df_filtered$TG_length <- as.numeric(str_extract(df_filtered$lipid, "(?<=\\DG\\()(dO-|O-)?\\d{2}"))
  # Step 1: Extract the relevant part of the string
  df_filtered$extracted_string <- str_extract(df_filtered$lipid, "(dO-|O-)?\\d{2}")
  
  # Step 2: Remove "dO-" or "O-" and convert to numeric
  df_filtered$DG_length <- as.numeric(str_replace(df_filtered$extracted_string, "dO-|O-", ""))
  
  
  # Aggregate the data
  df_summarized <- df_filtered %>% 
    group_by(DG_length) %>% 
    summarize(total_mean1 = sum(mean1), total_mean2 = sum(mean2))
  
  # Calculate percentages
  df_summarized$percent_mean1 <- df_summarized$total_mean1 / sum(df_summarized$total_mean1)
  df_summarized$percent_mean2 <- df_summarized$total_mean2 / sum(df_summarized$total_mean2)
  directory_name <- "TG_DG_lengths"
  file_name <- paste0("DG_", filename, ".csv")
  full_path <- file.path(directory_name, file_name)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
  }
  
  # Save the summarized data frame to a CSV file
  write.csv(df_summarized, full_path, row.names = FALSE)
  # Create the line plot
  p<- ggplot(df_summarized, aes(x = DG_length)) +
    geom_line(aes(y = percent_mean1), color = "blue") +
    geom_line(aes(y = percent_mean2), color = "red") +
    labs(title = Title1, subtitle = Title2, x = "DG Length", y = "Percentage") +
    theme_minimal()
  
  ggsave(filename, plot = p, width = 12, height = 8, units = "in")
  
  
}







plot_TG_distribution_individual_chain_length <- function(df, Title1, Title2,filename) {
  # Filter the dataframe
  df_filtered <- df %>% filter(type == "TAG")# %>% filter(FDR < 0.1)
  
  # Extract the first two numbers (two digits) after [TG(
  # Extract TG_length
  df_filtered <- df_filtered %>%
    mutate(TG_length = str_extract(lipid, "(?<=_FA).+$"))
  
  # Modify TG_length values based on conditions
  df_filtered <- df_filtered %>%
    mutate(TG_length = ifelse(TG_length %in% c("14:0 | [TG(66:7)]_FA14:0", "14:0 | [TG(66:7)]_FA14:0_2"), "14:0", TG_length))
  
  # Aggregate the data
  df_summarized <- df_filtered %>% 
    group_by(TG_length) %>% 
    summarize(total_mean1 = sum(mean1), total_mean2 = sum(mean2))
  
  # Calculate percentages
  df_summarized$percent_mean1 <- df_summarized$total_mean1 / sum(df_summarized$total_mean1)
  df_summarized$percent_mean2 <- df_summarized$total_mean2 / sum(df_summarized$total_mean2)
  
  directory_name <- "TG_DG_lengths_FA"
  file_name <- paste0("TG_", filename, ".csv")
  full_path <- file.path(directory_name, file_name)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
  }
  
  # Save the summarized data frame to a CSV file
  write.csv(df_summarized, full_path, row.names = FALSE)
  
  # Create the line plot
  p<- ggplot(df_summarized, aes(x = TG_length)) +
    geom_line(aes(y = percent_mean1), color = "blue") +
    geom_line(aes(y = percent_mean2), color = "red") +
    labs(title = Title1, subtitle = Title2, x = "TG Length", y = "Percentage") +
    theme_minimal()
  
  # ggsave(filename, plot = p, width = 12, height = 8, units = "in")
  
  
}


plot_DG_distribution_individual_chain_length <- function(df, Title1, Title2,filename) {
  # Filter the dataframe
  # df_filtered <- df %>% filter(type == "DAG") #%>% filter(FDR < 0.1)
  # 
  # # Extract the first two numbers (two digits) after DG()
  # # df_filtered$TG_length <- as.numeric(str_extract(df_filtered$lipid, "(?<=\\DG\\()(dO-|O-)?\\d{2}"))
  # # Step 1: Extract the relevant part of the string
  # df_filtered$extracted_string <- str_extract(df_filtered$lipid, "(dO-|O-)?\\d{2}")
  # 
  # # Step 2: Remove "dO-" or "O-" and convert to numeric
  # df_filtered$DG_length <- as.numeric(str_replace(df_filtered$extracted_string, "dO-|O-", ""))
  # 
  # 
  
  df_filtered <- df %>%
    filter(type == "DAG") %>% #filter(FDR < 0.1) %>%
    mutate(DG_length = str_extract(lipid, "(?<=_C).+$")) %>%
    arrange(DG_length)
  
  
  # Aggregate the data
  df_summarized <- df_filtered %>% 
    group_by(DG_length) %>% 
    summarize(total_mean1 = sum(mean1), total_mean2 = sum(mean2))
  
  # Calculate percentages
  df_summarized$percent_mean1 <- df_summarized$total_mean1 / sum(df_summarized$total_mean1)
  df_summarized$percent_mean2 <- df_summarized$total_mean2 / sum(df_summarized$total_mean2)
  directory_name <- "TG_DG_lengths_FA"
  file_name <- paste0("DG_", filename, ".csv")
  full_path <- file.path(directory_name, file_name)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
  }
  
  # Save the summarized data frame to a CSV file
  write.csv(df_summarized, full_path, row.names = FALSE)
  # Create the line plot
  p<- ggplot(df_summarized, aes(x = DG_length)) +
    geom_line(aes(y = percent_mean1), color = "blue") +
    geom_line(aes(y = percent_mean2), color = "red") +
    labs(title = Title1, subtitle = Title2, x = "DG Length", y = "Percentage") +
    theme_minimal()
  
  # ggsave(filename, plot = p, width = 12, height = 8, units = "in")
  
  
}


















file_list = list.files(path="processed_results_2", pattern=NULL, all.files=FALSE,
                       full.names=FALSE)
# 
# filename <- paste(dir_name, "/", Title1, Title2, "_TAG_DAG_Histograms.pdf", sep = "")
# ggsave(filename, plot = p, width = 12, height = 8, units = "in")

filename <-paste0("TG_distribution_length/",Title1,Title2,"ALL_LIPIDS.pdf",sep="")
df <- excel_file
for (jj in file_list){
  excel_file <- read_csv(paste0("processed_results_2/",jj,sep=""))


  # title_for_plot <- gsub(":", "_",excel_file$Title[1])
  Title1 <- gsub("[:| ]", "_", excel_file$Title_1[1])
  Title2 <- gsub("[:| ]", "_", excel_file$Title_2[1])
  Title1 <- gsub("__", "_", Title1)
  Title2 <- gsub("__", "_", Title2) 

  plot_TG_distribution_individual(excel_file, Title1, Title2, paste0("TG_distribution_length/",Title1,Title2,"ALL_LIPIDS.pdf",sep=""))
  
  plot_DG_distribution_individual(excel_file, Title1, Title2, paste0("DG_distribution_length/",Title1,Title2,"ALL_LIPIDS.pdf",sep=""))
  plot_DG_distribution_individual_chain_length(excel_file, Title1, Title2, paste0("DG_distribution_length/",Title1,Title2,"ALL_LIPIDS.pdf",sep=""))
  plot_TG_distribution_individual_chain_length(excel_file, Title1, Title2, paste0("TG_distribution_length/",Title1,Title2,"ALL_LIPIDS.pdf",sep=""))
  plot_FA_distribution_individual(excel_file, Title1, Title2, paste0("FA_distribution/",Title1,Title2,"ALL_LIPIDS.pdf",sep=""))
  plot_FA_distribution_with_saturation(excel_file, Title1, Title2, paste0(Title1,Title2,"Saturation_ALL_LIPIDS.pdf",sep=""))
  
}



# Example usage
# plot_TG_distribution_individual(your_dataframe, "Title 1", "Title 2")

