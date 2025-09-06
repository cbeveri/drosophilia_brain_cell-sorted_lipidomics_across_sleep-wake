library(tidyverse)
library(scales)
library(gridExtra)

# Create output folders if they don't exist
dir.create("Pie_charts", showWarnings = FALSE)
dir.create("Pie_charts/no_FA", showWarnings = FALSE)

# Define the function to plot and save pie charts
plot_pie_charts_updated <- function(df, Title1, Title2, file_tag) {
  # Define lipid classes and colors
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG", "DAG", "TAG | DAG", "DAG | CE", "TAG | DAG | CE")
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
  lipid_class_colors <- setNames(alpha(lipid_colors, 0.5), lipid_classes)
  
  # Remove 'CAR' entries that contain 'QUAL'
  df <- df %>% filter(!(type == "CAR" & str_detect(lipid, "QUAL")))
  
  # Identify column indices
  len1 <- as.numeric(df$Length1[1])
  len2 <- as.numeric(df$Length2[1])
  all_cols <- colnames(df)
  start_idx <- which(all_cols == "type") + 1
  cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
  cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
  blank_col <- all_cols[length(all_cols)]
  
  # Blank subtraction and value cleaning
  df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
  df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
  
  # Mean calculation
  df <- df %>%
    mutate(mean1 = rowMeans(select(., all_of(cols1)), na.rm = TRUE),
           mean2 = rowMeans(select(., all_of(cols2)), na.rm = TRUE))
  
  # Summarize by lipid class
  df_summarized <- df %>%
    group_by(type) %>%
    summarise(sum_mean1 = sum(mean1, na.rm = TRUE),
              sum_mean2 = sum(mean2, na.rm = TRUE), .groups = "drop")
  
  # Save full intensity CSV
  write.csv(df_summarized, file = paste0("Pie_charts/", file_tag, "_intensities.csv"), row.names = FALSE)
  
  # Plot pie charts (with FA)
  pie1 <- ggplot(df_summarized, aes(x = "", y = sum_mean1, fill = type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") + theme_void() +
    scale_fill_manual(values = lipid_class_colors) + labs(title = Title1)
  
  pie2 <- ggplot(df_summarized, aes(x = "", y = sum_mean2, fill = type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") + theme_void() +
    scale_fill_manual(values = lipid_class_colors) + labs(title = Title2)
  
  svg(paste0("Pie_charts/", file_tag, ".svg"), width = 10, height = 5)
  grid.arrange(pie1, pie2, ncol = 2)
  dev.off()
  
  # Plot pie charts (excluding FA)
  df_no_fa <- df_summarized %>% filter(type != "FA")
  write.csv(df_no_fa, file = paste0("Pie_charts/no_FA/", file_tag, "_no_FA_intensities.csv"), row.names = FALSE)
  
  pie1_no_fa <- ggplot(df_no_fa, aes(x = "", y = sum_mean1, fill = type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") + theme_void() +
    scale_fill_manual(values = lipid_class_colors) + labs(title = paste0(Title1, " (No FA)"))
  
  pie2_no_fa <- ggplot(df_no_fa, aes(x = "", y = sum_mean2, fill = type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") + theme_void() +
    scale_fill_manual(values = lipid_class_colors) + labs(title = paste0(Title2, " (No FA)"))
  
  svg(paste0("Pie_charts/no_FA/", file_tag, "_no_FA.svg"), width = 10, height = 5)
  grid.arrange(pie1_no_fa, pie2_no_fa, ncol = 2)
  dev.off()
}

# =============================
# === Process all CSV files ===
# =============================

file_list <- list.files("results", pattern = ".csv$", full.names = TRUE)

for (file in file_list) {
  excel_file <- read_csv(file)
  if (!grepl("full", file, ignore.case = TRUE)) next
  
  Title1 <- gsub("[:| ]", "_", excel_file$Title_1[1])
  Title2 <- gsub("[:| ]", "_", excel_file$Title_2[1])
  Title1 <- gsub("__", "_", Title1)
  Title2 <- gsub("__", "_", Title2)
  file_tag <- paste0(Title1, "_vs_", Title2)
  
  cat("Processing:", file_tag, "\n")
  
  plot_pie_charts_updated(excel_file, Title1, Title2, file_tag)
}

