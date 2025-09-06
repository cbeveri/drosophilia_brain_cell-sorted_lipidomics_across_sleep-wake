# ====================  FA Chain‑Length Stacked‑Bar Pipeline  ====================
# Author: Connor Beveridge, PhD
# Last edit: 2025‑05‑07
#
#  • Absolute totals  (BAR_FA_No_Normal)     – now with centred % labels
#  • Normalised (%)   (BAR_FA)
# ==============================================================================

# ---------- Libraries ----------
library(tidyverse)
library(scales)

# ---------- Folders ----------
dirs <- c("BAR_FA", "BAR_FA_No_Normal", "processed_results_2", "plots")
walk(dirs[!dir.exists(dirs)], dir.create)

# ---------- Helper ------------------------------------------------------------
chain_length_classes <- c("Short-Chain","Medium-Chain","Long-Chain",
                          "Very Long-Chain","Ultra Long-Chain")
chain_length_colors  <- c("#a6cee3","#1f78b4","#b2df8a",
                          "#33a02c","#fb9a99")

categorize_chain <- function(fa_type){
  n <- as.numeric(sub(".*\\((\\d+):\\d+\\).*", "\\1", fa_type))
  if (is.na(n))                 return(NA)
  if (n <  6)                   return("Short-Chain")
  if (n < 13)                   return("Medium-Chain")
  if (n < 22)                   return("Long-Chain")
  if (n < 30)                   return("Very Long-Chain")
  "Ultra Long-Chain"
}

# ---------- Main function -----------------------------------------------------
FA_plot_lengths <- function(df, Title1, Title2){
  # -- keep only fatty acids ---------------------------------------------------
  df <- df %>% filter(type == "FA")
  if (nrow(df) == 0) return(invisible())        # nothing to plot
  
  # -- columns for the two comparison groups ----------------------------------
  len1      <- as.numeric(df$Length1[1])
  len2      <- as.numeric(df$Length2[1])
  all_cols  <- colnames(df)
  start_idx <- which(all_cols == "type") + 1
  cols1     <- all_cols[start_idx:(start_idx + len1 - 1)]
  cols2     <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
  blank_col <- tail(all_cols, 1)
  
  # -- blank correction --------------------------------------------------------
  df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
  df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
  
  # -- ensure numerics ---------------------------------------------------------
  df[cols1] <- lapply(df[cols1], \(x) as.numeric(as.character(x)))
  df[cols2] <- lapply(df[cols2], \(x) as.numeric(as.character(x)))
  df[is.na(df)] <- 0
  
  # -- means + chain‑length category ------------------------------------------
  df <- df %>% mutate(
    mean1        = rowMeans(across(all_of(cols1))),
    mean2        = rowMeans(across(all_of(cols2))),
    chain_length = factor(map_chr(lipid, categorize_chain),
                          levels = chain_length_classes)
  )
  
  # -- aggregate ---------------------------------------------------------------
  df_all <- df %>% group_by(chain_length) %>% summarise(
    sum_mean1 = sum(mean1),
    sum_mean2 = sum(mean2), .groups = "keep")
  
  total1 <- sum(df_all$sum_mean1)
  total2 <- sum(df_all$sum_mean2)
  
  df_all <- df_all %>% mutate(
    pct1 = 100 * sum_mean1 / total1,
    pct2 = 100 * sum_mean2 / total2)
  
  # ---------------- absolute plot (with centred % labels) ---------------------
  orig_data <- bind_rows(
    df_all %>% transmute(group = Title1, chain_length, value = sum_mean1),
    df_all %>% transmute(group = Title2, chain_length, value = sum_mean2)) %>%
    group_by(group) %>%
    mutate(percent = 100 * value / sum(value)) %>%
    ungroup()
  
  p_abs <- ggplot(orig_data, aes(group, value, fill = chain_length)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", percent)),
              position = position_stack(vjust = 0.5),
              size = 3, colour = "black") +
    scale_fill_manual(values = setNames(chain_length_colors,
                                        chain_length_classes),
                      name = "Chain Length") +
    ylab("Total lipid content") +
    theme_classic() +
    labs(title = sprintf("%s vs %s (FA Chain Length)", Title1, Title2))
  
  ggsave(sprintf("BAR_FA_No_Normal/%s_vs_%s_STACKED_FA_ChainLength.pdf",
                 Title1, Title2), p_abs, width = 8, height = 6)
  
  # ---------------- normalised (%) plot ---------------------------------------
  norm_data  <- bind_rows(
    df_all %>% transmute(group = Title1, chain_length, percent = pct1),
    df_all %>% transmute(group = Title2, chain_length, percent = pct2))
  
  label_norm <- norm_data %>% group_by(group) %>% arrange(desc(chain_length)) %>%
    mutate(position = cumsum(percent) - 0.5 * percent)
  
  p_pct <- ggplot(norm_data, aes(group, percent, fill = chain_length)) +
    geom_col() +
    geom_text(data = label_norm,
              aes(y = position,
                  label = sprintf("%.1f%%", percent)),
              size = 3, colour = "black") +
    scale_fill_manual(values = setNames(chain_length_colors,
                                        chain_length_classes),
                      name = "Chain Length") +
    scale_y_continuous(labels = \(x) paste0(x, "%"), limits = c(0, 100)) +
    ylab("Proportion of total lipid content (%)") +
    theme_classic() +
    labs(title = sprintf("%s vs %s (FA Chain Length – Normalised)",
                         Title1, Title2))
  
  ggsave(sprintf("BAR_FA/%s_vs_%s_NORMALISED_STACKED_FA_ChainLength.pdf",
                 Title1, Title2), p_pct, width = 8, height = 6)
  
  # ---------------- CSV exports ----------------------------------------------
  write_csv(norm_data %>% arrange(group, chain_length),
            sprintf("BAR_FA/%s_vs_%s_FA_ChainLength_Percentages.csv",
                    Title1, Title2))
  write_csv(orig_data %>% select(group, chain_length, value, percent) %>%
              arrange(group, chain_length),
            sprintf("BAR_FA_No_Normal/%s_vs_%s_FA_ChainLength_Totals.csv",
                    Title1, Title2))
}

# ---------- Batch run ----------
file_list <- list.files("results", full.names = FALSE)
for (f in file_list){
  if (!grepl("full", f, ignore.case = TRUE)) next
  dat <- read_csv(file.path("results", f), show_col_types = FALSE)
  T1  <- dat$Title_1[1] %>% str_replace_all("[:| ]","_") %>% str_replace_all("__","_")
  T2  <- dat$Title_2[1] %>% str_replace_all("[:| ]","_") %>% str_replace_all("__","_")
  FA_plot_lengths(dat, T1, T2)
}




# # Load necessary libraries
# library(tidyverse)
# library(scales)
# 
# # Create necessary directories
# if (!dir.exists("BAR_FA")) {
#   dir.create("BAR_FA")
# }
# if (!dir.exists("processed_results_2")) {
#   dir.create("processed_results_2")
# }
# if (!dir.exists("plots")) {
#   dir.create("plots")
# }
# 
# # Stacked Bar Chart Function for FA Chain Lengths
# FA_plot_lengths <- function(df, Title1, Title2) {
#   # Define FA chain length categories
#   chain_length_classes <- c("Short-Chain", "Medium-Chain", "Long-Chain", "Very Long-Chain", "Ultra Long-Chain")
#   # Define color palette for chain length categories
#   chain_length_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99")
#   
#   # Function to categorize fatty acid chain lengths based on the number before the colon (e.g., 16 in FA(16:0))
#   categorize_chain_length <- function(fa_type) {
#     chain_length <- as.numeric(sub(".*\\((\\d+):\\d+\\).*", "\\1", fa_type))
#     if (is.na(chain_length)) return(NA)
#     if (chain_length < 6) {
#       return("Short-Chain")
#     } else if (chain_length < 13) {
#       return("Medium-Chain")
#     } else if (chain_length < 22) {
#       return("Long-Chain")
#     } else if (chain_length < 30) {
#       return("Very Long-Chain")
#     } else {
#       return("Ultra Long-Chain")
#     }
#   }
#   
#   # Filter FA lipids
#   df <- df %>% filter(type == "FA")
#   
#   # Extract group columns
#   len1 <- as.numeric(df$Length1[1])
#   len2 <- as.numeric(df$Length2[1])
#   all_cols <- colnames(df)
#   start_idx <- which(all_cols == "type") + 1
#   cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   blank_col <- all_cols[length(all_cols)]
#   
#   # Apply blank correction
#   df[cols1] <- pmax(df[cols1] - df[[blank_col]], 0)
#   df[cols2] <- pmax(df[cols2] - df[[blank_col]], 0)
#   
#   # Convert to numeric if necessary
#   df[cols1] <- apply(df[cols1], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   df[cols2] <- apply(df[cols2], 2, function(x) ifelse(!is.na(as.numeric(as.character(x))), x, 0))
#   
#   # Calculate means and categorize chain lengths
#   df <- df %>%
#     mutate(mean1 = rowMeans(select(., one_of(cols1))),
#            mean2 = rowMeans(select(., one_of(cols2))),
#            chain_length = factor(sapply(lipid, categorize_chain_length),
#                                  levels = chain_length_classes))
#   
#   # Group by chain length and sum
#   df_all <- df %>%
#     group_by(chain_length) %>%
#     summarise(sum_mean1 = sum(mean1),
#               sum_mean2 = sum(mean2), .groups = "keep")
#   
#   # Calculate totals for normalization
#   total_mean1 <- sum(df_all$sum_mean1)
#   total_mean2 <- sum(df_all$sum_mean2)
#   
#   # Calculate percentages
#   df_all <- df_all %>%
#     mutate(
#       percent_mean1 = sum_mean1 / total_mean1 * 100,
#       percent_mean2 = sum_mean2 / total_mean2 * 100
#     )
#   
#   # Create data for the new normalized plot
#   plot_data <- bind_rows(
#     df_all %>% 
#       select(chain_length, percent = percent_mean1) %>%
#       mutate(group = Title1),
#     df_all %>% 
#       select(chain_length, percent = percent_mean2) %>%
#       mutate(group = Title2)
#   )
#   
#   # Original non-normalized stacked bar plot (for backward compatibility)
#   orig_data <- bind_rows(
#     df_all %>% 
#       select(chain_length, value = sum_mean1) %>%
#       mutate(group = Title1),
#     df_all %>% 
#       select(chain_length, value = sum_mean2) %>%
#       mutate(group = Title2)
#   )
#   
#   plot_fa_original <- ggplot(orig_data, aes(x = group, y = value, fill = chain_length)) +
#     geom_col() +
#     theme_classic() +
#     ylab("Total lipid content") +
#     scale_fill_manual(values = setNames(chain_length_colors, chain_length_classes), name = "Chain Length") +
#     labs(title = paste(Title1, "vs", Title2, "(FA Chain Length)"))
# 
#   ggsave(filename = paste0("BAR_FA_No_Normal/", Title1, "_vs_", Title2, "_STACKED_FA_ChainLength.pdf"),
#          plot = plot_fa_original)
#   
#   # Create the normalized stacked bar plot
#   plot_fa_normalized <- ggplot(plot_data, aes(x = group, y = percent, fill = chain_length)) +
#     geom_col(position = "stack") +
#     theme_classic() +
#     ylab("Proportion of total lipid content (%)") +
#     scale_fill_manual(values = setNames(chain_length_colors, chain_length_classes), name = "Chain Length") +
#     labs(title = paste(Title1, "vs", Title2, "(FA Chain Length - Normalized)")) +
#     scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100))
#   
#   # Calculate label positions for each segment
#   label_data <- plot_data %>%
#     group_by(group) %>%
#     arrange(group, desc(chain_length)) %>%
#     mutate(
#       # Calculate position for each segment
#       position = cumsum(percent) - 0.5 * percent
#     )
#   
#   # Add percentage labels to all segments
#   plot_fa_normalized <- plot_fa_normalized +
#     geom_text(
#       data = label_data,
#       aes(y = position, label = sprintf("%.1f%%", percent)),
#       color = "black", size = 3
#     )
#   
#   # Save the normalized plot
#   ggsave(filename = paste0("BAR_FA/", Title1, "_vs_", Title2, "_NORMALIZED_STACKED_FA_ChainLength.pdf"),
#          plot = plot_fa_normalized, width = 8, height = 6)
#   
#   # Save the percentage data to a CSV file
#   write.csv(
#     plot_data %>% select(group, chain_length, percent) %>% arrange(group, chain_length),
#     file = paste0("BAR_FA/", Title1, "_vs_", Title2, "_FA_ChainLength_Percentages.csv"),
#     row.names = FALSE
#   )
# }
# 
# # Main loop over result files
# file_list <- list.files(path = "results", pattern = NULL, full.names = FALSE)
# 
# for (jj in file_list) {
#   if (!grepl("full", jj, ignore.case = TRUE)) next
#   
#   excel_file <- read_csv(paste0("results/", jj))
#   Title1 <- gsub("[:| ]", "_", excel_file$Title_1[1])
#   Title2 <- gsub("[:| ]", "_", excel_file$Title_2[1])
#   Title1 <- gsub("__", "_", Title1)
#   Title2 <- gsub("__", "_", Title2)
#   
#   FA_plot_lengths(excel_file, Title1, Title2)
# }

