# ────────────────────────────────────────────────────────────────────────────────
# New FA Bar Plots - Specific Chain Lengths with Error Bars
# ────────────────────────────────────────────────────────────────────────────────
# Creates bar plots for FA(2:X), FA(14:X), FA(15:X), FA(16:X) with all saturation levels
# Two plot types: multiple groups combined, and pairwise A vs B comparisons
# Both absolute and normalized versions with standard error bars
# ────────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(scales)

# Target chain lengths
target_chain_lengths <- c(2, 14, 15, 16)

# Create directory structure
create_fa_directories <- function() {
  dir.create("new_FA_bar", showWarnings = FALSE)
  dir.create("new_FA_bar/multiple_groups", showWarnings = FALSE)
  dir.create("new_FA_bar/multiple_groups/absolute", showWarnings = FALSE)
  dir.create("new_FA_bar/multiple_groups/normalized", showWarnings = FALSE)
  dir.create("new_FA_bar/pairwise", showWarnings = FALSE)
  dir.create("new_FA_bar/pairwise/absolute", showWarnings = FALSE)
  dir.create("new_FA_bar/pairwise/normalized", showWarnings = FALSE)
}

# Extract FA information from lipid names
extract_fa_info <- function(lipid_names) {
  fa_pattern <- "FA\\((\\d+):(\\d+)\\)"
  matches <- str_match(lipid_names, fa_pattern)
  
  data.frame(
    lipid = lipid_names,
    chain_length = as.numeric(matches[,2]),
    saturation = as.numeric(matches[,3]),
    fa_type = paste0("FA(", matches[,2], ":", matches[,3], ")"),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(chain_length), chain_length %in% target_chain_lengths)
}

# Process single file for FA data
process_fa_file <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(tibble())
  }
  
  dat <- read_csv(file_path, show_col_types = FALSE) %>%
    filter(type == "FA")
  
  if (nrow(dat) == 0) return(tibble())
  
  # Extract group info
  g1 <- dat$Title_1[1]
  g2 <- dat$Title_2[1]
  len1 <- as.numeric(dat$Length1[1])
  len2 <- as.numeric(dat$Length2[1])
  
  # Identify columns
  all_cols <- colnames(dat)
  start_idx <- which(all_cols == "type") + 1
  c1 <- all_cols[start_idx:(start_idx + len1 - 1)]
  c2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
  blank_col <- tail(all_cols, 1)
  
  # Subtract blank and ensure non-negative
  dat[c1] <- pmax(dat[c1] - dat[[blank_col]], 0)
  dat[c2] <- pmax(dat[c2] - dat[[blank_col]], 0)
  
  # Convert to numeric
  dat[c(c1, c2)] <- lapply(dat[c(c1, c2)], function(x) as.numeric(as.character(x)))
  dat[is.na(dat)] <- 0
  
  # Extract FA info
  fa_info <- extract_fa_info(dat$lipid)
  if (nrow(fa_info) == 0) return(tibble())
  
  # Combine data with FA info
  dat_fa <- dat %>%
    select(lipid, all_of(c(c1, c2))) %>%
    inner_join(fa_info, by = "lipid")
  
  # Reshape to long format
  result <- dat_fa %>%
    pivot_longer(cols = all_of(c(c1, c2)), names_to = "sample", values_to = "value") %>%
    mutate(
      group = ifelse(sample %in% c1, g1, g2),
      file_source = basename(file_path)
    ) %>%
    select(fa_type, chain_length, saturation, group, sample, value, file_source)
  
  return(result)
}

# Create multiple groups plot
create_multiple_groups_fa_plot <- function(files, save_name) {
  
  # Process all files
  all_data <- map_dfr(files, ~process_fa_file(file.path("results", 
                                                        ifelse(endsWith(.x, ".csv"), .x, paste0(.x, ".csv")))))
  
  if (nrow(all_data) == 0) {
    warning("No FA data found in files")
    return(invisible())
  }
  
  # Calculate summary statistics
  summary_data <- all_data %>%
    group_by(fa_type, chain_length, saturation, group) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      se_value = ifelse(n() > 1, sd(value, na.rm = TRUE) / sqrt(n()), 0),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    arrange(chain_length, saturation)
  
  # Create factor with ordered levels
  summary_data$fa_type <- factor(summary_data$fa_type, 
                                 levels = unique(summary_data$fa_type[order(summary_data$chain_length, summary_data$saturation)]))
  
  # Identify FA types where all values are 0
  zero_fa_types <- summary_data %>%
    group_by(fa_type) %>%
    summarise(all_zero = all(mean_value == 0), .groups = "drop") %>%
    filter(all_zero) %>%
    pull(fa_type)
  
  # Create "all" and "removed_0" versions
  summary_data_all <- summary_data
  summary_data_removed <- summary_data %>% filter(!fa_type %in% zero_fa_types)
  
  # Function to create plots
  create_plots <- function(data, suffix, plot_title_suffix) {
    if (nrow(data) == 0) return(list())
    
    # ABSOLUTE PLOT
    p_abs <- ggplot(data, aes(x = fa_type, y = mean_value, fill = group)) +
      geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
      geom_errorbar(aes(ymin = pmax(0, mean_value - se_value), ymax = mean_value + se_value),
                    position = position_dodge(width = 0.8), width = 0.4) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title = paste("FA Chain Lengths - Multiple Groups (Absolute)", plot_title_suffix, "-", save_name),
           x = "Fatty Acid",
           y = "Absolute Value",
           fill = "Group") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
    
    # NORMALIZED PLOT (each FA type normalized to its max across all groups)
    norm_data <- data %>%
      group_by(fa_type) %>%
      mutate(
        max_value = max(mean_value, na.rm = TRUE),
        norm_mean = ifelse(max_value > 0, (mean_value / max_value) * 100, 0),
        norm_se = ifelse(max_value > 0, (se_value / max_value) * 100, 0)
      ) %>%
      ungroup()
    
    p_norm <- ggplot(norm_data, aes(x = fa_type, y = norm_mean, fill = group)) +
      geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
      geom_errorbar(aes(ymin = pmax(0, norm_mean - norm_se), ymax = pmin(100, norm_mean + norm_se)),
                    position = position_dodge(width = 0.8), width = 0.4) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title = paste("FA Chain Lengths - Multiple Groups (Normalized)", plot_title_suffix, "-", save_name),
           x = "Fatty Acid",
           y = "Relative % (Max = 100%)",
           fill = "Group") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 105))
    
    # Save plots
    ggsave(file.path("new_FA_bar/multiple_groups/absolute", paste0(save_name, "_absolute_", suffix, ".pdf")), 
           p_abs, width = 14, height = 8)
    ggsave(file.path("new_FA_bar/multiple_groups/absolute", paste0(save_name, "_absolute_", suffix, ".png")), 
           p_abs, width = 14, height = 8, dpi = 300)
    
    ggsave(file.path("new_FA_bar/multiple_groups/normalized", paste0(save_name, "_normalized_", suffix, ".pdf")), 
           p_norm, width = 14, height = 8)
    ggsave(file.path("new_FA_bar/multiple_groups/normalized", paste0(save_name, "_normalized_", suffix, ".png")), 
           p_norm, width = 14, height = 8, dpi = 300)
    
    return(list(absolute = p_abs, normalized = p_norm))
  }
  
  # Create both versions
  plots_all <- create_plots(summary_data_all, "all", "(All)")
  plots_removed <- create_plots(summary_data_removed, "removed_0", "(Non-zero)")
  
  cat("Multiple groups FA plots saved for:", save_name, "\n")
  cat("- Removed", length(zero_fa_types), "FA types with all zero values\n")
  
  return(invisible(list(
    all = plots_all, 
    removed_0 = plots_removed, 
    data = summary_data_all,
    zero_fa_types = zero_fa_types
  )))
}

# Create pairwise A vs B plots
create_pairwise_fa_plots <- function() {
  
  # Get all full files
  all_files <- list.files("results", pattern = "_full\\.csv$", full.names = FALSE)
  
  for (file in all_files) {
    file_base <- gsub("\\.csv$", "", file)
    
    # Process single file
    data <- process_fa_file(file.path("results", file))
    
    if (nrow(data) == 0) {
      cat("No FA data in:", file, "\n")
      next
    }
    
    # Calculate summary statistics
    summary_data <- data %>%
      group_by(fa_type, chain_length, saturation, group) %>%
      summarise(
        mean_value = mean(value, na.rm = TRUE),
        se_value = ifelse(n() > 1, sd(value, na.rm = TRUE) / sqrt(n()), 0),
        n_samples = n(),
        .groups = "drop"
      ) %>%
      arrange(chain_length, saturation)
    
    # Create factor with ordered levels
    summary_data$fa_type <- factor(summary_data$fa_type, 
                                   levels = unique(summary_data$fa_type[order(summary_data$chain_length, summary_data$saturation)]))
    
    # Identify FA types where all values are 0
    zero_fa_types <- summary_data %>%
      group_by(fa_type) %>%
      summarise(all_zero = all(mean_value == 0), .groups = "drop") %>%
      filter(all_zero) %>%
      pull(fa_type)
    
    # Create "all" and "removed_0" versions
    summary_data_all <- summary_data
    summary_data_removed <- summary_data %>% filter(!fa_type %in% zero_fa_types)
    
    # Function to create plots
    create_pairwise_plots <- function(plot_data, suffix, plot_title_suffix) {
      if (nrow(plot_data) == 0) return(list())
      
      # ABSOLUTE PLOT
      p_abs <- ggplot(plot_data, aes(x = fa_type, y = mean_value, fill = group)) +
        geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
        geom_errorbar(aes(ymin = pmax(0, mean_value - se_value), ymax = mean_value + se_value),
                      position = position_dodge(width = 0.8), width = 0.4) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom") +
        labs(title = paste("FA Chain Lengths - Pairwise (Absolute)", plot_title_suffix),
             subtitle = file_base,
             x = "Fatty Acid",
             y = "Absolute Value",
             fill = "Group") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
      
      # NORMALIZED PLOT
      norm_data <- plot_data %>%
        group_by(fa_type) %>%
        mutate(
          max_value = max(mean_value, na.rm = TRUE),
          norm_mean = ifelse(max_value > 0, (mean_value / max_value) * 100, 0),
          norm_se = ifelse(max_value > 0, (se_value / max_value) * 100, 0)
        ) %>%
        ungroup()
      
      p_norm <- ggplot(norm_data, aes(x = fa_type, y = norm_mean, fill = group)) +
        geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
        geom_errorbar(aes(ymin = pmax(0, norm_mean - norm_se), ymax = pmin(100, norm_mean + norm_se)),
                      position = position_dodge(width = 0.8), width = 0.4) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom") +
        labs(title = paste("FA Chain Lengths - Pairwise (Normalized)", plot_title_suffix),
             subtitle = file_base,
             x = "Fatty Acid",
             y = "Relative % (Max = 100%)",
             fill = "Group") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 105))
      
      # Save plots
      safe_name <- gsub("[^A-Za-z0-9_-]", "_", file_base)
      
      ggsave(file.path("new_FA_bar/pairwise/absolute", paste0(safe_name, "_absolute_", suffix, ".pdf")), 
             p_abs, width = 12, height = 8)
      ggsave(file.path("new_FA_bar/pairwise/absolute", paste0(safe_name, "_absolute_", suffix, ".png")), 
             p_abs, width = 12, height = 8, dpi = 300)
      
      ggsave(file.path("new_FA_bar/pairwise/normalized", paste0(safe_name, "_normalized_", suffix, ".pdf")), 
             p_norm, width = 12, height = 8)
      ggsave(file.path("new_FA_bar/pairwise/normalized", paste0(safe_name, "_normalized_", suffix, ".png")), 
             p_norm, width = 12, height = 8, dpi = 300)
      
      return(list(absolute = p_abs, normalized = p_norm))
    }
    
    # Create both versions
    plots_all <- create_pairwise_plots(summary_data_all, "all", "(All)")
    plots_removed <- create_pairwise_plots(summary_data_removed, "removed_0", "(Non-zero)")
    
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", file_base)
    cat("Pairwise FA plots saved for:", safe_name)
    if (length(zero_fa_types) > 0) {
      cat(" - Removed", length(zero_fa_types), "zero FA types")
    }
    cat("\n")
  }
}

# ────────────────────────────────────────────────────────────────────────────────
# MAIN EXECUTION
# ────────────────────────────────────────────────────────────────────────────────

# Create directories
create_fa_directories()

# 1. MULTIPLE GROUPS PLOTS - Sex-specific files
glia_files <- c(
  "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ Female vs Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ FemaleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ Female vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ FemaleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ maleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ maleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ maleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ Female vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ FemaleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ maleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ maleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
  "Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ maleBlank1_full"
)

neuron_files <- c(
  "Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ Female vs Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ FemaleBlank1_full",
  "Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ Female vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
  "Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ FemaleBlank1_full",
  "Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ maleBlank1_full",
  "Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ maleBlank1_full",
  "Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ Female vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
  "Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ FemaleBlank1_full",
  "Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ maleBlank1_full",
  "Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ maleBlank1_full",
  "Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
  "Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ maleBlank1_full"
)

supernatant_files <- c(
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_FemalevsSampleType_Supernatant60BrainSleepCondition_GoingtosleepSex_FemaleBlank1_full",
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_FemalevsSampleType_Supernatant60BrainSleepCondition_WakingupSex_FemaleBlank1_full",
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_DeprivedSex_FemaleBlank1_full",
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_DeprivedSex_maleBlank1_full",
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_GoingtosleepSex_maleBlank1_full",
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_WakingupSex_maleBlank1_full"
)

# Create multiple groups plots
create_multiple_groups_fa_plot(glia_files, "Glia_Sex_Specific")
create_multiple_groups_fa_plot(neuron_files, "Neuron_Sex_Specific")
create_multiple_groups_fa_plot(supernatant_files, "Supernatant_Sex_Specific")

# 2. PAIRWISE A vs B PLOTS - All files in results folder
create_pairwise_fa_plots()

cat("\n=== FA PLOTTING COMPLETE ===\n")
cat("Multiple groups plots: 6 sets (3 cell types × 2 versions each)\n")
cat("Pairwise plots: All _full.csv files × 2 versions each\n")
cat("Plot versions: '_all' (includes zero bars) and '_removed_0' (excludes zero bars)\n")
cat("Target FA types: FA(2:X), FA(14:X), FA(15:X), FA(16:X)\n")
cat("Error bars: Standard Error (fixed for normalized plots)\n")
cat("Normalization: Each FA type max = 100%\n")
cat("Output: new_FA_bar/ directory\n")
cat("Total plots created: ~550+ (4× original due to all/removed_0 × absolute/normalized)\n")





# # ────────────────────────────────────────────────────────────────────────────────
# # New FA Bar Plots - Specific Chain Lengths with Error Bars
# # ────────────────────────────────────────────────────────────────────────────────
# # Creates bar plots for FA(2:X), FA(14:X), FA(15:X), FA(16:X) with all saturation levels
# # Two plot types: multiple groups combined, and pairwise A vs B comparisons
# # Both absolute and normalized versions with standard error bars
# # ────────────────────────────────────────────────────────────────────────────────
# 
# library(tidyverse)
# library(scales)
# 
# # Target chain lengths
# target_chain_lengths <- c(2, 14, 15, 16)
# 
# # Create directory structure
# create_fa_directories <- function() {
#   dir.create("new_FA_bar", showWarnings = FALSE)
#   dir.create("new_FA_bar/multiple_groups", showWarnings = FALSE)
#   dir.create("new_FA_bar/multiple_groups/absolute", showWarnings = FALSE)
#   dir.create("new_FA_bar/multiple_groups/normalized", showWarnings = FALSE)
#   dir.create("new_FA_bar/pairwise", showWarnings = FALSE)
#   dir.create("new_FA_bar/pairwise/absolute", showWarnings = FALSE)
#   dir.create("new_FA_bar/pairwise/normalized", showWarnings = FALSE)
# }
# 
# # Extract FA information from lipid names
# extract_fa_info <- function(lipid_names) {
#   fa_pattern <- "FA\\((\\d+):(\\d+)\\)"
#   matches <- str_match(lipid_names, fa_pattern)
#   
#   data.frame(
#     lipid = lipid_names,
#     chain_length = as.numeric(matches[,2]),
#     saturation = as.numeric(matches[,3]),
#     fa_type = paste0("FA(", matches[,2], ":", matches[,3], ")"),
#     stringsAsFactors = FALSE
#   ) %>%
#     filter(!is.na(chain_length), chain_length %in% target_chain_lengths)
# }
# 
# # Process single file for FA data
# process_fa_file <- function(file_path) {
#   if (!file.exists(file_path)) {
#     warning(paste("File not found:", file_path))
#     return(tibble())
#   }
#   
#   dat <- read_csv(file_path, show_col_types = FALSE) %>%
#     filter(type == "FA")
#   
#   if (nrow(dat) == 0) return(tibble())
#   
#   # Extract group info
#   g1 <- dat$Title_1[1]
#   g2 <- dat$Title_2[1]
#   len1 <- as.numeric(dat$Length1[1])
#   len2 <- as.numeric(dat$Length2[1])
#   
#   # Identify columns
#   all_cols <- colnames(dat)
#   start_idx <- which(all_cols == "type") + 1
#   c1 <- all_cols[start_idx:(start_idx + len1 - 1)]
#   c2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
#   blank_col <- tail(all_cols, 1)
#   
#   # Subtract blank and ensure non-negative
#   dat[c1] <- pmax(dat[c1] - dat[[blank_col]], 0)
#   dat[c2] <- pmax(dat[c2] - dat[[blank_col]], 0)
#   
#   # Convert to numeric
#   dat[c(c1, c2)] <- lapply(dat[c(c1, c2)], function(x) as.numeric(as.character(x)))
#   dat[is.na(dat)] <- 0
#   
#   # Extract FA info
#   fa_info <- extract_fa_info(dat$lipid)
#   if (nrow(fa_info) == 0) return(tibble())
#   
#   # Combine data with FA info
#   dat_fa <- dat %>%
#     select(lipid, all_of(c(c1, c2))) %>%
#     inner_join(fa_info, by = "lipid")
#   
#   # Reshape to long format
#   result <- dat_fa %>%
#     pivot_longer(cols = all_of(c(c1, c2)), names_to = "sample", values_to = "value") %>%
#     mutate(
#       group = ifelse(sample %in% c1, g1, g2),
#       file_source = basename(file_path)
#     ) %>%
#     select(fa_type, chain_length, saturation, group, sample, value, file_source)
#   
#   return(result)
# }
# 
# # Create multiple groups plot
# create_multiple_groups_fa_plot <- function(files, save_name) {
#   
#   # Process all files
#   all_data <- map_dfr(files, ~process_fa_file(file.path("results", 
#                                                         ifelse(endsWith(.x, ".csv"), .x, paste0(.x, ".csv")))))
#   
#   if (nrow(all_data) == 0) {
#     warning("No FA data found in files")
#     return(invisible())
#   }
#   
#   # Calculate summary statistics
#   summary_data <- all_data %>%
#     group_by(fa_type, chain_length, saturation, group) %>%
#     summarise(
#       mean_value = mean(value, na.rm = TRUE),
#       se_value = sd(value, na.rm = TRUE) / sqrt(n()),
#       n_samples = n(),
#       .groups = "drop"
#     ) %>%
#     arrange(chain_length, saturation)
#   
#   # Create factor with ordered levels
#   summary_data$fa_type <- factor(summary_data$fa_type, 
#                                  levels = unique(summary_data$fa_type[order(summary_data$chain_length, summary_data$saturation)]))
#   
#   # ABSOLUTE PLOT
#   p_abs <- ggplot(summary_data, aes(x = fa_type, y = mean_value, fill = group)) +
#     geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
#     geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
#                   position = position_dodge(width = 0.8), width = 0.4) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = "bottom") +
#     labs(title = paste("FA Chain Lengths - Multiple Groups (Absolute) -", save_name),
#          x = "Fatty Acid",
#          y = "Absolute Value",
#          fill = "Group") +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
#   
#   # NORMALIZED PLOT (each FA type normalized to its max across all groups)
#   norm_data <- summary_data %>%
#     group_by(fa_type) %>%
#     mutate(
#       max_value = max(mean_value, na.rm = TRUE),
#       norm_mean = (mean_value / max_value) * 100,
#       norm_se = (se_value / max_value) * 100
#     ) %>%
#     ungroup()
#   
#   p_norm <- ggplot(norm_data, aes(x = fa_type, y = norm_mean, fill = group)) +
#     geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
#     geom_errorbar(aes(ymin = norm_mean - norm_se, ymax = norm_mean + norm_se),
#                   position = position_dodge(width = 0.8), width = 0.4) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = "bottom") +
#     labs(title = paste("FA Chain Lengths - Multiple Groups (Normalized) -", save_name),
#          x = "Fatty Acid",
#          y = "Relative % (Max = 100%)",
#          fill = "Group") +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, 105))
#   
#   # Save plots
#   ggsave(file.path("new_FA_bar/multiple_groups/absolute", paste0(save_name, "_absolute.pdf")), 
#          p_abs, width = 14, height = 8)
#   ggsave(file.path("new_FA_bar/multiple_groups/absolute", paste0(save_name, "_absolute.png")), 
#          p_abs, width = 14, height = 8, dpi = 300)
#   
#   ggsave(file.path("new_FA_bar/multiple_groups/normalized", paste0(save_name, "_normalized.pdf")), 
#          p_norm, width = 14, height = 8)
#   ggsave(file.path("new_FA_bar/multiple_groups/normalized", paste0(save_name, "_normalized.png")), 
#          p_norm, width = 14, height = 8, dpi = 300)
#   
#   cat("Multiple groups FA plots saved for:", save_name, "\n")
#   return(invisible(list(absolute = p_abs, normalized = p_norm, data = summary_data)))
# }
# 
# # Create pairwise A vs B plots
# create_pairwise_fa_plots <- function() {
#   
#   # Get all full files
#   all_files <- list.files("results", pattern = "_full\\.csv$", full.names = FALSE)
#   
#   for (file in all_files) {
#     file_base <- gsub("\\.csv$", "", file)
#     
#     # Process single file
#     data <- process_fa_file(file.path("results", file))
#     
#     if (nrow(data) == 0) {
#       cat("No FA data in:", file, "\n")
#       next
#     }
#     
#     # Calculate summary statistics
#     summary_data <- data %>%
#       group_by(fa_type, chain_length, saturation, group) %>%
#       summarise(
#         mean_value = mean(value, na.rm = TRUE),
#         se_value = sd(value, na.rm = TRUE) / sqrt(n()),
#         n_samples = n(),
#         .groups = "drop"
#       ) %>%
#       arrange(chain_length, saturation)
#     
#     # Create factor with ordered levels
#     summary_data$fa_type <- factor(summary_data$fa_type, 
#                                    levels = unique(summary_data$fa_type[order(summary_data$chain_length, summary_data$saturation)]))
#     
#     # ABSOLUTE PLOT
#     p_abs <- ggplot(summary_data, aes(x = fa_type, y = mean_value, fill = group)) +
#       geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
#       geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
#                     position = position_dodge(width = 0.8), width = 0.4) +
#       theme_classic() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1),
#             legend.position = "bottom") +
#       labs(title = paste("FA Chain Lengths - Pairwise (Absolute)"),
#            subtitle = file_base,
#            x = "Fatty Acid",
#            y = "Absolute Value",
#            fill = "Group") +
#       scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
#     
#     # NORMALIZED PLOT
#     norm_data <- summary_data %>%
#       group_by(fa_type) %>%
#       mutate(
#         max_value = max(mean_value, na.rm = TRUE),
#         norm_mean = (mean_value / max_value) * 100,
#         norm_se = (se_value / max_value) * 100
#       ) %>%
#       ungroup()
#     
#     p_norm <- ggplot(norm_data, aes(x = fa_type, y = norm_mean, fill = group)) +
#       geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
#       geom_errorbar(aes(ymin = norm_mean - norm_se, ymax = norm_mean + norm_se),
#                     position = position_dodge(width = 0.8), width = 0.4) +
#       theme_classic() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1),
#             legend.position = "bottom") +
#       labs(title = paste("FA Chain Lengths - Pairwise (Normalized)"),
#            subtitle = file_base,
#            x = "Fatty Acid",
#            y = "Relative % (Max = 100%)",
#            fill = "Group") +
#       scale_y_continuous(expand = c(0, 0), limits = c(0, 105))
#     
#     # Save plots
#     safe_name <- gsub("[^A-Za-z0-9_-]", "_", file_base)
#     
#     ggsave(file.path("new_FA_bar/pairwise/absolute", paste0(safe_name, "_absolute.pdf")), 
#            p_abs, width = 12, height = 8)
#     ggsave(file.path("new_FA_bar/pairwise/absolute", paste0(safe_name, "_absolute.png")), 
#            p_abs, width = 12, height = 8, dpi = 300)
#     
#     ggsave(file.path("new_FA_bar/pairwise/normalized", paste0(safe_name, "_normalized.pdf")), 
#            p_norm, width = 12, height = 8)
#     ggsave(file.path("new_FA_bar/pairwise/normalized", paste0(safe_name, "_normalized.png")), 
#            p_norm, width = 12, height = 8, dpi = 300)
#     
#     cat("Pairwise FA plots saved for:", safe_name, "\n")
#   }
# }
# 
# # ────────────────────────────────────────────────────────────────────────────────
# # MAIN EXECUTION
# # ────────────────────────────────────────────────────────────────────────────────
# 
# # Create directories
# create_fa_directories()
# 
# # 1. MULTIPLE GROUPS PLOTS - Sex-specific files
# glia_files <- c(
#   "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ Female vs Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ FemaleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ Female vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ FemaleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ maleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ maleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ maleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ Female vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ FemaleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ maleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ maleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
#   "Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ male vs Sample Type_ Glia__Sleep Condition_ Waking up__Sex_ maleBlank1_full"
# )
# 
# neuron_files <- c(
#   "Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ Female vs Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ FemaleBlank1_full",
#   "Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ Female vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
#   "Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ FemaleBlank1_full",
#   "Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ maleBlank1_full",
#   "Sample Type_ Neuron__Sleep Condition_ Deprived__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ maleBlank1_full",
#   "Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ Female vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
#   "Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ FemaleBlank1_full",
#   "Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ maleBlank1_full",
#   "Sample Type_ Neuron__Sleep Condition_ Going to sleep__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ maleBlank1_full",
#   "Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ FemaleBlank1_full",
#   "Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ male vs Sample Type_ Neuron__Sleep Condition_ Waking up__Sex_ maleBlank1_full"
# )
# 
# supernatant_files <- c(
#   "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_FemalevsSampleType_Supernatant60BrainSleepCondition_GoingtosleepSex_FemaleBlank1_full",
#   "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_FemalevsSampleType_Supernatant60BrainSleepCondition_WakingupSex_FemaleBlank1_full",
#   "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_DeprivedSex_FemaleBlank1_full",
#   "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_DeprivedSex_maleBlank1_full",
#   "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_GoingtosleepSex_maleBlank1_full",
#   "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_WakingupSex_maleBlank1_full"
# )
# 
# # Create multiple groups plots
# create_multiple_groups_fa_plot(glia_files, "Glia_Sex_Specific")
# create_multiple_groups_fa_plot(neuron_files, "Neuron_Sex_Specific")
# create_multiple_groups_fa_plot(supernatant_files, "Supernatant_Sex_Specific")
# 
# # 2. PAIRWISE A vs B PLOTS - All files in results folder
# create_pairwise_fa_plots()
# 
# cat("\n=== FA PLOTTING COMPLETE ===\n")
# cat("Multiple groups plots: 3 sets (Glia, Neuron, Supernatant)\n")
# cat("Pairwise plots: All _full.csv files in results/\n")
# cat("Target FA types: FA(2:X), FA(14:X), FA(15:X), FA(16:X)\n")
# cat("Error bars: Standard Error\n")
# cat("Normalization: Each FA type max = 100%\n")
# cat("Output: new_FA_bar/ directory\n")