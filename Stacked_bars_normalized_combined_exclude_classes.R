# ────────────────────────────────────────────────────────────────────────────────
# Hardcoded sex-specific file analysis for Glia, Neuron, and Supernatant groups
# ────────────────────────────────────────────────────────────────────────────────

# GLIA - Sex-specific files only (hardcoded, duplicates removed)
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
  # Removed all SampleType_Glia* files as they are duplicates of the above
)

lipid_combined_bar(
  files = glia_files,
  save_name = "Glia_Sex_Specific_All_Groups",
  exclude_classes = c("FA")
)

# NEURON - Sex-specific files only (hardcoded, duplicates removed)
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
  # Removed all SampleType_Neuron* files as they are duplicates of the above
)

lipid_combined_bar(
  files = neuron_files,
  save_name = "Neuron_Sex_Specific_All_Groups", 
  exclude_classes = c("FA")
)

# SUPERNATANT - Sex-specific files only (hardcoded)
supernatant_files <- c(
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_FemalevsSampleType_Supernatant60BrainSleepCondition_GoingtosleepSex_FemaleBlank1_full",
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_FemalevsSampleType_Supernatant60BrainSleepCondition_WakingupSex_FemaleBlank1_full",
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_DeprivedSex_FemaleBlank1_full",
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_DeprivedSex_maleBlank1_full",
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_GoingtosleepSex_maleBlank1_full",
  "SampleType_Supernatant60BrainSleepCondition_DeprivedSex_malevsSampleType_Supernatant60BrainSleepCondition_WakingupSex_maleBlank1_full"
)

lipid_combined_bar(
  files = supernatant_files,
  save_name = "Supernatant_Sex_Specific_All_Groups",
  exclude_classes = c("FA")
)

# Print summary of what was processed
cat("=== PROCESSING SUMMARY ===\n")
cat("Glia files processed:", length(glia_files), "(duplicates removed)\n")
cat("Neuron files processed:", length(neuron_files), "(duplicates removed)\n") 
cat("Supernatant files processed:", length(supernatant_files), "\n")
cat("Excluded lipid class: FA\n")
cat("Output saved to new_plots/ directory\n")
cat("\nNote: Removed duplicate comparisons with different filename formats\n")# ────────────────────────────────────────────────────────────────────────────────
# Combined Lipid Class Plots (FA Format Style)
# ────────────────────────────────────────────────────────────────────────────────
# Arguments
#   files        : vector of basenames of CSV result files *without* the ".csv"
#                  (they are read from the folder results/ and ".csv" is added)
#   save_name    : prefix for all output files
#   exclude_classes : vector of lipid classes to exclude (optional)
#
# Output
#   • new_plots/normalized/     <save_name>_normalized.{pdf,png}
#   • new_plots/absolute/       <save_name>_absolute.{pdf,png}
#
# Requirements:
#   - tidyverse loaded
# ────────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(scales)

# Define lipid classes and colors
lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG", "DAG", "TAG | DAG", "DAG | CE", "TAG | DAG | CE")
lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072")

lipid_combined_bar <- function(files, save_name, exclude_classes = NULL) {
  
  # Create output directories
  dir.create("new_plots", showWarnings = FALSE)
  dir.create("new_plots/normalized", showWarnings = FALSE)
  dir.create("new_plots/absolute", showWarnings = FALSE)
  
  # Filter out excluded classes
  if (!is.null(exclude_classes)) {
    valid_classes <- lipid_classes[!lipid_classes %in% exclude_classes]
    valid_colors <- lipid_colors[!lipid_classes %in% exclude_classes]
  } else {
    valid_classes <- lipid_classes
    valid_colors <- lipid_colors
  }
  
  # Create named color vector
  lipid_class_colors <- setNames(alpha(valid_colors, 0.7), valid_classes)
  
  # Helper function to summarise one CSV file
  summarise_file <- function(base) {
    csv <- file.path("results",
                     ifelse(endsWith(base, ".csv"), base, paste0(base, ".csv")))
    
    if (!file.exists(csv)) {
      warning(paste("File not found:", csv))
      return(tibble())
    }
    
    dat <- read_csv(csv, show_col_types = FALSE)
    if (nrow(dat) == 0) return(tibble())
    
    # Get group names and lengths
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
    
    # Subtract blank and ensure non-negative values
    dat[c1] <- pmax(dat[c1] - dat[[blank_col]], 0)
    dat[c2] <- pmax(dat[c2] - dat[[blank_col]], 0)
    
    # Convert to numeric and handle NAs
    dat[c(c1, c2)] <- lapply(dat[c(c1, c2)], function(x) as.numeric(as.character(x)))
    dat[is.na(dat)] <- 0
    
    # Calculate means
    dat <- dat %>% mutate(
      mean1 = rowMeans(across(all_of(c1))),
      mean2 = rowMeans(across(all_of(c2)))
    )
    
    # Filter for valid lipid classes
    if (!is.null(exclude_classes)) {
      dat <- dat %>% filter(!type %in% exclude_classes)
    }
    
    # Aggregate by lipid type
    agg <- dat %>%
      group_by(type) %>%
      summarise(sum_mean1 = sum(mean1),
                sum_mean2 = sum(mean2), .groups = "drop")
    
    # Calculate totals for percentage
    tot1 <- sum(agg$sum_mean1)
    tot2 <- sum(agg$sum_mean2)
    
    # Return combined data
    bind_rows(
      agg %>% transmute(group = g1,
                        lipid_class = factor(type, levels = valid_classes),
                        value = sum_mean1,
                        percent = 100 * sum_mean1 / tot1),
      agg %>% transmute(group = g2,
                        lipid_class = factor(type, levels = valid_classes),
                        value = sum_mean2,
                        percent = 100 * sum_mean2 / tot2)
    )
  }
  
  # Process all files and combine
  comb_df <- map_dfr(files, summarise_file) %>%
    distinct() %>%  # Remove perfect duplicates
    filter(!is.na(lipid_class))  # Remove any NA lipid classes
  
  if (nrow(comb_df) == 0) {
    warning("No data found in any files")
    return(invisible())
  }
  
  # Helper function to calculate label positions
  centre_labels <- function(df, val_col) {
    df %>%
      group_by(group) %>% 
      arrange(desc(lipid_class)) %>%
      mutate(position = cumsum(.data[[val_col]]) - 0.5 * .data[[val_col]]) %>%
      ungroup()
  }
  
  # ----------------------------- ABSOLUTE plot --------------------------------
  abs_df <- comb_df
  abs_lab <- centre_labels(abs_df, "value")
  
  p_abs <- ggplot(abs_df, aes(group, value, fill = lipid_class)) +
    geom_col() +
    geom_text(data = abs_lab,
              aes(y = position, label = sprintf("%.1f%%", percent)),
              size = 3, colour = "black") +
    scale_fill_manual(values = lipid_class_colors,
                      name = "Lipid Class") +
    scale_y_continuous(expand = c(0, 0)) +
    ylab("Total lipid content") +
    xlab("Groups") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Combined Lipid Class Distribution (absolute)")
  
  # Save absolute plots
  ggsave(file.path("new_plots", "absolute", paste0(save_name, "_absolute.pdf")), 
         p_abs, width = 12, height = 8)
  ggsave(file.path("new_plots", "absolute", paste0(save_name, "_absolute.png")), 
         p_abs, width = 12, height = 8, dpi = 300)
  
  # --------------------------- normalized plot --------------------------------
  norm_df <- comb_df %>% select(group, lipid_class, percent)
  norm_lab <- centre_labels(norm_df, "percent")
  
  p_norm <- ggplot(norm_df, aes(group, percent, fill = lipid_class)) +
    geom_col() +
    geom_text(data = norm_lab,
              aes(y = position, label = sprintf("%.1f%%", percent)),
              size = 3, colour = "black") +
    scale_fill_manual(values = lipid_class_colors,
                      name = "Lipid Class") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 100),
                       expand = c(0, 0)) +
    ylab("Proportion of total lipid content (%)") +
    xlab("Groups") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Combined Lipid Class Distribution (normalized)")
  
  # Save normalized plots
  ggsave(file.path("new_plots", "normalized", paste0(save_name, "_normalized.pdf")), 
         p_norm, width = 12, height = 8)
  ggsave(file.path("new_plots", "normalized", paste0(save_name, "_normalized.png")), 
         p_norm, width = 12, height = 8, dpi = 300)
  
  # Print summary
  cat("Generated plots for", length(unique(comb_df$group)), "groups\n")
  cat("Included lipid classes:", paste(levels(comb_df$lipid_class), collapse = ", "), "\n")
  if (!is.null(exclude_classes)) {
    cat("Excluded lipid classes:", paste(exclude_classes, collapse = ", "), "\n")
  }
  
  return(invisible(list(absolute_plot = p_abs, normalized_plot = p_norm, data = comb_df)))
}

# ────────────────────────────────────────────────────────────────────────────────
# Example usage:
# ────────────────────────────────────────────────────────────────────────────────


# ────────────────────────────────────────────────────────────────────────────────
# Sex-specific analysis for Glia, Neuron, and Supernatant groups
# ────────────────────────────────────────────────────────────────────────────────

# Get all full files and remove .csv extension
all_full_files <- list.files("results", pattern = "_full\\.csv$", full.names = FALSE)
all_full_files <- gsub("\\.csv$", "", all_full_files)

# Filter for sex-specific files only (containing "male" or "Female")
sex_specific_files <- all_full_files[grepl("male|Female", all_full_files, ignore.case = TRUE)]

# 1. GLIA - Sex-specific groups only
glia_files <- sex_specific_files[grepl("Glia", sex_specific_files, ignore.case = TRUE)]
lipid_combined_bar(
  files = glia_files,
  save_name = "Glia_Sex_Specific_All_Groups",
  exclude_classes = c("FA")
)

# 2. NEURON - Sex-specific groups only  
neuron_files <- sex_specific_files[grepl("Neuron", sex_specific_files, ignore.case = TRUE)]
lipid_combined_bar(
  files = neuron_files,
  save_name = "Neuron_Sex_Specific_All_Groups", 
  exclude_classes = c("FA")
)

# 3. SUPERNATANT - Sex-specific groups only
supernatant_files <- sex_specific_files[grepl("Supernatant", sex_specific_files, ignore.case = TRUE)]
lipid_combined_bar(
  files = supernatant_files,
  save_name = "Supernatant_Sex_Specific_All_Groups",
  exclude_classes = c("FA")
)

# Print summary of what was processed
cat("=== PROCESSING SUMMARY ===\n")
cat("Total sex-specific files found:", length(sex_specific_files), "\n")
cat("Glia files processed:", length(glia_files), "\n")
cat("Neuron files processed:", length(neuron_files), "\n") 
cat("Supernatant files processed:", length(supernatant_files), "\n")
cat("Excluded lipid class: FA\n")
