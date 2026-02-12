#Goal1 : comparing the metabolite proflie of GBM patients in primary and recurrent stages(Paired analysis)
#Note : the metabolites just have the HMDB code have been included in the analysis (line 76).
#Gaol2:  Patient-Specific Metabolite Profiling in GBM (just based on logFC)
# Load necessary libraries
# ---------------------------
# required libraries
# ---------------------------
library(tidyverse)
library(limma)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggrepel)
library(UpSetR)
library(pheatmap)
library(grid)
library(gridExtra)
library(EnhancedVolcano)



# ---------------------------
# Load and preprocess metabolite data
# ---------------------------
metabolites_total <- read.csv("metabolites_total.csv", stringsAsFactors = FALSE)

# Extract unique metabolites per case and find common metabolites
sets <- lapply(metabolites_total, function(x) unique(na.omit(x)))
common_metabolites <- Reduce(intersect, sets)


common_metabs <- read.csv("data1/common_metabolites.csv", stringsAsFactors = FALSE)[[1]]

case_files <- list(
  Case1 = "data1/case1.csv",
  Case2 = "data1/case2.csv",
  Case3 = "data1/case3.csv",
  Case4 = "data1/case4.csv"
)

# Filter each case to keep only common metabolites
filtered_cases <- lapply(case_files, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  df_filtered <- df[df$Metabolite %in% common_metabs, ]
  return(df_filtered)
})



# Convert to long format
long_data_list <- lapply(names(filtered_cases), function(case_name) {
  df <- filtered_cases[[case_name]]
  df_long <- df %>%
    pivot_longer(cols = c("primary", "recurrent"),
                 names_to = "Stage",
                 values_to = "Value") %>%
    mutate(Case = case_name)
  return(df_long)
})

combined_data <- bind_rows(long_data_list)

# Filter metabolites with data in all 4 cases
complete_metabolites <- combined_data %>%
  group_by(Metabolite, Case) %>%
  summarize(non_na = sum(!is.na(Value)), .groups = "drop") %>%
  group_by(Metabolite) %>%
  summarize(num_cases_with_data = sum(non_na >= 2)) %>%
  filter(num_cases_with_data == 4) %>%
  pull(Metabolite)

filtered_combined <- combined_data %>%
  filter(Metabolite %in% complete_metabolites)

# Remove metabolites without HMDB.ID and log2-transform
filtered_combined_clean <- filtered_combined %>%
  filter(!is.na(HMDB.ID) & HMDB.ID != "") %>%
  mutate(Value = log2(Value + 1)) %>%
  mutate(Metabolite = str_trim(Metabolite) %>% tolower())

# ---------------------------
# LIMMA analysis (paired)
# ---------------------------
data <- filtered_combined_clean
#colnames(data)[colnames(data) == "Stage"] <- "Stage"
data$Case <- as.factor(data$Case)
data$Stage <- factor(data$Stage, levels = c("primary", "recurrent"))

# Wide format expression matrix
expr_wide <- data %>%
  mutate(Sample = paste0(Case, "_", Stage)) %>%
  select(Metabolite, Sample, Value) %>%
  pivot_wider(names_from = Sample, values_from = Value)

metabolite_names <- expr_wide$Metabolite
expr_mat <- expr_wide %>% select(-Metabolite) %>% as.matrix()
rownames(expr_mat) <- metabolite_names

# Sample metadata
sample_names <- colnames(expr_mat)
targets <- data.frame(
  Sample = sample_names,
  Case = gsub("_.*", "", sample_names),
  Stage = gsub(".*_", "", sample_names)
)
targets$Case <- factor(targets$Case)
targets$Stage <- factor(targets$Stage, levels = c("primary", "recurrent"))

# Design matrix and fit
design <- model.matrix(~Case + Stage, data = targets)
fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)

# Extract results for recurrent vs primary
results <- topTable(fit, coef = "Stagerecurrent", number = Inf, adjust.method = "fdr")
results <- results %>%
  tibble::rownames_to_column("Metabolite") %>%
  select(Metabolite, AveExpr, logFC, P.Value, adj.P.Val)

# Join HMDB IDs
meta_hmdb <- filtered_combined_clean %>%
  distinct(Metabolite, HMDB.ID)
results <- results %>% left_join(meta_hmdb, by = "Metabolite")
row.names(results) <- results$Metabolite

# ---------------------------
# Selecting metabolites 
# ---------------------------
filtered_results <- results %>%
  filter(abs(logFC) > 0.2 & P.Value < 0.05)

# Custom labels for Volcano plot
results$custom_label <- ifelse(
  results$Metabolite %in% filtered_results$Metabolite,
  paste0(results$Metabolite, 
         "\nlog2FC=", round(results$logFC, 2), 
         "\np=", signif(results$P.Value, 3)),
  ""
)


################################################################################
#-------------------------------------------------------------------------
#Subject Specific Analysis
#--------------------------------------------------------------------------
# -----------------------------
# Computing patient-specific log2FC
# -----------------------------

sample_cols <- colnames(expr_wide)[-1]

sample_info <- data.frame(
  Sample = sample_cols,
  Case = str_extract(sample_cols, "^[^_]+"),
  stage = str_extract(sample_cols, "[^_]+$")
)

case_list <- unique(sample_info$Case)

logFC_df <- data.frame(Metabolite = expr_wide$Metabolite)

for (case in case_list) {
  primary_col <- sample_info$Sample[sample_info$Case == case & sample_info$stage == "primary"]
  recurrent_col  <- sample_info$Sample[sample_info$Case == case & sample_info$stage == "recurrent"]
  
  if (length(primary_col) == 1 && length(recurrent_col) == 1) {
    logFC_df[[paste0("log2FC_", case)]] <- expr_wide[[recurrent_col]] - expr_wide[[primary_col]]
  }
}


# -----------------------------
# Creating logFC heatmap 
#(This is not used in the manuscript)
# -----------------------------

rownames(logFC_df) <- logFC_df$Metabolite
logFC_matrix <- logFC_df[, -1]
logFC_matrix <- as.matrix(logFC_matrix)

logFC_matrix_filtered <- logFC_matrix[apply(abs(logFC_matrix), 1, function(x) any(x > 1)), ]

# Creating high-diff metabolite list 


lfc_thresh <- 1

high_diff_list <- lapply(colnames(logFC_matrix), function(case) {
  vec <- logFC_matrix[, case]
  keep_idx <- which(abs(vec) > lfc_thresh)
  
  if (length(keep_idx) == 0) {
    return(data.frame(Case = character(0), Metabolite = character(0), logFC = numeric(0)))
  }
  
  data.frame(
    Case = case,
    Metabolite = names(vec)[keep_idx],
    logFC = vec[keep_idx],
    row.names = NULL
  )
})

high_diff_df <- do.call(rbind, high_diff_list)


top_metab_list <- split(high_diff_df$Metabolite, high_diff_df$Case)


#-------------------------------------------------------------------------------
# Assume we have case1_top, case2_top, case3_top, case4_top

# Defining top metabolite sets for each case
case1_top <- unique(high_diff_df$Metabolite[high_diff_df$Case == colnames(logFC_matrix)[1]])
case2_top <- unique(high_diff_df$Metabolite[high_diff_df$Case == colnames(logFC_matrix)[2]])
case3_top <- unique(high_diff_df$Metabolite[high_diff_df$Case == colnames(logFC_matrix)[3]])
case4_top <- unique(high_diff_df$Metabolite[high_diff_df$Case == colnames(logFC_matrix)[4]])

# Defining a named list
case_list <- list(
  Case1 = case1_top,
  Case2 = case2_top,
  Case3 = case3_top,
  Case4 = case4_top
)

# Helper to get intersections
get_unique_intersection <- function(keep, drop) {
  intersected <- Reduce(intersect, case_list[keep])
  dropped <- Reduce(union, case_list[drop])
  setdiff(intersected, dropped)
}

# Creating named intersections
intersections <- list(
  "Only Case1" = get_unique_intersection("Case1", c("Case2", "Case3", "Case4")),
  "Only Case2" = get_unique_intersection("Case2", c("Case1", "Case3", "Case4")),
  "Only Case3" = get_unique_intersection("Case3", c("Case1", "Case2", "Case4")),
  "Only Case4" = get_unique_intersection("Case4", c("Case1", "Case2", "Case3")),
  "Case1 & Case2 only" = get_unique_intersection(c("Case1", "Case2"), c("Case3", "Case4")),
  "Case1 & Case3 only" = get_unique_intersection(c("Case1", "Case3"), c("Case2", "Case4")),
  "Case1 & Case4 only" = get_unique_intersection(c("Case1", "Case4"), c("Case2", "Case3")),
  "Case2 & Case3 only" = get_unique_intersection(c("Case2", "Case3"), c("Case1", "Case4")),
  "Case2 & Case4 only" = get_unique_intersection(c("Case2", "Case4"), c("Case1", "Case3")),
  "Case3 & Case4 only" = get_unique_intersection(c("Case3", "Case4"), c("Case1", "Case2")),
  "Case1 & Case2 & Case3 only" = get_unique_intersection(c("Case1", "Case2", "Case3"), "Case4"),
  "Case1 & Case2 & Case4 only" = get_unique_intersection(c("Case1", "Case2", "Case4"), "Case3"),
  "Case1 & Case3 & Case4 only" = get_unique_intersection(c("Case1", "Case3", "Case4"), "Case2"),
  "Case2 & Case3 & Case4 only" = get_unique_intersection(c("Case2", "Case3", "Case4"), "Case1"),
  "All 4 Cases" = Reduce(intersect, case_list)
)

# Converting to data.frame
intersect_df <- do.call(rbind, lapply(names(intersections), function(name) {
  if (length(intersections[[name]]) == 0) return(NULL)
  data.frame(Comparison = name, Metabolite = intersections[[name]])
}))

#extracting HMDB.ID for intersections
HMDB_Only_Case1 <- unique(filtered_combined_clean$HMDB.ID[filtered_combined_clean$Metabolite %in% intersections$`Only Case1`])
HMDB_Only_Case2 <- unique(filtered_combined_clean$HMDB.ID[filtered_combined_clean$Metabolite %in% intersections$`Only Case2`])
HMDB_Only_Case3 <- unique(filtered_combined_clean$HMDB.ID[filtered_combined_clean$Metabolite %in% intersections$`Only Case3`])
HMDB_Only_Case4 <- unique(filtered_combined_clean$HMDB.ID[filtered_combined_clean$Metabolite %in% intersections$`Only Case4`])

# Ensuring each list is named properly
names(top_metab_list) <- gsub("logFC_", "", names(top_metab_list))

#-------------------------------------------------------------------------------
# Building a lookup table of metabolite to HMDB.ID
meta_hmdb_map <- filtered_combined_clean %>%
  distinct(Metabolite, HMDB.ID)

# Building a unified data frame for all intersection categories
intersect_with_hmdb <- do.call(rbind, lapply(names(intersections), function(name) {
  metab_list <- intersections[[name]]
  if (length(metab_list) == 0) return(NULL)
  
  # Creating data frame for each group
  df <- data.frame(
    Comparison = name,
    Metabolite = metab_list,
    stringsAsFactors = FALSE
  )
  
  # Adding HMDB.ID using left_join
  df <- left_join(df, meta_hmdb_map, by = "Metabolite")
  return(df)
}))


head(intersect_with_hmdb)


