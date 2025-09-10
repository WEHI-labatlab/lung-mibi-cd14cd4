library(ggpattern)
library(ggpmisc)
library(dplyr)
library(tidyr)

source("src/proportion/calculate_cell_type_percentages.R")
source("src/plot_survival_curve.R")

# read the patient counts
patient_counts <- readRDS(file="data/NSCLCcohort/processed/patient_counts_compartment_3_all_FM_TRM_macro_CD14_df.rds")

run_name <- "B_LUAD_LLM_RFSsurvival_curves_3_compartments"

patient_counts <- patient_counts %>%
    filter(`Cancer Type` == "LUAD" & `MIBI clusters` != "M")

# group col name
group_cols <- c("Patient", "RFS (days)", "Relapse", "Compartment (3)")

# cell types of interest
cell_types <- c(
    "B cells"
)

# the parent subtype
parent_type <- c("Immune cells")

### Calculate the proportions on demand
output_folder <- "data/NSCLCcohort/output/Figure3/survival_curves"
if (!dir.exists(output_folder)){
    dir.create(output_folder, recursive = T)
}

# calculate the proportions
percentages_df <- calculate_cell_type_percentages(patient_counts, group_cols, cell_types, parent_type)

# we need to filter out the proportions with less than 5 cells
percentages_df <- percentages_df %>% 
    filter(parent_type_count >= 5)


# ------------------------------------------------------------------------------
# Median
# ------------------------------------------------------------------------------

# Boundary 
threshold <- median(percentages_df %>% filter(`Compartment (3)`=="Boundary") %>% pull(value))
percentages_df$B_high <- percentages_df$value > threshold

percentages_df$Relapse[!(percentages_df$Relapse %in% c("TRUE", "FALSE"))] <- NA
percentages_df$Relapse<- as.logical(percentages_df$Relapse)

plot_survival_curve(percentages_df %>% filter(!is.na(!!as.name("Relapse"))) %>% filter(`Compartment (3)`=="Boundary"), 
                    "RFS (days)", 
                    "Relapse", 
                    "B_high",
                    title = "RFS LLM B high low Boundary median",
                    calc.p.value = TRUE,
                    filename = file.path(output_folder, paste0(run_name,"_median_Boundary.png")))


# Stroma Core 
threshold <- median(percentages_df %>% filter(`Compartment (3)`=="Stroma Core") %>% pull(value))
percentages_df$B_high <- percentages_df$value > threshold

percentages_df$Relapse[!(percentages_df$Relapse %in% c("TRUE", "FALSE"))] <- NA
percentages_df$Relapse<- as.logical(percentages_df$Relapse)

plot_survival_curve(percentages_df %>% filter(!is.na(!!as.name("Relapse"))) %>% filter(`Compartment (3)`=="Stroma Core"), 
                    "RFS (days)", 
                    "Relapse", 
                    "B_high",
                    title = "RFS LLM B high low Stroma Core median",
                    calc.p.value = TRUE,
                    filename = file.path(output_folder, paste0(run_name, "_median_Stroma_Core.png")))


# Tumour Core 
threshold <- median(percentages_df %>% filter(`Compartment (3)`=="Tumour Core") %>% pull(value))
percentages_df$B_high <- percentages_df$value > threshold

percentages_df$Relapse[!(percentages_df$Relapse %in% c("TRUE", "FALSE"))] <- NA
percentages_df$Relapse<- as.logical(percentages_df$Relapse)

plot_survival_curve(percentages_df %>% filter(!is.na(!!as.name("Relapse"))) %>% filter(`Compartment (3)`=="Tumour Core"), 
                    "RFS (days)", 
                    "Relapse", 
                    "B_high",
                    title = "RFS LLM B high low Tumour Core median",
                    calc.p.value = TRUE,
                    filename = file.path(output_folder, paste0(run_name, "_median_Tumour_Core.png")))

# ------------------------------------------------------------------------------
# Lower and upper Quartile
# ------------------------------------------------------------------------------

# Boundary 

quantiles <- quantile(percentages_df %>% filter(`Compartment (3)`=="Boundary") %>% pull(value), probs = seq(0, 1, 0.25))
percentages_df$B_high <- "Medium"
percentages_df$B_high[percentages_df$value <= quantiles[2]] <- "Low"
percentages_df$B_high[percentages_df$value >= quantiles[4]] <- "High"

percentages_df$Relapse[!(percentages_df$Relapse %in% c("TRUE", "FALSE"))] <- NA
percentages_df$Relapse<- as.logical(percentages_df$Relapse)

plot_survival_curve(percentages_df %>% filter(!is.na(!!as.name("Relapse"))) %>% filter(`Compartment (3)`=="Boundary") %>% filter(B_high!="Medium"), 
                    "RFS (days)", 
                    "Relapse", 
                    "B_high",
                    title = "RFS LLM B high low Boundary quartile",
                    calc.p.value = TRUE,
                    filename = file.path(output_folder, paste0(run_name,"_quartile_Boundary.png")))


# Stroma Core 
quantiles <- quantile(percentages_df %>% filter(`Compartment (3)`=="Stroma Core") %>% pull(value), probs = seq(0, 1, 0.25))
percentages_df$B_high <- "Medium"
percentages_df$B_high[percentages_df$value <= quantiles[2]] <- "Low"
percentages_df$B_high[percentages_df$value >= quantiles[4]] <- "High"

percentages_df$Relapse[!(percentages_df$Relapse %in% c("TRUE", "FALSE"))] <- NA
percentages_df$Relapse<- as.logical(percentages_df$Relapse)

plot_survival_curve(percentages_df %>% filter(!is.na(!!as.name("Relapse"))) %>% filter(`Compartment (3)`=="Stroma Core") %>% filter(B_high!="Medium"), 
                    "RFS (days)", 
                    "Relapse", 
                    "B_high",
                    title = "RFS LLM B high low Stroma Core quartile",
                    calc.p.value = TRUE,
                    filename = file.path(output_folder, paste0(run_name, "_quartile_Stroma_Core.png")))


# Tumour Core 
quantiles <- quantile(percentages_df %>% filter(`Compartment (3)`=="Tumour Core") %>% pull(value), probs = seq(0, 1, 0.25))
percentages_df$B_high <- "Medium"
percentages_df$B_high[percentages_df$value <= quantiles[2]] <- "Low"
percentages_df$B_high[percentages_df$value >= quantiles[4]] <- "High"

percentages_df$Relapse[!(percentages_df$Relapse %in% c("TRUE", "FALSE"))] <- NA
percentages_df$Relapse<- as.logical(percentages_df$Relapse)

plot_survival_curve(percentages_df %>% filter(!is.na(!!as.name("Relapse"))) %>% filter(`Compartment (3)`=="Tumour Core") %>% filter(B_high!="Medium"), 
                    "RFS (days)", 
                    "Relapse", 
                    "B_high",
                    title = "RFS LLM B high low Tumour Core quartile",
                    calc.p.value = TRUE,
                    filename = file.path(output_folder, paste0(run_name, "_quartile_Tumour_Core.png")))

# ------------------------------------------------------------------------------
# Lower Quartile
# ------------------------------------------------------------------------------

# Boundary 

quantiles <- quantile(percentages_df %>% filter(`Compartment (3)`=="Boundary") %>% pull(value), probs = seq(0, 1, 0.25))
percentages_df$B_high <- "MediumHigh"
percentages_df$B_high[percentages_df$value <= quantiles[2]] <- "Low"


percentages_df$Relapse[!(percentages_df$Relapse %in% c("TRUE", "FALSE"))] <- NA
percentages_df$Relapse<- as.logical(percentages_df$Relapse)

plot_survival_curve(percentages_df %>% filter(!is.na(!!as.name("Relapse"))) %>% filter(`Compartment (3)`=="Boundary"), 
                    "RFS (days)", 
                    "Relapse", 
                    "B_high",
                    title = "RFS LLM B medhigh vs low Boundary quartile",
                    calc.p.value = TRUE,
                    filename = file.path(output_folder, paste0(run_name,"_lowerquartile_Boundary.png")))


# Stroma Core 
quantiles <- quantile(percentages_df %>% filter(`Compartment (3)`=="Stroma Core") %>% pull(value), probs = seq(0, 1, 0.25))
percentages_df$B_high <- "MediumHigh"
percentages_df$B_high[percentages_df$value <= quantiles[2]] <- "Low"

percentages_df$Relapse[!(percentages_df$Relapse %in% c("TRUE", "FALSE"))] <- NA
percentages_df$Relapse<- as.logical(percentages_df$Relapse)

plot_survival_curve(percentages_df %>% filter(!is.na(!!as.name("Relapse"))) %>% filter(`Compartment (3)`=="Stroma Core"), 
                    "RFS (days)", 
                    "Relapse", 
                    "B_high",
                    title = "RFS LLM B medhigh low Stroma Core quartile",
                    calc.p.value = TRUE,
                    filename = file.path(output_folder, paste0(run_name, "_lowerquartile_Stroma_Core.png")))


# Tumour Core 
quantiles <- quantile(percentages_df %>% filter(`Compartment (3)`=="Tumour Core") %>% pull(value), probs = seq(0, 1, 0.25))
percentages_df$B_high <- "MediumHigh"
percentages_df$B_high[percentages_df$value <= quantiles[2]] <- "Low"

percentages_df$Relapse[!(percentages_df$Relapse %in% c("TRUE", "FALSE"))] <- NA
percentages_df$Relapse<- as.logical(percentages_df$Relapse)

plot_survival_curve(percentages_df %>% filter(!is.na(!!as.name("Relapse"))) %>% filter(`Compartment (3)`=="Tumour Core"), 
                    "RFS (days)", 
                    "Relapse", 
                    "B_high",
                    title = "RFS LLM B medhigh low Tumour Core quartile",
                    calc.p.value = TRUE,
                    filename = file.path(output_folder, paste0(run_name, "_lowerquartile_Tumour_Core.png")))