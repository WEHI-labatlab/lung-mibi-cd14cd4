#' @author Big K
#' 
#' 
#'
#' @param 
#'
#' @return
#' @examples

library(ggplot2)
library(dplyr)
source("src/plot_survival_curve.R")

inflamed_meta_data_path <- "data/meta-data/LUSC_MyeloLympho-Inflamed.xlsx"
inflamed_meta_data_df <- data.frame(readxl::read_excel(inflamed_meta_data_path, 
                                                       sheet="Sheet1"), 
                                    check.names=FALSE)



# Read the percentages
patient_counts_path <- "data/processed/patient_counts_cancertype_compartment_simple_df.rds"
patient_counts <- readRDS(patient_counts_path)

LUSC_patient_counts <- patient_counts %>% filter(`Cancer Type` == "LUSC")

inflamed_meta_data_df <- inflamed_meta_data_df[match(LUSC_patient_counts$`Patient ID`, 
                                                     inflamed_meta_data_df$`Patient`), 
                                               colnames(inflamed_meta_data_df)]



time_name <- "OS (days)"
event_name <- "Vital status"
group_name <- "ClusterMIBI"
new_group_name <- "Inflamation"

title <- "LUSC"

# Convert the days into a numeric
inflamed_meta_data_df[, time_name] <- as.numeric(inflamed_meta_data_df[, time_name])

# Convert it the column to logical
alive_bool <- inflamed_meta_data_df$`Vital status` == "Alive"
alive_maybe_bool <- inflamed_meta_data_df$`Vital status` == "Alive at last record"
dead_bool <- inflamed_meta_data_df$`Vital status` == "Deceased"

inflamed_meta_data_df[, event_name][alive_bool|alive_maybe_bool] <- TRUE
inflamed_meta_data_df[, event_name][dead_bool] <- FALSE
inflamed_meta_data_df[, event_name] <- !as.logical(inflamed_meta_data_df[, event_name])

# Rename column to something more logical
inflamed_meta_data_df[, new_group_name] <- inflamed_meta_data_df[, group_name]


LUSC_patient_counts[, time_name] <- inflamed_meta_data_df[, time_name]
LUSC_patient_counts[, event_name] <- inflamed_meta_data_df[, event_name]
LUSC_patient_counts[, new_group_name] <- inflamed_meta_data_df[, new_group_name]

# ``````````````````````````````````````````````````````````````````````````````
# Separate patients by above and below median of Immune cells%
# ``````````````````````````````````````````````````````````````````````````````
above_median_colname <- "medianCD45"
percentage_of_interest <- "(Immune cells)%"

# get below and above
LUSC_patient_counts[, above_median_colname] <- LUSC_patient_counts[, percentage_of_interest] > median(LUSC_patient_counts[, percentage_of_interest])

# title
title <- "LUSC Overall Survival > median([CD45%])"

# plot the kaplan meier curve
plot_survival_curve(LUSC_patient_counts, 
                    time_name, 
                    event_name, 
                    above_median_colname,
                    title = title
)



# ``````````````````````````````````````````````````````````````````````````````
# Order by Immune cells% then take top5 and bottom 5.
# ``````````````````````````````````````````````````````````````````````````````
ordering <- order(LUSC_patient_counts[, percentage_of_interest])
LUSC_patient_counts$bin <- "Middle"
LUSC_patient_counts[head(ordering, 5),]$bin <- "Bottom5"
LUSC_patient_counts[tail(ordering, 5),]$bin <- "Top5"


time_name <- "OS (days)"
event_name <- "Vital status"

# title
title <- "LUSC Overall Survival Top/Bottom 5 (CD45)%"
# plot the kaplan meier curve
plot_survival_curve(LUSC_patient_counts, 
                    time_name, 
                    event_name, 
                    "bin",
                    title = title
)
