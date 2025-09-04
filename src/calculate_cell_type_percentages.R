
library(dplyr)
library(tidyr)

calculate_cell_type_percentages <- function(patient_counts,
                                            group_cols,
                                            cell_types, 
                                            parent_types, 
                                            filter_types = NULL,
                                            keep_group_cols = T
                                            ){
    
    if ((length(cell_types) != length(parent_types))){
        print("The length of the parent_types is not equal to the cell types.
                Either set it as a string or make the length the same as the cell types.
              Will take the first element of parent_types as the parent for all cell types.") 
        parent_types <- rep(parent_types[[1]], length(cell_types))
    }
    
    if ((length(cell_types) != length(filter_types))){
        print("The length of the parent_types is not equal to the cell types.
                Either set it as a string or make the length the same as the cell types.
              Will take the first element of parent_types as the parent for all cell types.") 
        filter_types <- rep(filter_types[[1]], length(cell_types))
    }
    
    ##### Extract the counts of the child, denominator and the filter counts.
    missing_cell_types <- cell_types[!(cell_types %in% colnames(patient_counts) )]
    patient_counts[missing_cell_types] <- 0
    
    missing_parent_types <- cell_types[!(parent_types %in% colnames(patient_counts) )]
    patient_counts[missing_parent_types] <- 0
    
    missing_filter_types <- cell_types[!(filter_types %in% colnames(patient_counts) )]
    patient_counts[missing_filter_types] <- 0
    
    # child type
    cell_type_counts <- patient_counts[c(group_cols, cell_types)] %>%
        pivot_longer(all_of(cell_types), names_to = "cell_type", values_to = "cell_type_count")
    
    # parent type
    parent_type_counts <- patient_counts[c(group_cols, parent_types)] %>%
        pivot_longer(all_of(make.unique(parent_types)), names_to = "parent_type", values_to = "parent_type_count")
    parent_type_counts$parent_type <- gsub("\\.[0-9]+", "", parent_type_counts$parent_type)
    
    # keep the grouping cols
    if (keep_group_cols){
        p_group_cols <- paste0("parent_", group_cols)
        colnames(parent_type_counts) <- c(p_group_cols, "parent_type", "parent_type_count")
    } else{
        parent_type_counts <- parent_type_counts[c("parent_type", "parent_type_count")]
    }
    
    long_percentages <- cbind(cell_type_counts, parent_type_counts)
    
    
    if (!is.null(filter_types)){
        # filter type
        filter_type_counts <- patient_counts[c(group_cols, filter_types)] %>%
            pivot_longer(all_of(make.unique(filter_types)), names_to = "filter_type", values_to = "filter_type_count")
        filter_type_counts$filter_type <- gsub("\\.[0-9]+", "", filter_type_counts$filter_type)
        
        if (keep_group_cols){
            f_group_cols <- paste0("filter_", group_cols)
            colnames(filter_type_counts) <- c(f_group_cols, "filter_type", "filter_type_count")
        } else{
            filter_type_counts <- filter_type_counts[c("filter_type", "filter_type_count")]
        }
        long_percentages <- cbind(long_percentages, filter_type_counts)
        
    } 
    
    long_percentages$value <- (long_percentages$cell_type_count/long_percentages$parent_type_count) * 100
    long_percentages$percentage <- paste0("(", long_percentages$cell_type, "/", long_percentages$parent_type, ")%")
    
    return(long_percentages[!duplicated(long_percentages), ])
}

