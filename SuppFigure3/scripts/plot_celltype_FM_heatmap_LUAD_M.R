

library(ggplot2)
library(dplyr)

source("src/proportion/calculate_cell_type_percentages.R")
# ------------------------------------------------------------------------------
# Read the data
# ------------------------------------------------------------------------------
patient_counts <- readRDS(file="data/NSCLCcohort/processed/patient_counts_compartment_3_all_FM_TRM_df.rds")


patient_counts <- patient_counts  %>%
    filter(`Cancer Type`=="LUAD" & `MIBI clusters` == "M")

# group col
group_cols <- c("Patient", "Compartment (3)")


name_of_plot <- paste0(paste("celltype_functional_marker_heatmaps_LUAD_M_compartment_3",   sep="_"), ".png")

output_folder <- "data/NSCLCcohort/output/Figure3/heatmap"
if (!dir.exists(output_folder)){
    dir.create(output_folder, recursive = T)
}


# myeloid
cell_types <- c(
    "Macrophages", 
    "Granulocytes", 
    "Mast cells", 
    "Dendritic cells"
)


# functional markers
non_specific_functional_markers <- c(
    "Ki67",
    "IFN-y",
    "CD16",
    "CD69",
    "CD103",
    "CD49a",
    "Tim3"
    # "MHC I (HLA Class1)",
    # "MHC II (HLA-DR)"
)

myeloid_functional_markers <- c(
    "CD14",
    "CD163",
    "CD206"
)

functional_markers <- c(non_specific_functional_markers,
                        myeloid_functional_markers)


# get all the cell type functional marker combinations
combinations <- expand.grid(cell_types = cell_types, functional_markers = functional_markers)
combinations_list <- paste0(combinations$cell_types, ": ", combinations$functional_markers, "+")

# calculate the proportions
percentages_df <- calculate_cell_type_percentages(patient_counts, group_cols, combinations_list, as.vector(combinations$cell_types))


# remove 
percentages_df <- percentages_df %>%
    filter(parent_type_count >= 5)


# check if the correct order is maintained
print(paste("Cell type and parent pairing:", all(mapply(grepl, percentages_df$parent_type, percentages_df$cell_type))))

# check if correct patient order
print(paste("Correct patient order:", all(percentages_df$Patient == percentages_df$parent_Patient)))



mean_percentages <- percentages_df %>%
    group_by_at(c("Compartment (3)", "percentage", "cell_type", "parent_type")) %>%
    summarise(mean_percentage = mean(value))

# get the corresponding cell type and the functional marker for that column
mean_percentages$functional_marker <- gsub("(.*: )|(\\+)", "", mean_percentages$cell_type)


# add the factor levels to the cell type and functional markers 
mean_percentages$parent_type <- factor(mean_percentages$parent_type, levels=rev(cell_types))
mean_percentages$functional_marker <- factor(mean_percentages$functional_marker, levels = functional_markers)



# standardise the percentages
mean_percentages <- mean_percentages %>%
    group_by_at(c("functional_marker")) %>%
    mutate(standardised_mean = scale(mean_percentage)[,1])

# add the grouping for the functional markers 
mean_percentages$functional_marker_specificity <- "all"

fm_myeloid_bool <- mean_percentages$functional_marker %in% myeloid_functional_markers 
mean_percentages$functional_marker_specificity[fm_myeloid_bool] <- "Myeloid"


mean_percentages$functional_marker_specificity <- factor(mean_percentages$functional_marker_specificity,
                                                         levels = c("all", "Myeloid"))


mean_percentages$`Compartment (3)` <- factor(mean_percentages$`Compartment (3)`, levels = c("Tumour Core", "Boundary", "Stroma Core"))

# plot the heatmap 
ggplot(mean_percentages, aes(x=functional_marker, y=parent_type, fill = standardised_mean)) + 
    geom_raster() +
    geom_tile(aes(fill = standardised_mean), colour = "black", linewidth = 0.2) + 
    destiny::scale_fill_cube_helix(discrete = F, reverse=T)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # facet_wrap(vars(!!as.name(group_col)), nrow = 2) +
    facet_grid(
        # cols = vars(functional_marker_specificity), 
        rows = vars(!!as.name("Compartment (3)")), 
        scales="free", space = "free") +
    theme(
        
        panel.background = element_blank(),
        # strip.text.y = element_blank(),
        legend.position="top",
        plot.margin = margin(t = 0,  # Top margin
                             r = -0.8,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0),
        panel.spacing = unit(0.1, "cm")
    ) +
    ggtitle(paste("Scaled across rows (scaled within func. marker)"))
ggsave(filename=file.path(output_folder, name_of_plot))



# calculate the pvalues
new_col_name <- "Compartment_3"
percentages_df[, new_col_name] <- percentages_df[, "Compartment (3)"]

formula <- as.formula(paste0("value~", new_col_name))

stat.test.excl.boundary <- percentages_df %>%
    filter(`Compartment (3)`!="Boundary") %>%
    group_by_at(c("percentage")) %>%
    group_modify(~{
        tryCatch({
            return(rstatix::wilcox_test(formula = formula, data = .x))},
            error=function(cond){
                return(data.frame())
            }
        )
    })
    
stat.test.excl.tumour_core <- percentages_df %>%
    filter(`Compartment (3)`!="Tumour Core") %>%
    group_by_at(c("percentage")) %>%
    group_modify(~{
        tryCatch({
            return(rstatix::wilcox_test(formula = formula, data = .x))},
            error=function(cond){
                return(data.frame())
            }
        )
    })

stat.test.excl.stroma_core <- percentages_df %>%
    filter(`Compartment (3)`!="Stroma Core") %>%
    group_by_at(c("percentage")) %>%
    group_modify(~{
        tryCatch({
            return(rstatix::wilcox_test(formula = formula, data = .x))},
            error=function(cond){
                return(data.frame())
            }
        )
    })

stat.test <- rbind(stat.test.excl.boundary, 
                   stat.test.excl.stroma_core,
                   stat.test.excl.tumour_core)

write.csv(stat.test, file=file.path(output_folder, paste0(paste("pvalues", name_of_plot, sep="_"), ".csv")))



