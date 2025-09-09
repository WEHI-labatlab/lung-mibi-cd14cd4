
preprocessed_df_ <- readRDS("../interim_data/processed_spatial_clusters_cell_meta_data_df.rds")

output_folder <- "data/NSCLCcohort/output/heatmaps/CNclusters"
if (!dir.exists(output_folder)){
    dir.create(output_folder)
}


# required columns 
hierarchy_level_col <- 'Hierarchy: Level: 4'
x_coordinates <- "Centroid X µm"
y_coordinates <- "Centroid Y µm"
groups <- "Image"

cluster_col <- "CN.clusters"

group_cols <- c()

celltypes <- c("B cells",
               "CD4 T cells", 
               "CD8 T cells", 
               "Treg cells",
               "NK cells", 
               "yd T cells", 
               "Macrophages", 
               "Granulocytes", 
               "Mast cells", 
               "Dendritic cells")

colors <- c(
    "red",
    "mediumorchid1",
    "mediumpurple1",
    "pink1",
    "salmon1",
    "palevioletred1",
    "turquoise1",
    "seagreen1",
    "greenyellow",
    "royalblue1"
)

# summarise the cluster memberships
cohort_cluster_membership_proportions <- preprocessed_df_ %>%
    group_by_at(c(cluster_col, hierarchy_level_col, group_cols)) %>%
    summarise(count=n()) %>%
    ungroup() %>%
    group_by_at(cluster_col)

# calculate the proportions in each cluster
cohort_cluster_membership_proportions <- cohort_cluster_membership_proportions %>% 
    group_by_at(c(cluster_col, group_cols)) %>% 
    group_modify(~{
        proportion <- .x$count / sum(.x$count)
        cbind(.x, proportion)
    }) %>% 
    filter(CN.clusters != -1)

# add the levels 
cohort_cluster_membership_proportions[, hierarchy_level_col] <- factor(cohort_cluster_membership_proportions %>% 
                                                                           pull(hierarchy_level_col),
                                                                       levels = celltypes)

cohort_cluster_membership_proportions <- cohort_cluster_membership_proportions %>%
    group_by(CN.clusters) %>%
    mutate(standardised_proportion = scale(proportion)[,1])

# plot the heat map
ggplot(cohort_cluster_membership_proportions, aes(x=!!as.name(hierarchy_level_col), y=as.factor(CN.clusters), fill = standardised_proportion)) + 
    geom_raster() +
    geom_tile(aes(fill = standardised_proportion), colour = "black", linewidth = 0.2) + 
    scale_fill_distiller(palette = "YlGnBu")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # facet_grid(rows = vars(!!as.name(parent_col_name)), scales="free_y", space = "free") +
    theme(
        panel.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position="top",
        plot.margin = margin(t = 0,  # Top margin
                             r = -0.8,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0),
        panel.spacing = unit(0.001, "cm")
    ) 
ggsave(file.path(output_folder, "CN_membership_heatmap.png"))

ggplot(data=cohort_cluster_membership_proportions, aes(x=as.factor(CN.clusters), y=proportion, fill=!!as.name(hierarchy_level_col))) +
    geom_bar(stat="identity")+
    geom_text(aes(label=round(proportion,2)), size = 3, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = colors)+
    ggtitle("cohort") + 
    theme(axis.title.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          # axis.ticks.y=element_blank(),
          strip.text.y = element_blank(),
          panel.background = element_blank())
ggsave(file.path(output_folder, "CN_membership_proportion_barplot.png"))

ggplot(data=cohort_cluster_membership_proportions, aes(x=as.factor(CN.clusters), y=count, fill=!!as.name(hierarchy_level_col))) +
    geom_bar(stat="identity")+
    geom_text(aes(label=count), size = 3, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = colors)+
    ggtitle("cohort") + 
    theme(axis.title.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          # axis.ticks.y=element_blank(),
          strip.text.y = element_blank(),
          panel.background = element_blank())
ggsave(file.path(output_folder, "CN_membership_histogram.png"))