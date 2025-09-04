library(tidyr)
library(dplyr)

# read the cell df
cell_df <- readRDS("../../data/CM_processed_cell_meta_data_df.rds")

# get the intensity columns
intensity_cols <- colnames(cell_df)[grepl("Mean", colnames(cell_df))]

# cell type columns
celltype_col_name <- "Hierarchy: Level: 4"
parent_col_name <- "Hierarchy: Level: 2"

plot_name <- "main_cell_type_marker_heatmap"
output_folder <- "../plots/"
if (!dir.exists(output_folder)){
    dir.create(output_folder, recursive = T)
}


# extract the columns for the heatmap
celltype_intensities_df <- cell_df[, c(celltype_col_name, parent_col_name, intensity_cols)]

# remove Cell: Mean and Nucleus: Mean from column names
colnames(celltype_intensities_df) <- gsub("(: Cell: Mean)|(: Nucleus: Mean)", "", colnames(celltype_intensities_df))

# order of the markers
markers <- c("panCK", "Vimentin", "CD45", "CD20", "CD3", "CD4", "CD8a",  "Foxp3", "CD56", 
             "dTCR", "CD68", "MHC II (HLA-DR)", "CD11c", "CD66b", "Tryptase")
# cell types
celltypes <- c("Epithelial cells", 
               "Stromal cells",
               "B cells",
               "CD4 T cells", 
               "CD8 T cells", 
               "Treg cells",
               "NK cells", 
               "yd T cells", 
               "Macrophages", 
               "Granulocytes", 
               "Mast cells", 
               "Dendritic cells")

# parent type order
parenttypes <- c("Epithelial cells",
                 "Stromal cells",
                 "Lymphoid cells",
                 "Myeloid cells"
)

# get the mean intensities 
mean_intensities <- celltype_intensities_df %>%
    group_by_at(c(celltype_col_name, parent_col_name)) %>%
    summarise(count=n(), across(where(is.numeric), mean)) %>%
    ungroup() %>%
    arrange(!!as.name(parent_col_name), count)  %>%
    as.data.frame()

mean_intensities <- mean_intensities %>% filter(!!as.name(celltype_col_name) %in% celltypes)

long_mean_intensities <- mean_intensities %>%
    pivot_longer(-c(parent_col_name, celltype_col_name, "count"), names_to = "marker", values_to = "mean")

# scale the data
long_mean_intensities <- long_mean_intensities %>%
    group_by(marker) %>%
    mutate(standardised_mean = scale(mean)[,1])

long_mean_intensities[, parent_col_name] <- factor(long_mean_intensities %>% pull(parent_col_name),
                                                   levels=parenttypes)
long_mean_intensities[, celltype_col_name] <- factor(long_mean_intensities %>% pull(celltype_col_name),
                                                     levels=rev(celltypes))
long_mean_intensities[, "marker"] <- factor(long_mean_intensities %>% pull("marker"),
                                            levels=markers)

mean_intensities[, parent_col_name] <- factor(mean_intensities %>% pull(parent_col_name),
                                              levels=parenttypes)
mean_intensities[, celltype_col_name] <- factor(mean_intensities %>% pull(celltype_col_name),
                                                levels=rev(celltypes))

# actual plot
raster_plot <- ggplot(long_mean_intensities %>% filter(marker %in% markers), aes(x=marker, y=!!as.name(celltype_col_name), fill = standardised_mean)) + 
    geom_raster() +
    geom_tile(aes(fill = standardised_mean), colour = "black", linewidth = 0.5) + 
    scale_fill_distiller(palette = "YlGnBu")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(rows = vars(!!as.name(parent_col_name)), scales="free_y", space = "free") +
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


bar_plot <- ggplot(mean_intensities, aes(x=!!as.name(celltype_col_name), y=count)) + 
    geom_bar(stat = "identity") + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.y = element_blank(),
          plot.margin = margin(t = 0,  # Top margin
                               r = 10,  # Right margin
                               b = 0,  # Bottom margin
                               l = -0.8),
          panel.background = element_blank(),
          panel.spacing = unit(0.001, "cm")) + 
    coord_flip(ylim = c(0, 30000)) + 
    facet_grid(rows = vars(!!as.name(parent_col_name)), scales="free_y", space = "free") 




aligned <- cowplot::align_plots(raster_plot, bar_plot, align="h", axis="tblr")


p <- gridExtra::grid.arrange(aligned[[1]], aligned[[2]], ncol=2, nrow=1, widths=c(4, 1))
ggsave(file.path(output_folder, paste0(plot_name, ".png")), p, width = 12, height = 12, dpi = 600)
# theme(panel.spacing.x=unit(0.1, "lines"),panel.spacing.y=unit(1, "lines"))

