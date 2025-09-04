
library(ggplot2) # For creating plots
library(Rcpp) # For running UMAP
library(uwot)

cell_df <- readRDS("../../data/CM_processed_cell_meta_data_df.rds")

plot_name <- "main_cell_type_UMAP"

output_folder <- "../plots/"
if (!dir.exists(output_folder)){
    dir.create(output_folder, recursive = T)
}

# random state
set.seed(42)

celltype_col_name <- "Hierarchy: Level: 4"

# get the intensity columns
celltype_intensities_df <- cell_df[, colnames(cell_df)[grepl("Mean", colnames(cell_df))]]

# remove Cell: Mean and Nucleus: Mean from column names
colnames(celltype_intensities_df) <- gsub("(: Cell: Mean)|(: Nucleus: Mean)", 
                                          "", colnames(celltype_intensities_df))

# order of the markers
markers <- c("panCK", "Vimentin", "CD45", "CD3", "CD4", "CD8a", "CD20", "Foxp3", 
             "CD56", "dTCR", "CD68", "MHC II (HLA-DR)", "CD11c", "CD66b", "Tryptase")


celltypes <- c(
    "B cells",
    "CD4 T cells", 
    "CD8 T cells", 
    "Treg cells",
    "NK cells", 
    "yd T cells", 
    "Macrophages", 
    "Granulocytes", 
    "Mast cells", 
    "Dendritic cells",
    "Epithelial cells",
    "Stromal cells"
    )
# (T cells for other plots can be magenta1)

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
    "royalblue1",
    "grey",
    "gold"
    )


celltype_intensities_df <- celltype_intensities_df[, markers]

# filter out only the immune subtypes
filter_bool <- cell_df[, celltype_col_name] %in% celltypes
celltype_intensities_df <- celltype_intensities_df[filter_bool, ]

# run umap
system.time({
    cell.umap <- umap(celltype_intensities_df)
})

cell.umap <- as.data.frame(cell.umap)

cell.umap$cell_type <- cell_df[, celltype_col_name][filter_bool]
cell.umap$cell_type <- factor(cell.umap$cell_type, levels=celltypes)


ggplot(as.data.frame(cell.umap), aes(x=V1, y=V2, color=cell_type)) +
    geom_point(size = 0.5, stroke = 0)+
    theme(
        panel.background = element_blank()) + 
    scale_color_manual(values = colors)+
    guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(file.path(output_folder, paste0(plot_name, ".png")), width = 12, height = 12, dpi = 600)
