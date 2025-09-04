library(pheatmap)

plot_immune_subgr_heatmap <- function(percentage_df, 
                                      kmean_column, 
                                      percentage_matrix, 
                                      cohort, 
                                      colors,
                                      output_folder){
  # annotation for columns of a heatmap
  row_annot <- data.frame(kmean = percentage_df[[kmean_column]],
                          row.names = rownames(percentage_matrix))
  ann_colors = list(kmean = colors)
  
  p <- pheatmap(t(percentage_matrix[order(row_annot$kmean),]), 
                annotation_col = row_annot, annotation_colors = ann_colors,
                main = paste0(length(colors)," k-mean clustering order ",cohort," cohort"), 
                cluster_cols = F, cluster_rows = F, angle_col = "315",
                color = rev(colorRampPalette(brewer.pal(n = 7, name = "YlGnBu"))(100)))
  p
  
  plot_name <- paste0(output_folder, cohort,"_Heatmap_ImmuneSubgroups_kmean_",length(colors), ".jpg")
  ggsave(plot = p, filename = plot_name, 
         width = 1500, height = 800, units = "px", dpi = 150)
  
  # to preserve the order of patients for the cell type heatmap
  order_of_patients <- data.frame(cluster = row_annot$kmean, patient = rownames(percentage_matrix))
  order_of_patients <- order_of_patients[order(order_of_patients$cluster),]
  return(order_of_patients)
}



plot_cell_type_heatmap <- function(percentages_df,
                                   CancerType,
                                   order_of_patients,
                                   cohort,
                                   colors,
                                   output_folder){
  # we need to filter out the percentages with less than 5 cells
  percentages_df <- percentages_df %>% filter(parent_type_count >= 5)
  
  percentages_df$percentage <- factor(percentages_df$percentage, levels = unique(percentages_df$percentage))
  percentages_df_LUAD <- percentages_df %>% 
    select(Patient, ImmuneSubgroups, value, percentage)
  
  wide_df <- spread(percentages_df_LUAD, key = percentage, value = value)
  wide_df <- wide_df[match(order_of_patients$patient, wide_df$Patient),]
  
  pheatmap_data <- data.frame(t(wide_df[, startsWith(colnames(wide_df), "(")]))
  colnames(pheatmap_data) <- wide_df$Patient
  
  row_annot <- data.frame(kmean = wide_df$ImmuneSubgroups,row.names = wide_df$Patient)
  ann_colors = list(kmean = colors)
  
  p <- pheatmap(pheatmap_data, 
                annotation_col = row_annot, 
                annotation_colors = ann_colors, breaks = seq(0,90, length.out = 70), #breaks fixing the scale, so different range of input data would create the same color scale
                main = paste0(cohort, " cohort ", CancerType), 
                cluster_cols = F, cluster_rows = F, angle_col = "315",
                color = rev(colorRampPalette(brewer.pal(n = 5, name = "YlGnBu"), bias = 0.5)(70)))
  
  p
  plot_name <- paste0(output_folder, cohort,"_",CancerType,
                      "_Heatmap_cell_types_percentage", ".jpg")
  ggsave(plot = p, filename = plot_name, 
         width = 1400, height = 800, units = "px", dpi = 150)
}
