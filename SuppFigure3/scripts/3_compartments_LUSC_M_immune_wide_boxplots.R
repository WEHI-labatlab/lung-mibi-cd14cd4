# library(ggpattern)
# library(ggpmisc)
library(dplyr)
library(tidyr)

output_folder <- "data/NSCLCcohort/output/Figure3_boxplots"
if (!dir.exists(output_folder)){
    dir.create(output_folder)
}

# read the patient counts
patient_counts <- readRDS(file="data/NSCLCcohort/processed/patient_counts_compartment_3_all_FM_TRM_df.rds")
patient_counts <- patient_counts %>% filter(`Cancer Type`=="LUSC", `MIBI clusters`=="M")

compartment_col_name <- "Compartment (3)"

run_name <- "LUSC_M_3_compartments_Myeloid_cell_over_Immune_cells_boxlpot"

celltypes <- c(
    "Macrophages", 
    "Granulocytes", 
    "Mast cells", 
    "Dendritic cells")

colors <- c(
    "turquoise1",
    "seagreen1",
    "greenyellow",
    "royalblue1"
)

denominators <- c(
    "Immune cells",
    "Immune cells",
    "Immune cells",
    "Immune cells"
    )


child_counts <- patient_counts[c("Patient", "Cancer Type",  compartment_col_name, celltypes)] %>%
    pivot_longer(all_of(celltypes), names_to = "cell_type", values_to = "cell_type_count")
denominator_counts <- patient_counts[c("Patient","Cancer Type",  compartment_col_name, denominators)] %>%
    pivot_longer(all_of(make.unique(denominators)), names_to = "denominator_type", values_to = "denominator_count")
denominator_counts$denominator_type <- gsub("\\.[1-9]", "", denominator_counts$denominator_type)

long_percentages_ <- cbind(child_counts, denominator_counts[-c(1:3)])

long_percentages_$value <- (long_percentages_$cell_type_count/long_percentages_$denominator_count) * 100
long_percentages_$percentage <- paste0("(", long_percentages_$cell_type, "/", long_percentages_$denominator_type, ")%")


# need to filter out the patients which have less than K denominator cell types
min_number_denominator <- 5

filtered_out <- long_percentages_ %>%
    filter(denominator_count < min_number_denominator)

filtered_out_summary <- filtered_out %>%
    group_by_at(c("percentage", compartment_col_name)) %>%
    summarise(n())

long_percentages <- long_percentages_ %>%
    filter(denominator_count >= min_number_denominator)

kept_summary <- long_percentages %>%
    group_by_at(c("percentage", compartment_col_name)) %>%
    summarise(n())


write.table(filtered_out, file=file.path(output_folder, paste0(run_name, "_filtered_out_if_denominator_count_less_than_", min_number_denominator,".csv")),  sep = ",", row.names = F)
write.table(kept_summary, file=file.path(output_folder, paste0(run_name, "_kept_if_denominator_count_more_or_equal_", min_number_denominator, ".csv")),  sep = ",", row.names = F)



comparisons <- list(c("Tumour Core", "Boundary"), c("Boundary", "Stroma Core"), c("Tumour Core", "Stroma Core"))
long_percentages$`Compartment (3)` <- factor(long_percentages$`Compartment (3)`,
                                             levels = c("Tumour Core", "Boundary", "Stroma Core"))


dodge <- 1
#, fill=proportion, alpha=`Compartment (3)`

long_percentages$percentage <- factor(long_percentages$percentage, levels = unique(long_percentages$percentage))


ggplot(long_percentages %>% filter(!is.na(!!as.name(compartment_col_name))), aes(x = !!as.name(compartment_col_name), y=value)) +
    # stat_boxplot(geom = "errorbar", position = position_dodge(1))+ 
    # The actual box plot
    geom_boxplot(aes(fill=percentage, alpha=!!as.name(compartment_col_name)), position = position_dodge(dodge))+
    # The jitter plot
    geom_jitter(aes(fill=percentage, shape=!!as.name(compartment_col_name)), alpha=1, size=3, position = position_dodge(dodge)) +
    scale_shape_manual(values = c( 19, 10 ,1 , 2))+
    scale_fill_manual(values=colors) +
    scale_color_manual(values=c("black", "grey32", "grey67", "grey")) +
    scale_alpha_manual(values=c(1, 0.5, 0.1, 0.05)) +
    # ggbreak::scale_y_break(c(45, 70)) + 
    theme(
        axis.title.y=element_blank(),
        # axis.text.y=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        # axis.ticks.y=element_blank(),
        # strip.text.y = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        # panel.spacing = unit(2, "lines")
    ) +
    facet_wrap(~ percentage, nrow = 1, scales = "free_x") +
    ggpubr::stat_compare_means(comparisons = comparisons, hide.ns = F, method="wilcox.test", paired = F, exact = F)
ggsave(file.path(output_folder, paste0(run_name, "_pvalue.png")))


ggplot(long_percentages, aes(x = !!as.name(compartment_col_name), y=value)) +
    # stat_boxplot(geom = "errorbar", position = position_dodge(1))+ 
    # The actual box plot
    geom_boxplot(aes(fill=percentage, alpha=!!as.name(compartment_col_name)), position = position_dodge(dodge))+
    # The jitter plot
    geom_jitter(aes(fill=percentage, shape=!!as.name(compartment_col_name)), alpha=1, size=2, position = position_dodge(dodge)) +
    scale_shape_manual(values = c( 19, 10 ,1 ))+
    scale_fill_manual(values=colors) +
    scale_color_manual(values=c("black", "grey32", "grey67")) +
    scale_alpha_manual(values=c(1, 0.5, 0.1)) +
    # ggbreak::scale_y_break(c(42, 73)) + 
    scale_y_continuous(limits = c(0, 100))+ 
    # ggbreak::scale_y_break(c(80, 90)) + 
    theme(
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_blank(),
        # axis.text.y=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        # axis.ticks.y=element_blank(),
        # strip.text.y = element_blank(),
        # strip.text.x = element_blank(),
        # axis.line.y = element_line(),
        # axis.line.x = element_line(),
        # panel.spacing = unit(2, "lines")
    ) +
    facet_wrap(~ percentage, nrow = 1, scales = "free_x") 
ggsave(file.path(output_folder, paste0(run_name, "_boxplots.png")))


