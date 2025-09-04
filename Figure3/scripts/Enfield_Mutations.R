install.packages("fst")
library(fst)

df <- read_fst("../interim_data/20221109_TRACERx421_mutation_table.fst")
head(df)  # Check the structure

# Define your list of patient IDs
patient_list <- c("CRUK0001","CRUK0002", "CRUK0003", "CRUK0004","CRUK0005","CRUK0006","CRUK0007", 
                  "CRUK0008","CRUK0009","CRUK0010", "CRUK0012", "CRUK0015", "CRUK0016","CRUK0018",
                  "CRUK0021","CRUK0022","CRUK0023","CRUK0024","CRUK0025", "CRUK0026","CRUK0027","CRUK0028",
                  "CRUK0029","CRUK0034","CRUK0035","CRUK0036","CRUK0037","CRUK0038","CRUK0041","CRUK0042", 
                  "CRUK0045", "CRUK0046","CRUK0048","CRUK0050","CRUK0051","CRUK0052","CRUK0055","CRUK0057", "CRUK0061")  # Replace with actual IDs


# Filter the data frame
df_kras_filtered <- df[grepl("^KRAS", df$Hugo_Symbol) & df$patient_id %in% patient_list, ]
df_egfr_filtered <- df[grepl("^EGFR", df$Hugo_Symbol)& df$patient_id %in% patient_list, ]
df_stk11_filtered <- df[grepl("^STK11", df$Hugo_Symbol)& df$patient_id %in% patient_list, ]

table(df_kras_filtered$AAChange)
# choosing the common mutations with CM data, which are G12 mutations
df_kras_filtered <- df_kras_filtered[df_kras_filtered$AAChange != "p.Q61L",]


# CM
speCancerMIBI_GeoMx <- readRDS("~/cd14_paper/ilaria/claires_part/rdata/speCancerMIBI_GeoMx.rds")
CM_LUAD_immune_subgroups <- readRDS("~/cd14_paper/data/CM_LUAD_immune_subgroups.rds")
library(SpatialExperiment)
x <- colData(speCancerMIBI_GeoMx)
x <- x[x$CancerType == "LUAD",]

x$immuneSubgr <- CM_LUAD_immune_subgroups$cluster[match(x$SampleName, CM_LUAD_immune_subgroups$patient)]
x_filt <- x[,c("SampleName","KRAS", "KRASmut","immuneSubgr", "MibiCluster")]
x_filt <- unique(x_filt)
dim(x_filt)
table(x_filt$KRAS, x_filt$immuneSubgr, useNA = "ifany")
cm_data <- table(x_filt$KRAS, x_filt$immuneSubgr)

CM_LUAD_immune_subgroups$kras <- x_filt$KRAS[match(CM_LUAD_immune_subgroups$patient, x_filt$SampleName)]

# update immune subgroups in enfield
Enfield_LUAD_immune_subgroups <- readRDS("~/cd14_paper/data/Enfield_LUAD_immune_subgroups.rds")

Enfield_LUAD_immune_subgroups$Mutation <- df_kras_filtered$Hugo_Symbol[match(Enfield_LUAD_immune_subgroups$patient,df_kras_filtered$patient_id)]
Enfield_LUAD_immune_subgroups$Mutation <- ifelse(is.na(Enfield_LUAD_immune_subgroups$Mutation), "KRAS_Negative",
                                                 Enfield_LUAD_immune_subgroups$Mutation)
enf_data <- table(Enfield_LUAD_immune_subgroups$Mutation, Enfield_LUAD_immune_subgroups$cluster)

# do stats
fisher.test(cm_data + enf_data)
