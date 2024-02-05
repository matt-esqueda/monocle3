# Load the data into Monocle3

library(monocle3)
library(dplyr)

expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_rowData.rds"))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


# Pre-process the data 
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)


# Reduce dimensionality and visualize the cells
cds <- reduce_dimension(cds)
plot_cells(cds)

plot_cells(cds, color_cells_by = "cao_cell_type")

plot_cells(cds, genes = c("cpna-2", "egl-21", "ram-2", "inos-1"))

cds <- reduce_dimension(cds, reduction_method = "tSNE")
plot_cells(cds, reduction_method = "tSNE", color_cells_by = "cao_cell_type")


# Check and remove batch effects
plot_cells(cds, color_cells_by = "plate", label_cell_groups = FALSE)

cds <- align_cds(cds, num_dim = 100, alignment_group = "plate")
cds <- reduce_dimension(cds)
plot_cells(cds, color_cells_by = "plate", label_cell_groups = FALSE)


# Group cells into clusters
cds <- cluster_cells(cds, resolution = 1e-5)
plot_cells(cds)

plot_cells(cds, color_cells_by = "partition", group_cells_by = "partition")

plot_cells(cds, color_cells_by = "cao_cell_type")

plot_cells(cds, color_cells_by = "cao_cell_type", label_groups_by_cluster = FALSE)


# Find marker genes expressed by each cluster
marker_test_res <- top_markers(cds, group_cells_by = "partition",
                               reference_cells = 1000, cores = 8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_markers_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_markers_ids,
                    group_cells_by = "partition",
                    ordering_type = "maximal_on_diag",
                    max.size = 3)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_markers_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_markers_ids,
                    group_cells_by = "partition",
                    ordering_type = "cluster_row_col",
                    max.size = 3)


# Annotate your cells according to type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))

colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                 "1"="Body wall muscle",
                                                 "2"="Germline",
                                                 "3"="Motor neurons",
                                                 "4"="Seam cells",
                                                 "5"="Sex myoblasts",
                                                 "6"="Socket cells",
                                                 "7"="Marginal_cell",
                                                 "8"="Coelomocyte",
                                                 "9"="Am/PH sheath cells",
                                                 "10"="Ciliated neurons",
                                                 "11"="Intestinal/rectal muscle",
                                                 "12"="Excretory gland",
                                                 "13"="Chemosensory neurons",
                                                 "14"="Interneurons",
                                                 "15"="Unclassified eurons",
                                                 "16"="Ciliated neurons",
                                                 "17"="Pharyngeal gland cells",
                                                 "18"="Unclassified neurons",
                                                 "19"="Chemosensory neurons",
                                                 "20"="Ciliated neurons",
                                                 "21"="Ciliated neurons",
                                                 "22"="Inner labial neuron",
                                                 "23"="Ciliated neurons",
                                                 "24"="Ciliated neurons",
                                                 "25"="Ciliated neurons",
                                                 "26"="Hypodermal cells",
                                                 "27"="Mesodermal cells",
                                                 "28"="Motor neurons",
                                                 "29"="Pharyngeal gland cells",
                                                 "30"="Ciliated neurons",
                                                 "31"="Excretory cells",
                                                 "32"="Amphid neuron",
                                                 "33"="Pharyngeal muscle")

plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type")


# Choose cells for a subset
cds_subset <- choose_cells(cds)

pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3)

plot_cells(cds_subset, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

cds_subset <- cluster_cells(cds_subset, resolution=1e-2)
plot_cells(cds_subset, color_cells_by="cluster")

colData(cds_subset)$assigned_cell_type <- as.character(clusters(cds_subset)[colnames(cds_subset)])
colData(cds_subset)$assigned_cell_type <- dplyr::recode(colData(cds_subset)$assigned_cell_type,
                                                        "1"="Sex myoblasts",
                                                        "2"="Somatic gonad precursors",
                                                        "3"="Vulval precursors",
                                                        "4"="Sex myoblasts",
                                                        "5"="Vulval precursors",
                                                        "6"="Somatic gonad precursors",
                                                        "7"="Sex myoblasts",
                                                        "8"="Sex myoblasts",
                                                        "9"="Ciliated neurons",
                                                        "10"="Vulval precursors",
                                                        "11"="Somatic gonad precursor",
                                                        "12"="Distal tip cells",
                                                        "13"="Somatic gonad precursor",
                                                        "14"="Sex myoblasts",
                                                        "15"="Vulval precursors")

plot_cells(cds_subset, group_cells_by="cluster", color_cells_by="assigned_cell_type")

colData(cds)[colnames(cds_subset),]$assigned_cell_type <- colData(cds_subset)$assigned_cell_type
cds <- cds[,colData(cds)$assigned_cell_type != "Failed QC" | is.na(colData(cds)$assigned_cell_type )]
plot_cells(cds, group_cells_by="partition", 
           color_cells_by="assigned_cell_type", 
           labels_per_group=5)