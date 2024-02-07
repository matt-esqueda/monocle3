# Load the data into Monocle3

library(monocle3)
library(dplyr)

expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = 
                   "~ bg.300.loading + bg.400.loading + bg.500.1.loading + 
                 bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading") 


# Reduce the dimensionality and visualize the cells
cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "cell.type")

ciliated_genes <- c("che-1",
                    "hlh-17",
                    "nhr-6",
                    "dmd-6",
                    "ceh-36",
                    "ham-1")

plot_cells(cds,
           genes = ciliated_genes,
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE)


# Cluster the cells
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")


# Learn the Trajectory Graph
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE)


# Order the cells in pseudotime
plot_cells(cds,
           color_cells_by = "embryo.time.bin",
           label_cell_groups = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 1.5)


# Choose root nodes
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5)

# a helper function to identify the root principal points
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                            (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5,
           label_principal_points = TRUE)


# Subset cells by branch
cds_sub <- choose_graph_segments(cds)

# working with 3D trajectories
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes = get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by = "partition")

cds_3d_plot_obj
