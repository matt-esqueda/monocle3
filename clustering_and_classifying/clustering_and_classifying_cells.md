Clustering and Classifying Cells
================
Matthew Esqueda
2024-02-05

Single-cell experiments are often performed on tissues containing many
cell types. Monocle 3 provides functions to group cells according to
their gene expression profiles into clusters. Often, cells form clusters
that correspond to once cell type or a set ig highly related cell types.
Monocle uses techniques that are widely accepted in single-cell RNA-seq
analysis, similar to `Seurat`, `scanpy`, and other tools.

Demo the main functions used for clustering with the C. elegans data
from [Cao & Packer et
al.](https://www.science.org/doi/10.1126/science.aam8940). This study
described how to do single-cell RNA-seq with combinatorial indexing in a
protocol called “sci-RNA-seq”, used to produce the first single-cell
RNA-seq analysis of a whole animal, so there are many cell types
represented in the data.

## Load the data

``` r
library(monocle3)
library(dplyr)

expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_rowData.rds"))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
```

## Pre-process the data

Normalize the data, use either `PCA` (the standard for RNA-seq) or
`Latent Semantic Indexing` (common in ATAC-seq), and remove any batch
effects. Here, use the standard PCA method, specifying the number of
principal components you want Monocle to compute.

``` r
cds <- preprocess_cds(cds, num_dim = 100)
```

Check if using enough PCs to capture most of the variation in gene
expression across all the cells in the data set. Look at the fraction of
variation explained by each PC using `plot_pc_variance_explained()`.

``` r
plot_pc_variance_explained(cds)
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Using more than 100 PCs would capture only a small amount of additional
variation, each additional PC makes downstream steps slower.

## Reduce dimensionality and visualize the cells

Use either `t-SNE`, which is very popular in single-cell RNA-seq, or
`UMAP`, which is increasingly common. UMAP is used by default, it is
faster and better suited for clustering and trajectory analysis in
RNA-seq. To reduce the dimensionality of the data down into the X,Y
plane so we can plit it easily, call `reduce_dimension()`.

``` r
cds <- reduce_dimension(cds)
```

To plot the data, use Monocle’s main plotting function, `plot_cells()`.

``` r
plot_cells(cds)
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Each point in the plot represents a different cell in the
`cell_data_set` object `cds`. The cells form many groups, some with
thousands of cells, some with only a few. Cao & Packer annotated each
cell according to type, manually by looking at which genes it expresses.
Color the cells in the UMAP plot by the authors’ orignal annotations
using the `color_cells_by` arg to `plot_cells()`.

``` r
plot_cells(cds, color_cells_by = "cao_cell_type")
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Many of the cell types land very close to one another in the UMAP plot.

Except for a few cases, `color_by_cells` can be the name of any col in
`colData(cds)`. Note that when `color_cells_by` is a categorical
variable, labels are added to the plot, with each label positioned
roughly in the middle of all the cells that have that label.

Can also color cells according to how much of a gene or set of genes
they express

``` r
plot_cells(cds, genes = c("cpna-2", "egl-21", "ram-2", "inos-1"))
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Using t-SNE to visualize the data with `reduction_method="tSNE`.

``` r
cds <- reduce_dimension(cds, reduction_method = "tSNE")
```

Call `plot_cells()`, passing `reduction_method="tSNE` to it as well

``` r
plot_cells(cds, reduction_method = "tSNE", color_cells_by = "cao_cell_type")
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Can you use UMAP and t-SNE on the same `cds` obj - one won’t overwrite
the results of the other. Must specify which one you want in downstream
functions like `plot_cells`.

## Check and remove batch effects

It is important to check for batch effects, systematic differences in
the transcriptome of cells measured in different experimental batches.
These could be technical in nature, such as those introduced during the
single-cell RNA-seq protocol, or biological, such as those that might
arise from different litters of mice.

Check for batch effects when performing dimensionality reduction, Add a
col to the `colData` that encodes which batch each cell is from, then
color the cells by batch. Cao & Packer et al included a “plate”
annotation in their data, which specifies which sci-RNA-seq plate each
cell originated from.

``` r
plot_cells(cds, color_cells_by = "plate", label_cell_groups = FALSE)
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Dramatic bathc effects are not evident in this data. If so, there would
be groups of cells that really ony come from one plate. Nevertheless,
try and remove what batch effect is present by running `align_cds()`.

``` r
cds <- align_cds(cds, num_dim = 100, alignment_group = "plate")
cds <- reduce_dimension(cds)
plot_cells(cds, color_cells_by = "plate", label_cell_groups = FALSE)
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

When ruun the `alignment_group` arg, `align_cds()` tries to remove batch
effects using mutual nearest neighbor alignment. Monocle 3 calls the
package bachelor.

## Group cells into clusters

Monocle uses a technique called `community detection` to group cells for
identifying the cell types represented in the data. This approach was
introduced by Levine et al as part of the pheniGraph algorithm. Cluster
cells using `cluster_cells()`.

``` r
cds <- cluster_cells(cds, resolution = 1e-5)
plot_cells(cds)
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

The `cluster_cells()` also divides the cells into larger, more well
separated groups called partitions, using a stastical test from Alex
Wolf et al, introduced as part of their `PAGA` algorithm.

``` r
plot_cells(cds, color_cells_by = "partition", group_cells_by = "partition")
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Color the cells according to their cell type annotation, and each
cluster is labeled according to the most common annotation within it.

``` r
plot_cells(cds, color_cells_by = "cao_cell_type")
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Choose to label whole partitions instead of clusters by
`passing group_cells_by="partition"`. Plot the top 2 labels per cluster
by passing `labels_per_group=2` to `plot_cells()`. Disable this labeling
policy, making `plot_cells()` behave like it did before calling
`cluster_cells()`.

``` r
plot_cells(cds, color_cells_by = "cao_cell_type", label_groups_by_cluster = FALSE)
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Find marker genes expressed by each cluster

Once cells have been clustered, ask what genes makes them different from
one another with `top_markers()`.

``` r
marker_test_res <- top_markers(cds, group_cells_by = "partition",
                               reference_cells = 1000, cores = 8)
```


The df `marker_test_res` contains a number of metrics for how
specifically expressed each gene is in each partition. Group the cells
according to cluster, partition, or any categorical variable in
`colData(cds)`. Rank the table according to one or more of the
specificity metrics and take the top gene for each cluster. Rank
according to `pseudo_R2`.

``` r
top_specific_markers <- marker_test_res %>%
                            filter(fraction_expressing >= 0.10) %>%
                            group_by(cell_group) %>%
                            top_n(1, pseudo_R2)

top_specific_markers_ids <- unique(top_specific_markers %>% pull(gene_id))
```

Plot the expression and fraction of cells that express each marker in
each group with `plot_genes_by_group`.

``` r
plot_genes_by_group(cds,
                    top_specific_markers_ids,
                    group_cells_by = "partition",
                    ordering_type = "maximal_on_diag",
                    max.size = 3)
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

It is often informative to look at more than one marker, change the
first argument to `top_n()`.

``` r
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
```

![](clustering_and_classifying_cells_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

<!-- ## Annotate cells according to type -->
<!-- Assign cell types based on clustering by creating a new col in `colData(cds)` and initialize it with the values of `partitions(cds)` (can also use cluster(cds) depending on the dataset). -->
<!-- ```{r} -->
<!-- colData(cds)$assigned_cell_type <- as.character(partitions(cds)) -->
<!-- ``` -->
<!-- Use the `dplyr` package's `recode()` to remap each cluster to a differnt cell type. -->
<!-- ```{r} -->
<!-- colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type, -->
<!--                                                  "1"="Body wall muscle", -->
<!--                                                  "2"="Germline", -->
<!--                                                  "3"="Motor neurons", -->
<!--                                                  "4"="Seam cells", -->
<!--                                                  "5"="Sex myoblasts", -->
<!--                                                  "6"="Socket cells", -->
<!--                                                  "7"="Marginal_cell", -->
<!--                                                  "8"="Coelomocyte", -->
<!--                                                  "9"="Am/PH sheath cells", -->
<!--                                                  "10"="Ciliated neurons", -->
<!--                                                  "11"="Intestinal/rectal muscle", -->
<!--                                                  "12"="Excretory gland", -->
<!--                                                  "13"="Chemosensory neurons", -->
<!--                                                  "14"="Interneurons", -->
<!--                                                  "15"="Unclassified eurons", -->
<!--                                                  "16"="Ciliated neurons", -->
<!--                                                  "17"="Pharyngeal gland cells", -->
<!--                                                  "18"="Unclassified neurons", -->
<!--                                                  "19"="Chemosensory neurons", -->
<!--                                                  "20"="Ciliated neurons", -->
<!--                                                  "21"="Ciliated neurons", -->
<!--                                                  "22"="Inner labial neuron", -->
<!--                                                  "23"="Ciliated neurons", -->
<!--                                                  "24"="Ciliated neurons", -->
<!--                                                  "25"="Ciliated neurons", -->
<!--                                                  "26"="Hypodermal cells", -->
<!--                                                  "27"="Mesodermal cells", -->
<!--                                                  "28"="Motor neurons", -->
<!--                                                  "29"="Pharyngeal gland cells", -->
<!--                                                  "30"="Ciliated neurons", -->
<!--                                                  "31"="Excretory cells", -->
<!--                                                  "32"="Amphid neuron", -->
<!--                                                  "33"="Pharyngeal muscle") -->
<!-- ``` -->
<!-- See how the annotations looks: -->
<!-- ```{r} -->
<!-- plot_cells(cds, group_cells_by = "partition", color_cells_by = "assigned_cell_type") -->
<!-- ``` -->
<!-- Partition 7 has some substructure and it;s not obvious just from looking at the output of `top_markers()` what cell type or types it corresponds to. Isolate it with `choose_cells()` for further analysis. -->
<!-- ```{r} -->
<!-- cds_subset <- choose_cells(cds) -->
<!-- ``` -->
<!-- Use `graph_test()` to identify genes that are differentially expressed in different subsets of cells from this partition. -->
<!-- ```{r} -->
<!-- pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8) -->
<!-- pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value <0.05)) -->
<!-- ``` -->
<!-- Take all the genes that vary across this set of cells and group thos that have similar patterns of expression into modules. -->
<!-- ```{r} -->
<!-- gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3) -->
<!-- ``` -->
<!-- Plot these modules; aggregate expression values to reveal which cells express which modules. -->
<!-- ```{r} -->
<!-- plot_cells(cds_subset,genes = gene_module_df, -->
<!--            show_trajectory_graph=FALSE, -->
<!--            label_cell_groups=FALSE) -->
<!-- ``` -->
<!-- Explore the genes in each module or conduct [gene ontology enrichment analysis](https://cole-trapnell-lab.github.io/monocle3/docs/clustering/) on them to glean insights about which cell types are present. This might provide an idea about what the cell types in the partition are. Recluster the cells at finer resolution and then see how they overlap with the clusters in the partition. -->
<!-- ```{r} -->
<!-- cds_subset <- cluster_cells(cds_subset, resolution = 1e-2) -->
<!-- plot_cells(cds_subset, color_cells_by = "cluster") -->
<!-- ``` -->
<!-- Based on how the patterns line up, make the following assignments. -->
<!-- ```{r} -->
<!-- colData(cds_subset)$assigned_cell_type <- as.character(clusters(cds_subset)[colnames(cds_subset)]) -->
<!-- colData(cds_subset)$assigned_cell_type <- dplyr::recode(colData(cds_subset)$assigned_cell_type, -->
<!--                                                         "1"="Sex myoblasts", -->
<!--                                                         "2"="Somatic gonad precursors", -->
<!--                                                         "3"="Vulval precursors", -->
<!--                                                         "4"="Sex myoblasts", -->
<!--                                                         "5"="Vulval precursors", -->
<!--                                                         "6"="Somatic gonad precursors", -->
<!--                                                         "7"="Sex myoblasts", -->
<!--                                                         "8"="Sex myoblasts", -->
<!--                                                         "9"="Ciliated neurons", -->
<!--                                                         "10"="Vulval precursors", -->
<!--                                                         "11"="Somatic gonad precursor", -->
<!--                                                         "12"="Distal tip cells", -->
<!--                                                         "13"="Somatic gonad precursor", -->
<!--                                                         "14"="Sex myoblasts", -->
<!--                                                         "15"="Vulval precursors") -->
<!-- plot_cells(cds_subset, group_cells_by="cluster", color_cells_by="assigned_cell_type") -->
<!-- ``` -->
<!-- Transfer the annotations from the `cds_subset` obj back to the full dataset. Filter out low-quality cells at this stage. -->
<!-- ```{r} -->
<!-- colData(cds)[colnames(cds_subset),]$assigned_cell_type <- colData(cds_subset)$assigned_cell_type -->
<!-- cds <- cds[,colData(cds)$assigned_cell_type != "Failed QC" | is.na(colData(cds)$assigned_cell_type )] -->
<!-- plot_cells(cds, group_cells_by="partition",  -->
<!--            color_cells_by="assigned_cell_type",  -->
<!--            labels_per_group=5) -->
<!-- ``` -->