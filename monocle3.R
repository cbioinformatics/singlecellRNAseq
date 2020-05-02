require(monocle3)
require(ggplot2)
require(dplyr)
require(htmlwidgets)

### Specifying input parameters
date = gsub("2020-","20",Sys.Date(),perl=TRUE);
date = gsub("-","",date)
Dim = "2D"
upn = "Sawamoto_monocle"
input.file <- "integrated_data.rds"
nPC <- 15
cluster.res <- 0.7
output.dir <- sprintf("%s/%s.%s.Seurat.to.Monocle3.%s", getwd(), upn, date, Dim)
dir.create(file.path(output.dir), showWarnings = FALSE)

### Reading in Seurat object
print("Reading in Seurat objects")
seurat <- readRDS(file = input.file)

### Building the necessary parts for a basic cds
# part one, gene annotations
gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information
cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]], row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix
New_matrix <- seurat@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

### Construct the basic cds object
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)
print("Construction Done")
#print(str(cds_from_seurat))

### Construct and assign the made up partition
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
print("Assigned partition")

### Assign the cluster info
#list_cluster <- seurat@meta.data[[sprintf("ClusterNames_%s_%sPC", cluster.res, nPC)]]
list_cluster <- seurat@meta.data[["seurat_clusters"]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
print("Assigned cluster info")

### Assign UMAP coordinate
cds_from_seurat@reduce_dim_aux@listData[["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings
print("Assigned UMAP coordinate")

### Assign feature loading for downstream module analysis
cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings
print("Assigned features")

### Pre-process the data
cds_from_seurat <- preprocess_cds(cds_from_seurat, num_dim=50)
print("Pre-processed the data")

### Reduce dimesionality and visualize the results
cds_from_seurat <- reduce_dimension(cds_from_seurat, umap.fast_sgd=TRUE)
pdf(sprintf("%s/reduce_dimension.pdf", output.dir),width=10,height=10)
clus <- plot_cells(cds_from_seurat,
		   color_cells_by = 'cluster',
		   label_groups_by_cluster=FALSE)
print(clus)
dev.off()
print("Reduced dimensionality")

### Cluster cells
cds_from_seurat <- cluster_cells(cds_from_seurat)
pdf(sprintf("%s/cluster_cells.pdf", output.dir),width=10,height=10)
clus <- plot_cells(cds_from_seurat,
                   color_cells_by = 'partition',
                   label_groups_by_cluster=FALSE)
print(clus)
dev.off()
print("Clustered cells")

### Learn graph, this step usually takes a significant period of time for larger samples
print("Learning graph, which can take a while depends on the sample")
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)
print("Learned graph")

### Ordering cells
root_group = colnames(cds_from_seurat)[clusters(cds_from_seurat)==1]
cds_from_seurat <- order_cells(cds_from_seurat, root_cells = root_group)
print("Ordered cells")

### Reading in cell type information (Not executed)
#print("Reading in cell type information")
#Celltypefile <- "celltype.predictions.July2019.xls"

#lineage.table <- read.table(file = Celltypefile, sep = "", header = T)
#lineage.table$cell.barcodes.i. <- paste(lineage.table$cell.barcodes.i., "-1", sep = "")
#lineage.table <- lineage.table[, 1:2]
#lineage.table$lineage1 <- gsub("_\\d+","", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("TCELLA","T-CELL", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("BCELLA","B-CELL", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("NKA","NK", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("DENDRA","DENDRITIC", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("DENDA","DENDRITIC", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("[0-9]+", "", lineage.table$lineage1)
#lineage.table$cell.barcodes.i. <- gsub("-1", "", lineage.table$cell.barcodes.i.)

#indices <- match(pData(cds_from_seurat)$barcode, lineage.table$cell.barcodes.i.)
#lineage.table.refined <- lineage.table[indices,]
#lineage.table.final <- lineage.table.refined[, 2]
#colData(cds_from_seurat)$celltype <- lineage.table.final

### Plot cluster info with trajectory
print("Plotting clusters")
if (Dim == "2D") {
  pdf(sprintf("%s/clusters.with.trajectory.%s.pdf", output.dir, Dim), width = 10, height = 10)
  clus <- plot_cells(cds_from_seurat, 
                color_cells_by = 'cluster',
                label_groups_by_cluster=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE)
  print(clus)
  dev.off()
}
if (Dim == "3D") {
  cluster <- plot_cells_3d(cds_from_seurat, 
                           color_cells_by = 'cluster',
                           label_groups_by_cluster=FALSE,
                           label_leaves=FALSE,
                           label_branch_points=FALSE)
  saveWidget(cluster, file = sprintf("%s/clusters.with.trajectory.%s.html", output.dir, Dim))
}

### Plot cell type info with trajectory (Not executed)
#print("Plotting cell type info")
#if (Dim == "2D") {
#  pdf(sprintf("%s/celltype.with.trajectory.%s.pdf", output.dir, Dim), width = 10, height = 10)
#  ctype <- plot_cells(cds_from_seurat, 
#                color_cells_by = 'celltype',
#                label_groups_by_cluster=FALSE,
#                label_leaves=FALSE,
#                label_branch_points=FALSE)
#  print(ctype)
#  dev.off()
#}
#if (Dim == "3D") {
#  celltype <- plot_cells_3d(cds_from_seurat, 
#                            color_cells_by = 'celltype',
#                            label_groups_by_cluster=FALSE,
#                            label_leaves=FALSE,
#                            label_branch_points=FALSE)
#  saveWidget(celltype, file = sprintf("%s/celltype.with.trajectory.%s.html", output.dir, Dim))
#}

### Now we plot pseudotime
print("Plotting pseudotime")
if (Dim == "2D") {
  pdf(sprintf("%s/pseudotime.with.trajectory.%s.pdf", output.dir, Dim), width = 10, height = 10)
  ptime <- plot_cells(cds_from_seurat, 
                color_cells_by = 'pseudotime',
                label_groups_by_cluster=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE)
  print(ptime)
  dev.off()
}
if (Dim == "3D") {
  pseudotime <- plot_cells_3d(cds_from_seurat, 
                              color_cells_by = 'pseudotime',
                              label_groups_by_cluster=FALSE,
                              label_leaves=FALSE,
                              label_branch_points=FALSE)
  saveWidget(pseudotime, file = sprintf("%s/pseudotime.with.trajectory.%s.html", output.dir, Dim))
}

print("saving cds to rds")
saveRDS(cds_from_seurat, 'monocle3.rds')

