# Function to create UMAP plot
createMap <- function(seuratObject = datasetSeurat, mapType, pca, resolution){
  datasetSeurat <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:pca, verbose = FALSE)
  datasetSeurat <- FindClusters(datasetSeurat, resolution = resolution, verbose = FALSE)
  if(mapType == "umap"){
    datasetSeurat <- RunUMAP(object = datasetSeurat, reduction = "pca", dims = 1:pca, verbose = FALSE)
  }else if(mapType == "tsne"){
    datasetSeurat <- RunTSNE(object = datasetSeurat, reduction = "pca", dims = 1:pca, verbose = FALSE,check_duplicates=FALSE)
  }
  datasetSeurat
}  