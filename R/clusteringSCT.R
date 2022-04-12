# Function that takes a dataset (either spatial or scRNA-seq), creates Seurat object, runs SCT normalization, PCA, etc, clustering

clusteringSCT <- function(datasetObject,metadataObject,selectFeatures = NULL,selectedFeatureNames,pcaValueOptimal = NULL, datasetName, labelValue = T,resolutionValue = 0.2){
  # Spatial data
  
  # Select only non pcdh genes for running the clustering of the data
  if (selectFeatures) {
    #selectFeatures <- grep(pattern = "^Pcdh|Il33",x = rownames(datasetObject),invert = T,value = T)
    selectFeatures <- selectedFeatureNames
  }
  dataset1 <- CreateSeuratObject(counts = datasetObject,meta.data = metadataObject)
  dataset1 <- SCTransform(object = dataset1, assay = "RNA", verbose = TRUE,residual.features = selectFeatures,do.scale = T,ncells = 1000)

  dataset1 <- RunPCA(dataset1, assay = "SCT", verbose = TRUE)
  
  if (is.null(pcaValueOptimal)) {
    pcaValueOptimal <- findOptimalPCA(seuratObj = dataset1)
  }else{
    pcaValueOptimal <- pcaValueOptimal
  }
  
  dataset1@reductions[["pcaValueOptimal"]] <- pcaValueOptimal
  dataset1 <- FindNeighbors(dataset1, reduction = "pca", dims = 1:pcaValueOptimal)
  dataset1 <- FindClusters(dataset1, verbose = FALSE, resolution = resolutionValue)
  dataset1 <- RunUMAP(dataset1, reduction = "pca", dims = 1:pcaValueOptimal)
  #resolutionValue <- dataset1@commands[["FindClusters"]]@params[["resolution"]]
  
  # Draw and plot dimplot (both with and without labels)
  #drawDimPlot(dataset = dataset1,mapType = "umap",datasetName = datasetName,pca = pcaValueOptimal,resolution = resolutionValue,labelValue = labelValue)
  #drawDimPlot(dataset = dataset1,mapType = "umap",datasetName = datasetName,pca = pcaValueOptimal,resolution = resolutionValue,labelValue = F)
  
  dataset1
}