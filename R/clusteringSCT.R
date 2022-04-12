#' Seurat clustering with SCT transformation
#' 
#' Function that takes a dataset (either spatial or scRNA-seq), creates Seurat object, runs SCT normalization, PCA, etc, clustering
#' 
#' @param datasetObject dataset object
#' @param metadataObject metadata object
#' @param selectFeatures choose if features are to be selected (default NULL)
#' @param selectedFeatureNames list of features to use on the clustering
#' @param pcaValueOptimal define a specific PCA value to use (default NULL). Otherwise it will calculate an optimal PCA value.
#' @param datasetName define the name of the dataset
#' @param resolutionValue define the resolution value to use for clustering (default 0.2)

clusteringSCT <- function(datasetObject,metadataObject,selectFeatures = NULL,selectedFeatureNames,pcaValueOptimal = NULL, datasetName,resolutionValue = 0.2){
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
  
  # Add resolution value to the seurat object
  spClustering_SCT_CTX1@reductions[["resolutionValue"]] <- resolutionValue
  #resolutionValue <- dataset1@commands[["FindClusters"]]@params[["resolution"]]
  
  # Draw and plot dimplot (both with and without labels)
  #drawDimPlot(dataset = dataset1,mapType = "umap",datasetName = datasetName,pca = pcaValueOptimal,resolution = resolutionValue,labelValue = labelValue)
  #drawDimPlot(dataset = dataset1,mapType = "umap",datasetName = datasetName,pca = pcaValueOptimal,resolution = resolutionValue,labelValue = F)
  
  dataset1
}