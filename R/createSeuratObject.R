# Function that creates a Seurat object and used on seuratAnalysis function
# Crease Seurat object, either for single analysis or considering integration if integration value is set as true
createSeurat <- function(projectName, dataset, integration = FALSE, ngenes, scaleFactor,nPCAS,selectionFeatures){
  # Create seurat object
  datasetSeurat <- CreateSeuratObject(counts = dataset, project = projectName)
  # Add name of the dataset
  datasetSeurat$stim <- projectName
  # Filter cells that express less than ngenes
  datasetSeurat <- subset(datasetSeurat, subset = nFeature_RNA >= ngenes)
  # Apply log normalization
  datasetSeurat <- NormalizeData(datasetSeurat, normalization.method = "LogNormalize", scale.factor = scaleFactor, verbose = TRUE)
  
  # Run PCA if integration is not the goal
  if (integration == FALSE) {
    # Apply data scaling
    datasetSeurat <- ScaleData(datasetSeurat, verbose = TRUE, features = rownames(datasetSeurat))
    if(selectionFeatures == TRUE){
      # If selection of features is needed, run on selected features
      selectedFeatures <- grep(pattern = "^Pcdh|Il33",x = datasetSeurat@assays[["RNA"]]@data@Dimnames[[1]],value = T,invert = T)
      datasetSeurat <- RunPCA(datasetSeurat, features = selectedFeatures, npcs = nPCAS,verbose = TRUE)
    }else{
      # If selection of features is not needed, run findvariablefeatures
      datasetSeurat <- FindVariableFeatures(datasetSeurat, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
      datasetSeurat <- RunPCA(datasetSeurat, npcs = nPCAS, features = VariableFeatures(object = datasetSeurat), verbose = TRUE)
    }
  }
  # Return clustering object
  datasetSeurat
}
