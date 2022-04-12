# This function runs a full Seurat analysis: normalization, scaling, clustering, PCA analysis
seuratAnalysis <- function(datasetName,dataset, metadataAvailable = TRUE, mapTypeValue = "umap",doSubset = FALSE, metadata, resolutionValue = 0.1, integrateValue = FALSE, checkPCAs = FALSE,specificPCA = FALSE, ngenes = 200,scaleFactor = 10000,nPCAS=30, selectionFeatures=FALSE){
    # Create Seurat object
    datasetSeurat <- createSeurat(projectName = datasetName,dataset = dataset,integration = integrateValue, ngenes, scaleFactor,nPCAS,selectionFeatures)
    # Create UMAP/Tsne depending on PCA value defined
    if (specificPCA != FALSE) {
      PCAvalue <- specificPCA
      # Generate UMAP
      datasetSeurat <- createMap(seuratObject = datasetSeurat,mapType = mapTypeValue,pca = PCAvalue,resolution = resolutionValue)
      
    }else{
      # Check optimal PCA value
      PCAvalue <- findOptimalPCA(datasetSeurat)
      # Generate UMAP
      datasetSeurat <- createMap(seuratObject = datasetSeurat, mapType = mapTypeValue,pca = PCAvalue,resolution = resolutionValue)
    }
    
    # Add Metadata file if available
    datasetSeurat <- addMetadataSeurat(metadataAvailable,metadata,datasetSeurat)
    
    # Return final clustered file
    datasetSeurat
}