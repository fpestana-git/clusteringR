#' Seurat clustering with log or SCT normalization
#' 
#' Function that takes a dataset (either spatial or scRNA-seq), creates Seurat object, runs log normalization, PCA, etc, clustering
#' 
#' @param datasetObject dataset object
#' @param datasetName name of the dataset
#' @param metadataAvailable whether there is prior metadata available (default TRUE)
#' @param normalizationMethod which normalization method to use (default LOG). options log normalization or SCTransform normalization (SCT)
#' @param reductionValue which dimensionality reduction to use, umap or tsne (default umap)
#' @param metadataObject metadata object
#' @param resolutionValue define the resolution value to use for clustering (default 0.1)
#' @param integrateValue whether dataset integration is needed (default FALSE)
#' @param specificPCA whether PCA values needs to be calculated or not
#' @param pcaValueOptimal define a specific PCA value to use (default NULL calculates automatically the optimal PCA value)
#' @param ngenes minimum number of features expressed in a cell (default 200)
#' @param scaleFactor value for scaling (default 10000)
#' @param nPCAS number of principal components to use (default 30)
#' @param selectFeatures whether features need to be selected. Takes all features as default (default NULL)
#'
#' @export clusteringSeurat


# This function runs a full Seurat analysis: normalization, scaling, clustering, PCA analysis
clusteringSeurat <- function(datasetObject,datasetName, metadataAvailable = TRUE,normalizationMethod = "LOG", mapTypeValue = "umap", metadataObject, reductionValue = "umap",resolutionValue = 0.1, pcaValueOptimal = NULL, integrateValue = FALSE,specificPCA = FALSE, ngenes = 200,scaleFactor = 10000,nPCAS=30, selectFeatures = NULL){
  # Create Seurat object
    #datasetSeurat <- createSeurat(projectName = datasetName,dataset = datasetObject,integration = integrateValue, ngenes, scaleFactor,nPCAS,selectionFeatures)
    # Load selected features (if not NULL)
    if (!is.null(selectFeatures)) {
      selectFeatures <- selectedFeatureNames
    }
  
    # Create seurat object
    datasetSeurat <- CreateSeuratObject(counts = datasetObject, project = datasetName,min.cells = 1,min.features = 1)
    # Add name of the dataset
    datasetSeurat$stim <- datasetName

    # Apply different normalization methods 
    if (normalizationMethod == "LOG") {
      datasetSeurat <- NormalizeData(datasetSeurat, normalization.method = "LogNormalize", scale.factor = scaleFactor, verbose = TRUE)
      datasetSeurat <- ScaleData(datasetSeurat, verbose = TRUE, features = rownames(datasetSeurat))
      datasetSeurat <- FindVariableFeatures(datasetSeurat, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
    }else if (normalizationMethod == "SCT"){
      datasetSeurat <- SCTransform(object = datasetSeurat, assay = "RNA", verbose = TRUE,residual.features = selectFeatures,do.scale = T,ncells = 1000)
    }
    
    # Run PCA
    datasetSeurat <- RunPCA(datasetSeurat, npcs = nPCAS, verbose = TRUE)

    # Calculate PCA value automatically
    if (is.null(pcaValueOptimal)) {
      pcaValueOptimal <- findOptimalPCA(seuratObj = datasetSeurat)
    }else{
      pcaValueOptimal <- pcaValueOptimal
    }
    
    # Store the pca and resolution values on the Seurat object
    datasetSeurat@reductions[["pcaValueOptimal"]] <- pcaValueOptimal
    datasetSeurat@reductions[["resolutionValue"]] <- resolutionValue
    
    # Run clustering
    datasetSeurat <- FindNeighbors(object = datasetSeurat, reduction = "pca", dims = 1:pcaValueOptimal, verbose = FALSE)
    datasetSeurat <- FindClusters(object = datasetSeurat, resolution = resolutionValue, verbose = FALSE)
    
    # Run visualization
    if(reductionValue == "umap"){
      datasetSeurat <- RunUMAP(object = datasetSeurat, reduction = "pca", dims = 1:pcaValueOptimal, verbose = FALSE)
    }else if(reductionValue == "tsne"){
      datasetSeurat <- RunTSNE(object = datasetSeurat, reduction = "pca", dims = 1:pcaValueOptimal, verbose = FALSE,check_duplicates=FALSE)
    }
    
    # Add Metadata file if available
    datasetSeurat <- addMetadataSeurat(metadataAvailable,metadataObject,datasetSeurat)
    
    # Return final clustered file
    datasetSeurat
}