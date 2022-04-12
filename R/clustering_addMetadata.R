# Function that adds metadata to seurat object
addMetadataSeurat <- function(metadataAvailable,metadata,datasetSeurat){
  # Add metadata
  if (metadataAvailable == TRUE) {
    datasetSeurat <- AddMetaData(object = datasetSeurat,metadata = metadata)
  }
  datasetSeurat
}
