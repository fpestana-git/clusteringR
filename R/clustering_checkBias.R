# Function that checks if any of the metadata variables shows any bias 
checkBias <- function(metadataAvailable = TRUE, datasetSeurat,datasetName,PCAvalue,resolutionValue, metadataObject){
  # Check for any bias given by metadata variables above defined
  if (metadataAvailable == TRUE) {
    for (i in colnames(metadataObject)) { #i in metadataVariables
      print(i)
      if(!is.null(datasetSeurat@meta.data[[i]])){
        print(head(datasetSeurat@meta.data[[i]]))
        
        # Define new identity
        Idents(datasetSeurat) <- datasetSeurat@meta.data[[i]]
        
        # Generate plot
        drawDimPlot(dataset = datasetSeurat,datasetName = datasetName,pca = PCAvalue,resolution = resolutionValue,labelValue = T)
        Sys.sleep(0.5)
      }
    }}
}
