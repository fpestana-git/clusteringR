#' Seurat spatial object
#' 
#' Function that takes a spatial dataset and respective x and y coordinates to generate the spatial map of an individual sample
#' 
#' @param datasetObject dataset object
#' @param metadataObject metadata object that includes X and Y coordinates for each cell
#'
#' @export createSpatialSeurat


# This function creates a spatial object based on sample Counts and sample Coordinates
createSpatialSeurat <- function(datasetObject, metadataObject){
  #sampleCoordinates <- sampleCoordinates
  slide.seq <- CreateSeuratObject(counts = datasetObject, assay="Spatial")
  
  coord.df <- data.frame(x=metadataObject$X, y=metadataObject$Y, stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
  rownames(coord.df) <- rownames(metadataObject)
  
  slide.seq@images$image <-  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )
  # Update spatial metadata
  slide.seq@meta.data <- cbind(slide.seq@meta.data,metadataObject)
  
  slide.seq
}
