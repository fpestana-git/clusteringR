# Function that creates a spatial object based on sample Counts and sample Coordinates
createSpatialSeurat <- function(sampleCounts, sampleCoordinates){
  sampleCoordinates <- sampleCoordinates
  slide.seq <- CreateSeuratObject(counts = sampleCounts, assay="Spatial")
  
  coord.df <- data.frame(x=sampleCoordinates$X, y=sampleCoordinates$Y, stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
  rownames(coord.df) <- rownames(sampleCoordinates)
  
  slide.seq@images$image <-  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )
  # Update spatial metadata
  slide.seq@meta.data <- cbind(slide.seq@meta.data,sampleCoordinates)
  
  slide.seq
}
