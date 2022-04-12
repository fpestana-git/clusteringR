# This function helps to determine which is the optimal PCA value to use
findOptimalPCA <- function(seuratObj){
  seuratObject <- seuratObj
  # Determine percent of variation associated with each PC
  pct <- seuratObject@reductions$pca@stdev / sum(seuratObject@reductions$pca@stdev) * 100
  # Calculate cumulative percents for each PC
  cum <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cum > 90 & pct < 5)[1]
  
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2) # change to any other number
  
  print(paste0("The best pca to use should be: ", pcs))
  pcs
}