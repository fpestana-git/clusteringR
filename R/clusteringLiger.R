clusteringLiger <- function(datasets, lambdaValue = 3, kValue = 25, referenceDatasetName, resolutionValue = 0.4){
  liger <- createLiger(raw.data = datasets,remove.missing = T)# Normalize datasets
  # Normalize the datasets
  liger <- rliger::normalize(liger)
  # Try selecting variable genes automatically
  liger@var.genes = variableGenes
  liger <- scaleNotCenter(liger, remove.missing = TRUE,verbose = T)
  # run with different settings
  liger <- optimizeALS(liger,  lambda = lambdaValue , use.unshared = F, max.iters = 30, thresh=1e-6, k =kValue, vectorized.lambda = TRUE)
  liger <- quantile_norm(liger, ref_dataset = referenceDatasetName)
  liger <- louvainCluster(liger,resolution = resolutionValue)
  liger <- runUMAP(liger, n_neighbors = 35) # do I need to change the n_neighbors?
  
  liger
}