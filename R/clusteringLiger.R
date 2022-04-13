#' LIGER clustering using iNMF or UINMF algorythm
#' 
#' Function that takes a named list of datasets and integrates using LIGER
#' 
#' @param datasets list of named datasets
#' @param lambdaValue default 5
#' @param kValue default 25
#' @param referenceDatasetName name of dataset to use as reference dataset
#' @param resolutionValue default 0.4
#' @param setVariableGenes default FALSE
#' @param num.genesValue default NULL
#' @param var.threshValue default NULL
#' @param useiNMF default FALSE
#' @param useUINMF default FALSE
#' @param unshared.threshValue default NULL
#' 
#' @export clusteringLiger

clusteringLiger <- function(datasets, lambdaValue = 5, kValue = 25, referenceDatasetName, resolutionValue = 0.4,setVariableGenes = FALSE,num.genesValue = NULL, var.threshValue = NULL,useiNMF = FALSE, useUINMF = FALSE,unshared.threshValue = NULL){

    if(useUINMF){
      liger <- createLiger(raw.data = datasets)# Normalize datasets
      liger <- rliger::normalize(liger)
      liger@var.genes = variableGenes
      liger <- selectGenes(object = liger,unshared = T,unshared.datasets = list(3),unshared.thresh = unshared.threshValue) # only for new liger
      liger <- scaleNotCenter(object = liger,verbose = T)
      liger <- optimizeALS(liger,  lambda = lambdaValue , use.unshared = T, max.iters = 30, thresh=1e-6, k =kValue)
    }
    
    if(useiNMF){
      liger <- createLiger(raw.data = datasets,remove.missing = T)# Normalize datasets # for normal liger
      liger <- rliger::normalize(liger)
      liger@var.genes = variableGenes
      liger <- scaleNotCenter(liger, remove.missing = TRUE,verbose = T) # for normal liger
      liger <- optimizeALS(liger,  lambda = lambdaValue , use.unshared = F, max.iters = 50, thresh=1e-6, k =kValue, vectorized.lambda = TRUE) # for normal liger
      
    }
    
  # run with different settings
  liger <- quantile_norm(liger, ref_dataset = referenceDatasetName)
  liger <- louvainCluster(liger,resolution = resolutionValue)
  liger <- runUMAP(liger, n_neighbors = 30) # do I need to change the n_neighbors?
  
  liger
}