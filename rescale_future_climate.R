rescale_future_climate <- function(future_climate,current_climate,pca_climate_mod,current_pca_score) {
  current_mean <- colMeans(current_climate)
  sd <- apply(current_climate,2,sd)
  
  future_climate_rescaled <- t(apply(future_climate, 1,function(x) (x-current_mean)/sd))
  future_climate_rescaled <- apply(future_climate_rescaled,1, function(x) rowSums(x*t(pca_climate_mod$rotation[,1])))
  future_climate_rescaled <- (future_climate_rescaled - mean(current_pca_score))/sd(current_pca_score)
  
  return(future_climate_rescaled)
}

