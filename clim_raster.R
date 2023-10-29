library(terra)
library(spdep)
setwd("C:/Users/pakno/OneDrive - University of Toronto/Raster/TerraClim")
prec <- rast("TerraClimate19812010_ppt.nc")
sm <- rast("TerraClimate19812010_soil.nc")

#prec_min <- sqrt(min(prec))
#prec_max <- sqrt(max(prec))
#prec_sum <- sqrt(sum(prec))
prec_min <- min(prec)
prec_max <- max(prec)
sm_min <- min(sm)
sm_max <- max(sm)

prec_sd <- stdev(prec)
sm_sd <- stdev(sm)
#sm_min <- sqrt(min(sm))
#sm_max <- sqrt(max(sm))
#sm_mean <- sqrt(mean(sm))
#prec_sd <- sqrt(stdev(prec))
#sm_sd <- sqrt(stdev(sm))
tmin <- rast("TerraClimate19812010_tmin.nc")
tmax <- rast("TerraClimate19812010_tmax.nc")
tmonth <- (tmin+tmax)/2
AT_mean_5km <- mean(tmonth)
AT_min_5km <- min(tmin)
AT_max_5km <- max(tmax)

#AT_min_sd_5km <- stdev(tmin)
#AT_max_sd_5km <- stdev(tmax)
AT_mean_sd_5km <- stdev(tmonth)

#clim <- c(AT_min_5km,AT_max_5km,AT_min_sd_5km,AT_max_sd_5km,prec_min,prec_max,prec_sd,sm_min,sm_max,sm_sd)
#names(clim) <- c("tmin","tmax","tmin_sd","tmax_sd","Prec_min","Prec_max","Prec_sd","SM_min","SM_max","SM_sd")
#temp_current<- c(AT_min_5km,AT_max_5km,AT_min_sd_5km,AT_max_sd_5km)
#water_current <- c(prec_min,prec_max,prec_sd,sm_min,sm_max,sm_sd)
#temp_pca <- prcomp(as.data.frame(temp_current),center=T,scale.=T)
#water_pca <- prcomp(as.data.frame(water_current),center=T,scale.=T)
#BiodiversityR::PCAsignificance(vegan::rda(scale(as.data.frame(temp_current))~1),axes=4)
#BiodiversityR::PCAsignificance(vegan::rda(scale(as.data.frame(water_current))~1),axes=6)

clim <- c(AT_min_5km,AT_max_5km,AT_mean_sd_5km,prec_min,prec_max,prec_sd,sm_min,sm_max,sm_sd)
names(clim) <- c("tmin","tmax","tmean_sd","Prec_min","Prec_max","Prec_sd","SM_min","SM_max","SM_sd")
cor(na.omit(as.data.frame(clim)))

temp_current<- c(AT_min_5km,AT_max_5km,AT_mean_sd_5km)
water_current <- c(prec_min,prec_max,prec_sd,sm_min,sm_max,sm_sd)
temp_pca <- prcomp(as.data.frame(temp_current),center=T,scale.=T)
water_pca <- prcomp(na.omit(as.data.frame(water_current)),center=T,scale.=T)
BiodiversityR::PCAsignificance(vegan::rda(scale(as.data.frame(temp_current))~1),axes=2)
BiodiversityR::PCAsignificance(vegan::rda(scale(na.omit(as.data.frame(water_current)))~1),axes=4)


AT_min_future <- rast("TerraClimate2C_tmin.nc")
AT_max_future <- rast("TerraClimate2C_tmax.nc")
prec_future <- rast("TerraClimate2C_ppt.nc")
sm_future <- rast("TerraClimate2C_soil.nc")

prec_min_f <- log(min(prec_future+1))
prec_max_f <- log(max(prec_future+1))
#prec_sum_f <- sum(prec_future))
sm_min_f <- log(min(sm_future+1))
sm_max_f <- log(max(sm_future+1))
#sm_mean_f <- mean(sm_future))
prec_sd_f <- log(stdev(prec_future)+1)
sm_sd_f <- log(stdev(sm_future)+1)

#prec_min_f <- sqrt(min(prec_future))
#prec_max_f <- sqrt(max(prec_future))
#prec_sum_f <- sqrt(sum(prec_future))
#sm_min_f <- sqrt(min(sm_future))
#sm_max_f <- sqrt(max(sm_future))
#sm_mean_f <- sqrt(mean(sm_future))
#prec_sd_f <- sqrt(stdev(prec_future))
#sm_sd_f <- sqrt(stdev(sm_future))
tmonth_future <- (AT_min_future+AT_max_future)/2
#AT_mean_5km_f <- mean(tmonth_future)
AT_min_5km_f <- min(AT_min_future)
AT_max_5km_f <- max(AT_max_future)
#AT_min_sd_5km_f <- stdev(AT_min_future)
#AT_max_sd_5km_f <- stdev(AT_max_future)
AT_mean_sd_5km_f <- stdev(tmonth_future)

AT_min_future_4C <- rast("TerraClimate4C_tmin.nc")
AT_max_future_4C <- rast("TerraClimate4C_tmax.nc")
prec_future_4C <- rast("TerraClimate4C_ppt.nc")
sm_future_4C <- rast("TerraClimate4C_soil.nc")

prec_min_f_4C <- log(min(prec_future_4C+1))
prec_max_f_4C <- log(max(prec_future_4C+1))
#prec_sum_f_4C <- sum(prec_future_4C))
sm_min_f_4C <- log(min(sm_future_4C+1))
sm_max_f_4C <- log(max(sm_future_4C+1))
#sm_mean_f_4C <- mean(sm_future_4C))
prec_sd_f_4C <- log(stdev(prec_future_4C)+1)
sm_sd_f_4C <- log(stdev(sm_future_4C)+1)

#prec_min_f_4C <- sqrt(min(prec_future_4C))
#prec_max_f_4C <- sqrt(max(prec_future_4C))
#prec_sum_f_4C <- sqrt(sum(prec_future_4C))
#sm_min_f_4C <- sqrt(min(sm_future_4C))
#sm_max_f_4C <- sqrt(max(sm_future_4C))
#sm_mean_f_4C <- sqrt(mean(sm_future_4C))
#prec_sd_f_4C <- sqrt(stdev(prec_future_4C))
#sm_sd_f_4C <- sqrt(stdev(sm_future_4C))

tmonth_future_4C <- (AT_min_future_4C+AT_max_future_4C)/2
#AT_mean_5km_f_4C <- mean(tmonth_future_4C)
AT_min_5km_f_4C <- min(AT_min_future_4C)
AT_max_5km_f_4C <- max(AT_max_future_4C)
#AT_min_sd_5km_f_4C <- stdev(AT_min_future_4C)
#AT_max_sd_5km_f_4C <- stdev(AT_max_future_4C)
AT_mean_sd_5km_f_4C <- stdev(tmonth_future_4C)

#clim_f <- c(AT_min_5km_f,AT_max_5km_f,AT_min_sd_5km_f,AT_max_sd_5km_f,
            #prec_min_f,prec_max_f,prec_sd_f,sm_min_f,sm_max_f,sm_sd_f,
            #AT_min_5km_f_4C,AT_max_5km_f_4C,AT_min_sd_5km_f_4C,AT_max_sd_5km_f_4C,
            #prec_min_f_4C,prec_max_f_4C,prec_sd_f_4C,sm_min_f_4C,sm_max_f_4C,sm_sd_f_4C)

#names(clim_f) <- c("tmin.2C","tmax.2C","tmin_sd.2C","tmax_sd.2C","Prec_min.2C","Prec_max.2C","Prec_sd.2C","sM_min.2C","SM_max.2C","SM_sd.2C",
                   #"tmin.4C","tmax.4C","tmin_sd.4C","tmax_sd.4C","Prec_min.4C","Prec_max.4C","Prec_sd.4C","sM_min.4C","SM_max.4C","SM_sd.4C")

temp.2C <- c(AT_min_5km_f,AT_max_5km_f,AT_mean_sd_5km_f)
temp.4C <- c(AT_min_5km_f_4C,AT_max_5km_f_4C,AT_mean_sd_5km_f_4C)
water.2C <- c(prec_min_f,prec_max_f,prec_sd_f,sm_min_f,sm_max_f,sm_sd_f)
water.4C <- c(prec_min_f_4C,prec_max_f_4C,prec_sd_f_4C,sm_min_f_4C,sm_max_f_4C,sm_sd_f_4C)

#temp.2C <- c(AT_mean_5km_f,AT_mean_sd_5km_f)
#temp.4C <- c(AT_mean_5km_f_4C,AT_mean_sd_5km_f_4C)
#water.2C <- c(prec_sum_f,prec_sd_f,sm_mean_f,sm_sd_f)
#water.4C <- c(prec_sum_f_4C,prec_sd_f_4C,sm_mean_f_4C,sm_sd_f_4C)

###############################################################
gc()
final_df <- cbind(as.data.frame(AT_min_5km_f,xy=T)[,1:2],
           temp_pca=predict(temp_pca)[,1],
           water_pca=predict(water_pca)[,1],
           temp_pca_2C = predict(temp_pca,as.data.frame(temp.2C))[,1],
           temp_pca_4C = predict(temp_pca,as.data.frame(temp.4C))[,1],
           water_pca_2C = predict(water_pca,as.data.frame(water.2C))[,1],
           water_pca_4C =predict(water_pca,as.data.frame(water.4C))[,1])


final_raster <- rast(final_df)

plot(final_raster)
writeRaster(final_raster[[1]],"temp_pca.tiff",overwrite=T)
writeRaster(final_raster[[2]],"water_pca_log.tiff",overwrite=T)
writeRaster(final_raster[[3]],"temp_pca_2C.tiff",overwrite=T)
writeRaster(final_raster[[4]],"temp_pca_4C.tiff",overwrite=T)
writeRaster(final_raster[[5]],"water_pca_2C_log.tiff",overwrite=T)
writeRaster(final_raster[[6]],"water_pca_4C_log.tiff",overwrite=T)

#################################################################
library(terra)
soil_min_surface <- rast("SBIO6_0_5cm_MinT_coldestMonth.tif")
soil_min_soil <- rast("SBIO6_5_15cm_MinT_coldestMonth.tif")
soil_max_surface <- rast("SBIO5_0_5cm_MaxT_warmestMonth.tif")
soil_max_soil <- rast("SBIO5_5_15cm_MaxT_warmestMonth.tif")

soil_min_surface <- project(soil_min_surface,tmin)
soil_min_surface <- min(soil_min_surface)
soil_min_soil <- project(soil_min_soil,tmin)
soil_min_soil <- min(soil_min_soil)

soil_max_surface <- project(soil_max_surface,tmin)
soil_max_surface <- max(soil_max_surface)

soil_max_soil <- project(soil_max_soil,tmin)
soil_max_soil <- max(soil_max_soil)

min_temp_corr <- as.data.frame(c(soil_min_surface,soil_min_soil,AT_min_5km))
max_temp_corr <- as.data.frame(c(soil_max_surface,soil_max_soil,AT_max_5km))
cor(na.omit(min_temp_corr))
cor(na.omit(max_temp_corr))
