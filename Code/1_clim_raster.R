library(terra)
library(spdep)
prec <- rast("Raster/ClimRast/TerraClimate19812010_ppt.nc") #download the raster yourself before running the script
sm <- rast("Raster/ClimRast/TerraClimate19812010_soil.nc")

prec_min <- min(prec) #min monthly precipitation 
prec_max <- max(prec) #max
sm_min <- min(sm) #min monthly soil moisture
sm_max <- max(sm) #max 

prec_sd <- stdev(prec) #precipitation sd along a year
sm_sd <- stdev(sm) #soil moisture sd along a year
tmin <- rast("TerraClimate19812010_tmin.nc")
tmax <- rast("TerraClimate19812010_tmax.nc")
tmonth <- (tmin+tmax)/2 #as a proxy of average monthly temp
AT_mean_5km <- mean(tmonth) #mean annual temperature
AT_min_5km <- min(tmin) #minimum monthly temp
AT_max_5km <- max(tmax) #max
AT_mean_sd_5km <- stdev(tmonth) #sd (temperature seasonality)

clim <- c(AT_min_5km,AT_max_5km,AT_mean_sd_5km,prec_min,prec_max,prec_sd,sm_min,sm_max,sm_sd)
names(clim) <- c("tmin","tmax","tmean_sd","Prec_min","Prec_max","Prec_sd","SM_min","SM_max","SM_sd")
cor(na.omit(as.data.frame(clim)))

temp_current<- c(AT_min_5km,AT_max_5km,AT_mean_sd_5km)
water_current <- c(prec_min,prec_max,prec_sd,sm_min,sm_max,sm_sd)
temp_pca <- prcomp(as.data.frame(temp_current),center=T,scale.=T) #pca on temp
water_pca <- prcomp(na.omit(as.data.frame(water_current)),center=T,scale.=T) #pca on precipitation
BiodiversityR::PCAsignificance(vegan::rda(scale(as.data.frame(temp_current))~1),axes=2)
BiodiversityR::PCAsignificance(vegan::rda(scale(na.omit(as.data.frame(water_current)))~1),axes=4)

###Did the same for 2C and 4C scenario
AT_min_future <- rast("TerraClimate2C_tmin.nc")
AT_max_future <- rast("TerraClimate2C_tmax.nc")
prec_future <- rast("TerraClimate2C_ppt.nc")
sm_future <- rast("TerraClimate2C_soil.nc")

prec_min_f <- min(prec_future)
prec_max_f <- max(prec_future)
sm_min_f <- min(sm_future)
sm_max_f <- max(sm_future)
prec_sd_f <- stdev(prec_future)
sm_sd_f <- stdev(sm_future)
tmonth_future <- (AT_min_future+AT_max_future)/2
AT_min_5km_f <- min(AT_min_future)
AT_max_5km_f <- max(AT_max_future)
AT_mean_sd_5km_f <- stdev(tmonth_future)

AT_min_future_4C <- rast("TerraClimate4C_tmin.nc")
AT_max_future_4C <- rast("TerraClimate4C_tmax.nc")
prec_future_4C <- rast("TerraClimate4C_ppt.nc")
sm_future_4C <- rast("TerraClimate4C_soil.nc")

prec_min_f_4C <- min(prec_future_4C)
prec_max_f_4C <- max(prec_future_4C)
sm_min_f_4C <- min(sm_future_4C)
sm_max_f_4C <- max(sm_future_4C)
prec_sd_f_4C <- stdev(prec_future_4C)
sm_sd_f_4C <- stdev(sm_future_4C)

tmonth_future_4C <- (AT_min_future_4C+AT_max_future_4C)/2
AT_min_5km_f_4C <- min(AT_min_future_4C)
AT_max_5km_f_4C <- max(AT_max_future_4C)
AT_mean_sd_5km_f_4C <- stdev(tmonth_future_4C)

temp.2C <- c(AT_min_5km_f,AT_max_5km_f,AT_mean_sd_5km_f)
temp.4C <- c(AT_min_5km_f_4C,AT_max_5km_f_4C,AT_mean_sd_5km_f_4C)
water.2C <- c(prec_min_f,prec_max_f,prec_sd_f,sm_min_f,sm_max_f,sm_sd_f)
water.4C <- c(prec_min_f_4C,prec_max_f_4C,prec_sd_f_4C,sm_min_f_4C,sm_max_f_4C,sm_sd_f_4C)

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
writeRaster(final_raster[[1]],"Raster/temp_pca.tiff",overwrite=T)
writeRaster(final_raster[[2]],"Raster/water_pca.tiff",overwrite=T)
writeRaster(final_raster[[3]],"Raster/temp_pca_2C.tiff",overwrite=T)
writeRaster(final_raster[[4]],"Raster/temp_pca_4C.tiff",overwrite=T)
writeRaster(final_raster[[5]],"Raster/water_pca_2C.tiff",overwrite=T)
writeRaster(final_raster[[6]],"Raster/water_pca_4C.tiff",overwrite=T)

###check if soil temperature is correlated with air temperature, using data from the database SoilTemp
library(terra)
soil_min_surface <- rast("Raster/ClimRast/SBIO6_0_5cm_MinT_coldestMonth.tif")
soil_min_soil <- rast("Raster/ClimRast/SBIO6_5_15cm_MinT_coldestMonth.tif")
soil_max_surface <- rast("Raster/ClimRast/SBIO5_0_5cm_MaxT_warmestMonth.tif")
soil_max_soil <- rast("Raster/ClimRast/SBIO5_5_15cm_MaxT_warmestMonth.tif")

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

### check climatic conditinos under all scenarios for table 
temp <- c(AT_min_5km,AT_min_5km_f,AT_min_5km_f_4C,AT_max_5km,AT_max_5km_f,AT_max_5km_f_4C)

set <- terra::extract(temp,bentity.shp,fun=mean,na.rm=T) #get mean climatic conditinos for each polygon
set <- na.omit(set)
sd(set[,3]-set[,2])
sd(set[,4]-set[,2])
mean(set[,3]-set[,2])
mean(set[,4]-set[,2])

sd(set[,6]-set[,5])
sd(set[,7]-set[,5])
mean(set[,6]-set[,5])
mean(set[,7]-set[,5])

diff_min_2 <-AT_min_5km_f-AT_min_5km
diff_min_4 <- AT_min_5km_f_4C-AT_min_5km
diff_max_2 <- AT_max_5km_f-AT_max_5km
diff_max_4 <- AT_max_5km_f_4C-AT_max_5km

global(diff_min_2,"sd",na.rm=T)
global(diff_min_4,"sd",na.rm=T)
global(diff_max_2,"sd",na.rm=T)
global(diff_max_4,"sd",na.rm=T)
