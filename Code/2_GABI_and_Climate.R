require(terra)
require(tidyverse)

#these files are too large for upload. Please download the terraclimate data and use the clim_raster.R to make them
temp_pca <- rast("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/Raster/temp_pca.tiff")
water_pca <- rast("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/Raster/water_pca_untransformed.tiff")
temp_pca_2C <- rast("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/Raster/temp_pca_2C.tiff")
temp_pca_4C <- rast("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/Raster/temp_pca_4C.tiff")
water_pca_2C <- rast("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/Raster/water_pca_2C_untransformed.tiff")
water_pca_4C <- rast("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/Raster/water_pca_4C_untransformed.tiff")

clim <- c(temp_pca,water_pca,temp_pca_2C,temp_pca_4C,water_pca_2C,water_pca_4C)

native.records <- read.csv("Data/Native.csv")
native.records$valid_species_name <- gsub(" ",".",native.records$valid_species_name)
exotic.records <- read.csv("Data/Exotic_record_date.csv")

#############fix spelling
exotic.records$valid_species_name[exotic.records$valid_species_name == "Pheidole.indica "] <- "Pheidole.indica"
exotic.records$valid_species_name[exotic.records$valid_species_name == "Hyponera.punctatissima"] <- "Hypoponera.punctatissima"

#############
exotic.records$valid_species_name <- gsub(" ",".",exotic.records$valid_species_name)

bentity.match <- read.csv("Data/match_bentity.csv")
bentity.match <- bentity.match[bentity.match$Bentity.assigned != "Didn't assign",] #reassign some records to different poylgons.

nrow(subset(exotic.records,Improved.Classification != "Native"))
name.to.be.replaced <- exotic.records$bentity2_name[exotic.records$bentity2_name %in% bentity.match$Original.name]
exotic.records$bentity2_name[exotic.records$bentity2_name %in% bentity.match$Original.name] <- bentity.match[match(name.to.be.replaced,bentity.match$Original.name),"Bentity.assigned"]

exotic.records <- subset(exotic.records,Improved.Classification != "Intercepted" & Improved.Classification != "Indoor Introduced (UNCERTAIN)" & Improved.Classification != "Native" & bentity2_name != "Didn't assign")
exotic.records[exotic.records$Improved.Classification == "exotic","Improved.Classification"] <- "Exotic"
exotic.records$pop <- paste0(exotic.records$valid_species_name," ",exotic.records$bentity2_name)
unique(exotic.records$pop)

### check if the earliest indoor records precede the earliest outdoor records
set <- exotic.records %>% 
  filter(Improved.Classification == "Exotic" | Improved.Classification == "Indoor Introduced") %>% 
  group_by(valid_species_name,bentity2_name,Improved.Classification) %>% 
  arrange(Date) %>% 
  slice(1L) %>%
  ungroup() %>%
  group_by(pop) %>%
  summarize(diff = Date[Improved.Classification == "Indoor Introduced"] - Date[Improved.Classification == "Exotic"])

### The number of indoor and outdoor records for each species
status_count <- exotic.records %>% group_by(valid_species_name) %>% dplyr::count(Improved.Classification)
status_count_wide <- status_count %>% pivot_wider(id_cols=valid_species_name,names_from=Improved.Classification,values_from=n)
status_count_wide$n <- status_count_wide$`Indoor Introduced`+status_count_wide$Exotic

###data exploration + preparation for extracting climatic conditions
classification <- read.csv("Data/Classification_Mark.csv") #Mark data
classification$Species.name <- gsub(" ",".",classification$Species.name)
bentity.shp <- vect("Data/Bentity2_shapefile_fullres.shp") #Bentity polygon
bentity.df <- as.data.frame(bentity.shp)

na.native <- native.records[which(is.na(match(native.records$bentity2_name, bentity.df$BENTITY2_N))),] #all native records can be associated with the bentity system
na.exotic <- exotic.records[which(is.na(match(exotic.records$bentity2_name, bentity.df$BENTITY2_N))),] #96 records cannot be associated with the bentity system 

subset.native <- native.records[which(!is.na(match(native.records$bentity2_name, bentity.df$BENTITY2_N))),] #remove those that cannot be associated
subset.exotic <- exotic.records[which(!is.na(match(exotic.records$bentity2_name, bentity.df$BENTITY2_N))),] # same
subset.exotic$shp <- match(subset.exotic$bentity2_name,bentity.df$BENTITY2_N)

sp.list <- unique(c(subset.native$valid_species_name,subset.exotic$valid_species_name)) #alien species with known native distribution
no.native <- sp.list[sp.list %in% unique(subset.exotic$valid_species_name) & !sp.list %in% unique(subset.native$valid_species_name)] #14 species have no native distribution data

sp.list <- sp.list[sp.list %in% unique(subset.exotic$valid_species_name) & sp.list %in% unique(subset.native$valid_species_name)]
sp.list <- sp.list[sp.list %in% unique(status_count_wide$valid_species_name)]

subset.exotic$populations <- paste0(subset.exotic$bentity2_name,subset.exotic$valid_species_name)

no_year <- subset.exotic %>% group_by(bentity2_name,valid_species_name) %>% summarise(mean_year = mean(Date,na.rm=T)) #see if some populations have no dated records
length(which(is.nan(no_year$mean_year)))

###extract climatic conditions (native + exotic) for each population
overall.native <- NULL
overall.exotic <- NULL
na.species <- NULL

set <- terra::extract(clim,bentity.shp,fun=mean,na.rm=T) #get mean climatic conditinos for each polygon
set <- cbind(set,bentity2_name=bentity.shp$BENTITY2_N)
set$area <- expanse(bentity.shp) #extract area - useful for faster calculations later (instead of averaging across all grids, first calculat menas of each polygon, then weighted avg based on area of each poylgon)

colnames(set) <- c("ID","temp_pca","water_pca","temp_pca_2C","temp_pca_4C","water_pca_2C","water_pca_4C",
                  "bentity2_name","area") #make sure name is right

subset.exotic$num <- ifelse(subset.exotic$Improved.Classification == "Indoor Introduced",0,1)
clim_invasion_df <- list()

centroid_df <- as.data.frame(centroids(bentity.shp),geom="XY") #get coordinates
subset.exotic <- cbind(subset.exotic,centroid_df[match(subset.exotic$bentity2_name,centroid_df$BENTITY2_N),c("x","y")])

### Huge for-loop for climatic conditinos of each species 
### Note that instead of extracting climatic data from each grid for each species, we simply used polygon-level climate (based on averaging across grids) in the loop.
### These data have been extracted in L80-82
### This makes the script much more computationally effect (extracting grid-level data for 336 speceis is very slow!!)

for (sp in sp.list) {
  message("Analyzing ",sp,"; ",which(sp.list %in% sp),"/",length(sp.list))
  layer <- classification[classification$Species.name == sp,]
  sp.records.native <- subset(subset.native, valid_species_name==sp)
  sp.records.exotic <- subset(subset.exotic, valid_species_name==sp)
  
  sp.records.native <- sp.records.native[!duplicated(sp.records.native$bentity2_name),]
  sp.records.exotic <- sp.records.exotic %>% group_by(bentity2_name)  %>% arrange(Date) %>% mutate(num1=mean(num)) %>% slice(1L)
  
  native.clim.df <- set[match(sp.records.native$bentity2_name,set$bentity2_name),]
  exotic.clim.df <- set[match(sp.records.exotic$bentity2_name,set$bentity2_name),]
  
  native.climatic.niche <- native.clim.df %>% summarise(
                                                     mean_temp_pca_native = weighted.mean(temp_pca,w=area,na.rm=T),
                                                     mean_water_pca_native = weighted.mean(water_pca,w=area,na.rm=T))
  
  clim_invasion_df[[sp]] <- data.frame(genus=word(sp,1,1,"\\."),species=sp,polygon_name=unique(sp.records.exotic$bentity2_name),
                                       x= unique(sp.records.exotic$x),
                                       y = unique(sp.records.exotic$y),
                                       num=sp.records.exotic$num1,exotic.clim.df,date=sp.records.exotic$Date,
                                       sp.layer=layer$strata_classified,native.climatic.niche)
}

clim_invasion_df <- do.call(rbind,clim_invasion_df)
clim_invasion_df$population <- paste0(clim_invasion_df$bentity2_name,"_",clim_invasion_df$species)
nrow(na.omit(clim_invasion_df)) #number of populations analyzed (after omitting any records with NA)

clim_invasion_df$id <- paste0(clim_invasion_df$species,"_",clim_invasion_df$polygon_name)

exotic.records$record <- paste0(exotic.records$valid_species_name,"_",exotic.records$bentity2_name)
exotic.records$num <- ifelse(exotic.records$Improved.Classification == "Exotic",1,0) #numeric variable for indoor status (for later analyses)

oldest_record <- aggregate(Date~record,data=exotic.records,FUN=min)
classification_status <- aggregate(num~record,data=exotic.records,FUN=mean) #calculate mean allows us to see if multiple records are consistent in terms of indoor status

clim_invasion_df$date <- oldest_record[match(clim_invasion_df$id,oldest_record$record),"Date"]
clim_invasion_df$num <- classification_status[match(clim_invasion_df$id,classification_status$record),"num"]

clim_invasion_df$sp.layer <-  classification$strata_classified[match(clim_invasion_df$species,classification$Species.name)]
centroid <- terra::as.data.frame(centroids(bentity.shp),geom="XY")
clim_invasion_df <- cbind(clim_invasion_df,centroid[match(clim_invasion_df$polygon_name,centroid$BENTITY2_N),c("x","y")])
write.csv(clim_invasion_df,"Data/clim_invasion_df.csv")
