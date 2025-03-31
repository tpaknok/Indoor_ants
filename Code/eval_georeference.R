# This script is for evaluating the possibility to use geo-referenced records only for the analyses

library(terra)

### using optimized records
ant_kass <- read.csv("Data/1_GABI_raw_all_records.csv") #need to download this from Kass et al. 2022
ant_kass_species <- read.csv("Data/3A_forAnalysis_processed_database_SPECIES.csv") #need to download this from Kass et al. 2022
native <- read.csv("Data/Native.csv")
clim_invasion_df <- read.csv("Data/clim_invasion_df.csv")
clim_invasion_df <- na.omit(clim_invasion_df)

ant_kass$label <- paste0(ant_kass$valid_species_name,"_",ant_kass$bentity2_name)
ant_kass$native <- ifelse(ant_kass$label %in% native$label,1,0)
ant_kass_non_native <- subset(ant_kass,native == 0)
ant_kass_native <- subset(ant_kass,native == 1)

###
bentity.shp <- vect("Data/Bentity2_shapefile_fullres.shp")
native$valid_species_name <- gsub(" ",".",native$valid_species_name)
ant_kass_subset <- ant_kass_species[ant_kass_species$valid_species_name %in% unique(clim_invasion_df$species),]

ant_kass_subset$bentity2_name <- ant_kass[match(ant_kass_subset$gabi_acc_number,ant_kass$gabi_acc_number),"bentity2_name"]
ant_kass_subset$label <- paste0(ant_kass_subset$valid_species_name,"_",ant_kass_subset$bentity2_name)
ant_kass_subset$native <- ifelse(ant_kass_subset$label %in% native$label,1,0)
ant_kass_subset <- subset(ant_kass_subset,native==1)

library(spThin)
ant_kass_subset_sp_df <- NULL
for (i in 1:length(unique(ant_kass_subset$valid_species_name))) {
  message(i)
  ant_kass_subset_sp <- subset(ant_kass_subset,valid_species_name == unique(ant_kass_subset$valid_species_name)[[i]])
  thinned <-　thin(ant_kass_subset_sp,
                  lat.col="lat_opt",
                  long.col="lon_opt",
                  spec.col="valid_species_name",
                  thin.par=10,
                  write.files=F,
                  locs.thinned.list.return=T,
                  reps=1)
  ant_kass_subset_sp_df <- rbind(ant_kass_subset_sp_df,
                                 ant_kass_subset_sp[as.numeric(rownames(thinned[[1]])),])
}

###
temp_pca <- rast("Raster/temp_pca.tiff") #Run 1_clim_raster.R to get this raster
water_pca <- rast("Raster/water_pca_untransformed.tiff") #Run 1_clim_raster.R to get this raster

climate_pca <- c(temp_pca,water_pca)
climate_georef <- extract(climate_pca,ant_kass_subset_sp_df[,c("lon_opt","lat_opt")])

ant_kass_subset_sp_df <- cbind(ant_kass_subset_sp_df,climate_georef)

###

library(tidyverse)

native_dist2 <- ant_kass_subset_sp_df %>%
  group_by(valid_species_name) %>%
  summarize(record = n(),
            n.regions=n_distinct(bentity2_name),
            PCA1_temp = mean(temp_pca,na.rm=T),
            PCA1_water = mean(water_pca,na.rm=T))
    
native_dist <- native %>% 
  filter(Final == "Yes") %>% 
  group_by(valid_species_name) %>%
  count(bentity2_name) %>%
  summarize(n.regions=sum(n)) %>%
  mutate(included =  valid_species_name %in% unique(clim_invasion_df$species)) %>%
  filter(included == T) %>%
  left_join(native_dist2,by=join_by(valid_species_name)) %>%
  mutate(n.regions.y=replace_na(n.regions.y,0)) %>%
  mutate(prop = n.regions.y/n.regions.x)

PCA_df <- clim_invasion_df %>%
  group_by(species) %>%
  summarize(PCA1_temp = mean(temp_pca,na.rm=T),
            PCA1_water = mean(water_pca,na.rm=T)) %>%
  left_join(native_dist[,c("valid_species_name","PCA1_temp","PCA1_water","record","prop")],by=join_by(species == valid_species_name))

PCA_df <- subset(PCA_df,record >= 50 & prop >= 0.8)

cor(PCA_df$PCA1_temp.x,PCA_df$PCA1_temp.y,method="kendall")
cor(PCA_df$PCA1_water.x,PCA_df$PCA1_water.y,method="kendall")

plot(PCA_df$PCA1_temp.x,PCA_df$PCA1_temp.y)

###
subset(native_dist,record >= 50 & prop >= 0.8)
mean(native_dist2$n.regions.x/native_dist2$n.regions.y)

###

clim_invasion_df <- clim_invasion_df %>%
  left_join(PCA_df,by=join_by(species))

clim_invasion_df$num <- ifelse(clim_invasion_df$num ==0,0,1)
library(glmmTMB)
m1 <- glmmTMB(num~temp_pca*mean_temp_pca_native+(temp_pca||species)+(mean_temp_pca_native||polygon_name),
              data=clim_invasion_df,family="binomial") # model does converge. but note the very low random effect variance of temp_pca.
summary(m1)

library(glmmTMB)
m2 <- glmmTMB(num~temp_pca*PCA1_temp.y+(temp_pca||species)+(PCA1_temp.y||polygon_name),
              data=clim_invasion_df,family="binomial") # model does converge. but note the very low random effect variance of temp_pca.
summary(m2)

####
ant_kass <- read.csv(file.choose())
ant_kass_clean <- read.csv(file.choose())
ant_kass_species <- read.csv(file.choose())

native <- subset(native,Final=="Yes")
ant_kass <- ant_kass[!is.na(ant_kass$dec_long) & !is.na(ant_kass$dec_lat),]
ant_kass_subset <- ant_kass[ant_kass$valid_species_name %in% unique(clim_invasion_df$species),]
ant_kass_subset$label <- paste0(ant_kass_subset$valid_species_name,"_",ant_kass_subset$bentity2_name)
ant_kass_subset$native <- ifelse(ant_kass_subset$label %in% native$label,1,0)

ant_kass_subset <- ant_kass_clean[match(ant_kass_subset$gabi_acc_number,ant_kass_clean$accession_number),]

native_dist2 <- ant_kass_subset %>%
  filter(native==1) %>%
  group_by(valid_species_name) %>%
  count(bentity2_name) %>%
  mutate(n=1) %>%
  summarize(n.regions=sum(n))

native_dist <- native %>% 
  filter(Final == "Yes") %>% 
  group_by(valid_species_name) %>%
  count(bentity2_name) %>%
  summarize(n.regions=sum(n)) %>%
  left_join(native_dist2,by=join_by(valid_species_name)) %>%
  mutate(included = valid_species_name %in% unique(clim_invasion_df$species)) %>%
  filter(included == T) %>%
  mutate(n.regions.y=replace_na(n.regions.y,0)) %>%
  mutate(prop = n.regions.y/n.regions.x)

mean(native_dist$prop)

###
ant_kass_clean <- ant_kass_clean %>%
  left_join(ant_kass[,c("valid_species_name","gabi_acc_number")],by=join_by(gabi_acc_number))

ant_kass_clean <- ant_kass_clean[!is.na(ant_kass_clean$lon_gabi) & !is.na(ant_kass_clean$lat_gabi),]
ant_kass_clean <- subset(ant_kass_clean,summary=="True")
ant_kass_subset <- ant_kass_clean[ant_kass_clean$valid_species_name %in% unique(clim_invasion_df$species),]
ant_kass_subset$label <- paste0(ant_kass_subset$valid_species_name,"_",ant_kass_subset$bentity2_name)
ant_kass_subset$native <- ifelse(ant_kass_subset$label %in% native$label,1,0)

native_dist2 <- ant_kass_subset %>%
  filter(native==1) %>%
  group_by(valid_species_name) %>%
  count(bentity2_name) %>%
  mutate(n=1) %>%
  summarize(n.regions=sum(n))

native_dist <- native %>% 
  filter(Final == "Yes") %>% 
  group_by(valid_species_name) %>%
  count(bentity2_name) %>%
  summarize(n.regions=sum(n)) %>%
  left_join(native_dist2,by=join_by(valid_species_name)) %>%
  mutate(included = valid_species_name %in% unique(clim_invasion_df$species)) %>%
  filter(included == T) %>%
  mutate(n.regions.y=replace_na(n.regions.y,0)) %>%
  mutate(prop = n.regions.y/n.regions.x)
