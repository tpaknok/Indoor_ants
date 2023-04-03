setwd("C:/Users/pakno/OneDrive - University of Toronto")
source("single_var.R")
source("make_genus_tree.R")

library(phytools)
library(ape)
library(stringr)
library(tidyverse)
library(terra)
NMI_temp <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/NMI_temp.csv")
clim_invasion_df <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/clim_invasion_df_14022023.csv")
clim_invasion_df_f <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/clim_invasion_df_f_14022023.csv")
clim_invasion_df <- cbind(clim_invasion_df,clim_invasion_df_f[,6:17])
clim_invasion_df$num <- ifelse(clim_invasion_df$num ==0,0,1)

ant_tree <- read.tree("C:/Users/pakno/OneDrive - University of Toronto/ant_tree/backbone_MLtree_RaxML.tre")
is.ultrametric(ant_tree)
ant_tree <- chronos(ant_tree,lambda=0)
ant_tree$node.label <- 1:length(ant_tree$node.label)

ant_tree$tip.label <- gsub("[\"]","",ant_tree$tip.label)
ant_tree

ant_placement <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/ant_genus_placement.csv")
#genus_tree <- make_genus_tree(ant_tree,ant_placement,unique(NMI_temp$sp))
genus_tree <- make_genus_tree(ant_tree,ant_placement,unique(clim_invasion_df$species))

NMI_analysis <- subset(clim_invasion_df)
NMI_analysis$date <- scale(log10(NMI_analysis$date))

NMI_analysis$species <- gsub("[.]","_",NMI_analysis$species)
NMI_analysis$species <- gsub(" ","_",NMI_analysis$species)
NMI_analysis <- na.omit(NMI_analysis)
NMI_analysis_pca_temp <- prcomp(scale(NMI_analysis[,c("MAT","Temp_sd")]))
NMI_analysis$temp_pca <- apply(scale(NMI_analysis[,c("MAT","Temp_sd")]),1, function(x) rowSums(x*t(NMI_analysis_pca_temp$rotation[,1])))
NMI_analysis_pca_water <- prcomp(scale(NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")]))
NMI_analysis$water_pca <- apply(scale(NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")]),1, function(x) rowSums(x*t(NMI_analysis_pca_water$rotation[,1])))
NMI_analysis_niche_pca_temp <- prcomp(scale(NMI_analysis[,c("niche.MAT","niche.Temp_sd")]))
NMI_analysis$niche_temp_pca <- apply(scale(NMI_analysis[,c("niche.MAT","niche.Temp_sd")]),1, function(x) rowSums(x*t(NMI_analysis_niche_pca_temp$rotation[,1])))
NMI_analysis_niche_pca_water <- prcomp(scale(NMI_analysis[,c("niche.Prec","niche.Prec_sd","niche.SM","niche.SM_sd")]))
NMI_analysis$niche_water_pca <- apply(scale(NMI_analysis[,c("niche.Prec","niche.Prec_sd","niche.SM","niche.SM_sd")]),1, function(x) rowSums(x*t(NMI_analysis_niche_pca_water$rotation[,1])))

NMI_analysis <- cbind(NMI_analysis,scaled=scale(NMI_analysis[,c("water_pca","niche_water_pca","temp_pca","niche_temp_pca","date")]))

NMI_analysis$num <- ifelse(NMI_analysis$num >0, 1,0)

#####################Characterize environments
options(glmmTMB.cores=6)
NMI_analysis$factor <- as.factor(NMI_analysis$num)
library(glmmTMB)

m_MAT <- lm((NMI_analysis$MAT-NMI_analysis$niche.MAT)~factor,data=NMI_analysis)
summary(m_MAT)
car::Anova(m_MAT,white.adjust=T)

m_Temp_sd <-lm((NMI_analysis$Temp_sd-NMI_analysis$niche.Temp_sd)~factor,data=NMI_analysis)
summary(m_Temp_sd)
car::Anova(m_Temp_sd,white.adjust=T)

m_Prec <- lm((NMI_analysis$Prec-NMI_analysis$niche.Prec)~factor,data=NMI_analysis)
summary(m_Prec)
car::Anova(m_Prec,white.adjust=T)

m_Prec_sd <-lm((NMI_analysis$Prec_sd-NMI_analysis$niche.Prec_sd)~factor,data=NMI_analysis)
summary(m_Prec_sd)
car::Anova(m_Prec_sd,white.adjust=T)

m_SM <- lm((NMI_analysis$SM-NMI_analysis$niche.SM)~factor,data=NMI_analysis)
summary(m_SM)
car::Anova(m_SM,white.adjust=T)

m_SM_sd <-lm((NMI_analysis$SM_sd-NMI_analysis$niche.SM_sd)~factor,data=NMI_analysis)
summary(m_SM_sd)
car::Anova(m_SM_sd,white.adjust=T)

library(glmmTMB)
NMI_analysis$Soil <- 0
NMI_analysis$Soil[NMI_analysis$sp.layer == "Soil"] <- 1
NMI_analysis$Generalist <- 0
NMI_analysis$Generalist[NMI_analysis$sp.layer == "Generalist"] <- 1
NMI_analysis$Arboreal<- 0
NMI_analysis$Arboreal[NMI_analysis$sp.layer == "Arboreal"] <- 1
NMI_analysis$Soil <- scale(NMI_analysis$Soil,center=T,scale=F)
NMI_analysis$Generalist <- scale(NMI_analysis$Generalist,center=T,scale=F)
NMI_analysis$Arboreal <- scale(NMI_analysis$Arboreal,center=T,scale=F)

m <- glmmTMB(num~scaled.water_pca*scaled.niche_water_pca+scaled.temp_pca*scaled.niche_temp_pca+scaled.date+sp.layer+(scaled.water_pca+scaled.temp_pca||species)+(scaled.niche_water_pca+scaled.niche_temp_pca||polygon_name),
             data=NMI_analysis,family="binomial")
summary(m)
car::Anova(m)
#library(blme)
#m1 <- bglmer(num~scaled.water_pca+scaled.niche_water_pca+scaled.temp_pca+scaled.niche_temp_pca+scaled.date+sp.layer+(scaled.water_pca+scaled.temp_pca||species)+(scaled.niche_water_pca+scaled.niche_temp_pca||polygon_name),
#             data=NMI_analysis,family="binomial",fixef.prior = normal(cov = diag(9,8)),
#            control=glmerControl(optimizer = "bobyqa"))
#summary(m1)

#residuals(m1)

library(DHARMa)
library(geosphere)
sres <- simulateResiduals(m)
sres_spatial <- recalculateResiduals(sres,NMI_analysis$polygon_name)
polygon_unique <- NMI_analysis[!duplicated(NMI_analysis$polygon_name),c("polygon_name","x","y")]
polygon_unique <- polygon_unique[order(polygon_unique$polygon_name),]
polygon_dist <- distm(polygon_unique[,c("x","y")])
testSpatialAutocorrelation(sres_spatial,distMat=polygon_dist)

phyD <- cophenetic(keep.tip(genus_tree,NMI_analysis$species))
sres_phylo <- recalculateResiduals(sres,NMI_analysis$species)
phyD <- phyD[order(rownames(phyD)),order(colnames(phyD))]
testSpatialAutocorrelation(sres_phylo,distMat=phyD)
####################################################
library(rnaturalearth)

country <- ne_countries(scale=10)
country <- vect(country)
country_df <- as.data.frame(country)

bentity.shp <- vect("C:/Users/pakno/OneDrive - University of Toronto/Bentity2_shapefile_fullres.shp")

bentity.df <- as.data.frame(geom(bentity.shp))
bentity.df$id <- paste0(bentity.df$geom,"_",bentity.df$part)
bentity.df <- as.data.frame(bentity.df)

centroid <- centroids(bentity.shp,inside=T)
centroid_df <- terra::extract(country,centroid)
bentity.shp.df <- as.data.frame(bentity.shp)
bentity.shp.df$Final_country <- centroid_df$brk_name 

bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Tuamotu Islands"] <- "Fr. Polynesia"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Wake Atoll"] <- "U.S. Minor Outlying Is."
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Baker Island"] <- "U.S. Minor Outlying Is."
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Cocos (Keeling) Islands"] <- "Indian Ocean Ter."
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Howland Island"] <- "U.S. Minor Outlying Is."
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Lakshadweep"] <- "India"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Maldives"] <- "Maldives"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Marshall Islands"] <- "Marshall Is."
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Seychelles"] <- "Seychelles"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Tokelau"] <- "Tokelau"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Tromelin Island"] <- "French Southern Territories"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Tuvalu"] <- "Tuvalu"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Chagos Archipelago"] <- "Br. Indian Ocean Ter."
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Kermadec Islands"] <- "New Zealand"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Line Islands"] <- "Kiribati"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Trobriand Islands"] <- "Papua New Guinea"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Senkaku Islands"] <- "Japan"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Pontine Islands"] <- "Italy"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Krakatau Islands"] <- "Indonesia"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Juan Fernandez Islands"] <- "Chile"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Gambier Islands"] <- "Fr. Polynesia"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Galite Islands"] <- "Tunisia"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Cocos Island (Costa Rica)"] <- "Costa Rica"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Wallis and Futuna"] <- "Wallis and Futuna Is."
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Paracel Islands"] <- "China"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Coral Sea Islands"] <- "Coral Sea Is."
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Clipperton Island"] <- "Clipperton I."
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Ashmore and Cartier Islands"] <- "Ashmore and Cartier Is."

bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Cameroon line Archipelago"] <- NA
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Sudan"] <- NA

bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Ivory Coast"] <- "Cote d'Ivoire"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Western Sahara"] <- "W. Sahara"
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Israel and Palestine"] <- NA
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Serbia"] <- NA
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Hispaniola"] <- NA
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Borneo"] <- NA
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Guernsey"] <- "Jersey" #Jersey means Channel Islands
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Samoan Islands"] <- NA
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Cyprus"] <- NA
bentity.shp.df$Final_country[bentity.shp.df$BENTITY2_N == "Lesser Antilles"] <- NA

######Lesser sunda islands as Indonesia, Aland as Finland, Malyesia & Singapore = Maylesia, 
capacity <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/capacity.csv",encoding="UTF-8", stringsAsFactors=FALSE)

not_found <- bentity.shp.df$Final_country[!bentity.shp.df$Final_country %in% capacity$Country_global]
not_found <- sort(unique(not_found))

not_found <- capacity$Country_global[!capacity$Country_global %in% bentity.shp.df$Final_country]
not_found <- sort(unique(not_found))
######################################################
#not sure why predict function doesn't work - calculate manually heere

NMI_analysis$dummy <- 1
future_df <- NMI_analysis[,c("species","polygon_name","shp.ID","sp.layer","scaled.date","scaled.niche_temp_pca","scaled.niche_water_pca","dummy")]
future_temp <- t(apply(NMI_analysis[,c("MAT.2C","Temp_sd.2C")],1, function(x) (x-colMeans(NMI_analysis[,c("MAT","Temp_sd")]))/apply(NMI_analysis[,c("MAT","Temp_sd")],2,sd)))
colnames(future_temp) <- c("MAT","Temp_sd")
future_df$temp_pca <- apply(future_temp,1, function(x) rowSums(x*t(NMI_analysis_pca_temp$rotation[,1])))
future_df$scaled.temp_pca <- (future_df$temp_pca - mean(NMI_analysis$temp_pca))/sd(NMI_analysis$temp_pca)
future_df$scaled.temp_pca_current <- NMI_analysis$scaled.temp_pca

#not sure why predict function doesn't work - calculate manually heere
future_water <- t(apply(NMI_analysis[,c("Prec.2C","Prec_sd.2C","SM.2C","SM_sd.2C")],1, function(x) (x-colMeans(NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")]))/apply(NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")],2,sd)))
colnames(future_water) <- c("Prec","Prec_sd","SM","SM_sd")
future_df$water_pca <- apply(future_water,1, function(x) rowSums(x*t(NMI_analysis_pca_water$rotation[,1])))
future_df$scaled.water_pca <- (future_df$water_pca - mean(NMI_analysis$water_pca))/sd(NMI_analysis$water_pca)
future_df$scaled.water_pca_current <- NMI_analysis$scaled.water_pca

future_df$future_status <- predict(m,future_df,type="response")
#future_df$future_status_bayes <- predict(m1,future_df,type="response")

future_df$current_status_projection <- predict(m,type="response")
future_df$current_status <- NMI_analysis$num

future_df$diff <- future_df$future_status-future_df$current_status_projection
future_df$Final_country <- bentity.shp.df[match(future_df$polygon_name,bentity.shp.df$BENTITY2_N),"Final_country"]

library(tidyverse)

site_summary <- future_df %>% group_by(polygon_name,shp.ID) %>% summarize(future.sum=sum(future_status),
                                                                   current.sum = sum(current_status),
                                                                   current.sum.projection = sum(current_status_projection),
                                                                   diff = mean(diff),
                                                                   n = sum(dummy))

site_summary <- rbind(as.data.frame(site_summary),data.frame(polygon_name=bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary$polygon_name,"BENTITY2_N"],
                                                                               shp.ID = rownames(bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary$polygon_name,]),
                                                                               future.sum = 0,
                                                                               current.sum = 0,
                                                                               current.sum.projection = 0,
                                                                               diff = 0,
                                                                               n = 0))

site_summary$proj_diff <- site_summary$future.sum-site_summary$current.sum.projection

sp_summary <- future_df %>% group_by(species) %>% summarize(future.sum=sum(future_status),
                                                            current.sum = sum(current_status),
                                                            current.sum.projection = sum(current_status_projection),
                                                            diff = mean(diff),
                                                            n = sum(dummy))
sp_summary$proj_diff <- sp_summary$future.sum-sp_summary$current.sum.projection

invasive <- c("Acromyrmex_octospinosus","Anoplolepis_gracilipes","Brachyponera_chinensis","Lasius_neglectus","Linepithema_humile",
              "Monomorium_destructor","Monomorium_floricola","Monomorium_pharaonis","Myrmica_rubra","Nylanderia_pubens","Paratrechina_longicornis","Pheidole_megacephala","Solenopsis_invicta",
              "Solenopsis_geminata","Solenopsis_papuana","Solenopsis_richteri","Tapinoma_melanocephalum","Technomyrmex_albipes","Wasmannia_auropunctata")
invasive_df <- future_df[future_df$species %in% invasive,]

site_summary_invasive <- invasive_df %>% group_by(polygon_name,shp.ID) %>% summarize(future.sum=sum(future_status),
                                                                   current.sum = sum(current_status),
                                                                   current.sum.projection = sum(current_status_projection),
                                                                   diff = mean(diff),
                                                                   n = sum(dummy))

site_summary_invasive <- rbind(as.data.frame(site_summary_invasive),data.frame(polygon_name=bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary_invasive$polygon_name,"BENTITY2_N"],
                                                                                     shp.ID = rownames(bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary_invasive$polygon_name,]),
                                                                                     future.sum = 0,
                                                                                     current.sum = 0,
                                                                                     current.sum.projection = 0,
                                                                                     diff = 0,
                                                                                     n = 0))

site_summary_invasive$proj_diff <- site_summary_invasive$future.sum-site_summary_invasive$current.sum.projection

sp_summary_invasive <- invasive_df %>% group_by(species) %>% summarize(future.sum=sum(future_status),
                                                            current.sum = sum(current_status),
                                                            current.sum.projection = sum(current_status_projection),
                                                            diff = mean(diff),
                                                            n = sum(dummy))
sp_summary_invasive$proj_diff <- sp_summary_invasive$future.sum-sp_summary_invasive$current.sum.projection

site_summary$proj_diff <- site_summary$future.sum-site_summary$current.sum.projection

#future_df$Final_country <- bentity.shp.df[match(future_df$polygon_name,bentity.shp.df$BENTITY2_N),"Final_country"]
#Country_2C <- aggregate(cbind(future_status,current_status_projection)~Final_country+species, FUN=function(x) prod(1-x),data=future_df)
#Country_2C[,c("future_status","current_status_projection")] <- 1-Country_2C[,c("future_status","current_status_projection")]
#Country_2C$proj_diff <- Country_2C$future_status-Country_2C$current_status_projection

#Country_summary_2C <- Country_2C %>% group_by(Final_country) %>% summarize(future.sum=sum(future_status),
                                                                               #current.sum.projection = sum(current_status_projection),
                                                                               #proj_diff = sum(proj_diff))
#Country_summary_2C_invasive <- Country_2C[Country_2C$species %in% invasive,] %>% group_by(Final_country) %>% summarize(future.sum=sum(future_status),
                                                     #current.sum.projection = sum(current_status_projection),
                                                     #proj_diff = sum(proj_diff))

site_summary$Final_country <- bentity.shp.df[match(site_summary$polygon_name,bentity.shp.df$BENTITY2_N),"Final_country"]
site_summary_invasive$Final_country <- bentity.shp.df[match(site_summary_invasive$polygon_name,bentity.shp.df$BENTITY2_N),"Final_country"]

Country_summary_2C <- site_summary %>% group_by(Final_country) %>% summarize(proj_diff = max(proj_diff))
Country_summary_2C_invasive <- site_summary_invasive %>% group_by(Final_country) %>% summarize(proj_diff = max(proj_diff))

##################### 4C Projection

NMI_analysis$dummy <- 1
future_df <- NMI_analysis[,c("species","polygon_name","shp.ID","sp.layer","scaled.date","scaled.niche_temp_pca","scaled.niche_water_pca","dummy")]
future_temp <- t(apply(NMI_analysis[,c("MAT.4C","Temp_sd.4C")],1, function(x) (x-colMeans(NMI_analysis[,c("MAT","Temp_sd")]))/apply(NMI_analysis[,c("MAT","Temp_sd")],2,sd)))
colnames(future_temp) <- c("MAT","Temp_sd")
future_df$temp_pca <- apply(future_temp,1, function(x) rowSums(x*t(NMI_analysis_pca_temp$rotation[,1])))
future_df$scaled.temp_pca <- (future_df$temp_pca - mean(NMI_analysis$temp_pca))/sd(NMI_analysis$temp_pca)
future_df$scaled.temp_pca_current <- NMI_analysis$scaled.temp_pca

#not sure why predict function doesn't work - calculate manually heere
future_water <- t(apply(NMI_analysis[,c("Prec.4C","Prec_sd.4C","SM.4C","SM_sd.4C")],1, function(x) (x-colMeans(NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")]))/apply(NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")],2,sd)))
colnames(future_water) <- c("Prec","Prec_sd","SM","SM_sd")
future_df$water_pca <- apply(future_water,1, function(x) rowSums(x*t(NMI_analysis_pca_water$rotation[,1])))
future_df$scaled.water_pca <- (future_df$water_pca - mean(NMI_analysis$water_pca))/sd(NMI_analysis$water_pca)
future_df$scaled.water_pca_current <- NMI_analysis$scaled.water_pca

future_df$future_status <- predict(m,future_df,type="response")
#future_df$future_status_bayes <- predict(m1,future_df,type="response")

future_df$current_status_projection <- predict(m,type="response")
future_df$current_status <- NMI_analysis$num

future_df$diff <- future_df$future_status-future_df$current_status_projection

library(tidyverse)

indoor_df <- subset(future_df)

site_summary_4C <- indoor_df %>% group_by(polygon_name,shp.ID) %>% summarize(future.sum=sum(future_status),
                                                                      current.sum = sum(current_status),
                                                                      current.sum.projection = sum(current_status_projection),
                                                                      diff = mean(diff),
                                                                      n = sum(dummy))
site_summary_4C <- rbind(as.data.frame(site_summary_4C),data.frame(polygon_name=bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary_4C$polygon_name,"BENTITY2_N"],
           shp.ID = rownames(bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary_4C$polygon_name,]),
           future.sum = 0,
           current.sum = 0,
           current.sum.projection = 0,
           diff = 0,
           n = 0))

site_summary_4C$proj_diff <- site_summary_4C$future.sum-site_summary_4C$current.sum.projection

sp_summary_4C <- indoor_df %>% group_by(species) %>% summarize(future.sum=sum(future_status),
                                                               current.sum = sum(current_status),
                                                               current.sum.projection = sum(current_status_projection),
                                                               diff = mean(diff),
                                                               n = sum(dummy))
sp_summary_4C$proj_diff <- sp_summary_4C$future.sum-sp_summary_4C$current.sum.projection

invasive <- c("Acromyrmex_octospinosus","Anoplolepis_gracilipes","Brachyponera_chinensis","Lasius_neglectus","Linepithema_humile",
              "Monomorium_destructor","Monomorium_floricola","Monomorium_pharaonis","Myrmica_rubra","Nylanderia_pubens","Paratrechina_longicornis","Pheidole_megacephala","Solenopsis_invicta",
              "Solenopsis_geminata","Solenopsis_papuana","Solenopsis_richteri","Tapinoma_melanocephalum","Technomyrmex_albipes","Wasmannia_auropunctata")
indoor_df <- future_df[future_df$species %in% invasive,]

### Aggregate to country level unique invasive species
site_summary_invasive_4C <- indoor_df %>% group_by(polygon_name,shp.ID) %>% summarize(future.sum=sum(future_status),
                                                                      current.sum = sum(current_status),
                                                                      current.sum.projection = sum(current_status_projection),
                                                                      diff = mean(diff),
                                                                      n = sum(dummy))

site_summary_invasive_4C <- rbind(as.data.frame(site_summary_invasive_4C),data.frame(polygon_name=bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary_invasive_4C$polygon_name,"BENTITY2_N"],
                                                                   shp.ID = rownames(bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary_invasive_4C$polygon_name,]),
                                                                   future.sum = 0,
                                                                   current.sum = 0,
                                                                   current.sum.projection = 0,
                                                                   diff = 0,
                                                                   n = 0))

site_summary_invasive_4C$proj_diff <- site_summary_invasive_4C$future.sum-site_summary_invasive_4C$current.sum.projection
site_summary_invasive_4C$warming_diff <- site_summary_invasive_4C$future.sum-site_summary_invasive$future.sum
site_summary_4C$warming_diff <- site_summary_4C$future.sum-site_summary$future.sum

sp_summary_invasive_4C <- indoor_df %>% group_by(species) %>% summarize(future.sum=sum(future_status),
                                                               current.sum = sum(current_status),
                                                               current.sum.projection = sum(current_status_projection),
                                                               diff = mean(diff),
                                                               n = sum(dummy))
sp_summary_invasive_4C$proj_diff <- sp_summary_invasive_4C$future.sum-sp_summary_invasive_4C$current.sum.projection

#future_df$Final_country <- bentity.shp.df[match(future_df$polygon_name,bentity.shp.df$BENTITY2_N),"Final_country"]
#Country_4C <- aggregate(cbind(future_status,current_status_projection)~Final_country+species, FUN=function(x) prod(1-x),data=future_df)
#Country_4C[,c("future_status","current_status_projection")] <- 1-Country_4C[,c("future_status","current_status_projection")]
#Country_4C$proj_diff <- Country_4C$future_status-Country_4C$current_status_projection

site_summary_4C$Final_country <- bentity.shp.df[match(site_summary_4C$polygon_name,bentity.shp.df$BENTITY2_N),"Final_country"]
site_summary_invasive_4C$Final_country <- bentity.shp.df[match(site_summary_invasive_4C$polygon_name,bentity.shp.df$BENTITY2_N),"Final_country"]

Country_summary_4C <- site_summary_4C %>% group_by(Final_country) %>% summarize(proj_diff = max(proj_diff),
                                                                           warming_diff = max(warming_diff))
Country_summary_4C_invasive <- site_summary_4C %>% group_by(Final_country) %>% summarize(proj_diff = max(proj_diff),
                                                                                warming_diff = max(warming_diff))



###############################
library(ggplot2)

indoor.sr <- subset(clim_invasion_df,num==0) %>% group_by(shp.ID) %>% count(num)
bentity.df$indoor.sr <- unlist(indoor.sr[match(bentity.df$geom,indoor.sr$shp.ID),"n"])

outdoor.sr <- subset(clim_invasion_df,num==1) %>% group_by(shp.ID) %>% count(num)
bentity.df$outdoor.sr <- unlist(outdoor.sr[match(bentity.df$geom,outdoor.sr$shp.ID),"n"])

theme <- theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="bottom",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

p1 <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=indoor.sr))+
  geom_polygon(colour="black")+
  scale_fill_continuous(low="blue",high="red",na.value="white")+
  theme+
  coord_equal()
plot(p1)

p1b <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=outdoor.sr))+
  geom_polygon(colour="black")+
  labs(fill="Alien ant species richness")+
  scale_fill_continuous(low="blue",high="red",na.value="white")+
  theme+
  coord_equal()
plot(p1b)

#####################################

site_summary <- as.data.frame(site_summary)
site_summary_4C <- as.data.frame(site_summary_4C)
site_summary_invasive <- as.data.frame(site_summary_invasive)
site_summary_invasive_4C <- as.data.frame(site_summary_invasive_4C)

bentity.df$future_gain_2C <- site_summary[match(bentity.df$geom,(site_summary$shp.ID)),"proj_diff"]
bentity.df$future_gain_4C <- site_summary_4C[match(bentity.df$geom,site_summary_4C$shp.ID),"proj_diff"]
bentity.df$future_gain_invasive_2C <- site_summary_invasive[match(bentity.df$geom,(site_summary_invasive$shp.ID)),"proj_diff"]
bentity.df$future_gain_invasive_4C <- site_summary_invasive_4C[match(bentity.df$geom,site_summary_invasive_4C$shp.ID),"proj_diff"]
bentity.df$warming_diff <- site_summary_4C[match(bentity.df$geom,(site_summary_4C$shp.ID)),"warming_diff"]
bentity.df$warming_diff_invasive <- site_summary_invasive_4C[match(bentity.df$geom,site_summary_invasive_4C$shp.ID),"warming_diff"]

p3a <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_2C))+
  geom_polygon(colour="black")+
  labs(fill="Gain in alien ant species under 2°C warming")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,9))+
  theme+
  coord_equal()
plot(p3a)

p3b <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_4C))+
  geom_polygon(colour="black")+
  labs(fill="Gain in alien ant species under 4°C warming")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,9))+
  theme+
  coord_equal()
plot(p3b)

p3c <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_invasive_2C))+
  geom_polygon(colour="black")+
  labs(fill="Gain in invasive ant species under 2°C warming")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,2))+
  theme+
  coord_equal()
plot(p3c)

p3d <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_invasive_4C))+
  geom_polygon(colour="black")+
  labs(fill="Gain in invasive ant species under 4°C warming")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,2))+
  theme+
  coord_equal()
plot(p3d)

p3e <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=warming_diff))+
  geom_polygon(colour="black")+
  labs(fill="Differences in alien species richness gain between 4°C and 2°C warming")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,5.8))+
  theme+
  coord_equal()
plot(p3e)

p3f <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=warming_diff_invasive))+
  geom_polygon(colour="black")+
  labs(fill="Differences in invasive species richness gain between 4°C and 2°C warming")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,1.2))+
  theme+
  coord_equal()
plot(p3f)

library(ggpubr)
p3 <- ggarrange(p3a,p3c,p3b,p3d,p3e,p3f,nrow=3,ncol=2)
plot(p3)

ggsave("maps.tiff",width=24,height=12,compression="lzw")
