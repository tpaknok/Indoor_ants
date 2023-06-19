setwd("C:/Users/pakno/OneDrive - University of Toronto")
source("Indoor_ants/Indoor_ants/make_genus_tree.R")

library(phytools)
library(ape)
library(stringr)
library(tidyverse)
library(terra)
library(performance)
clim_invasion_df <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/clim_invasion_df_14022023.csv")
clim_invasion_df_f <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/clim_invasion_df_f_14022023.csv")
clim_invasion_df <- cbind(clim_invasion_df,clim_invasion_df_f[,6:17])
clim_invasion_df$num <- ifelse(clim_invasion_df$num ==0,0,1)

NMI_analysis <- subset(clim_invasion_df)
#indoor_sp <- unique(NMI_analysis[NMI_analysis$num == 0,"species"])
#NMI_analysis <- NMI_analysis[NMI_analysis$species %in% indoor_sp,]
#NMI_analysis <- na.omit(NMI_analysis)

ant_tree <- org_ant_tree <- read.tree("C:/Users/pakno/OneDrive - University of Toronto/ant_tree/backbone_MLtree_RaxML.tre")
is.ultrametric(ant_tree)
#ant_tree <- force.ultrametric(ant_tree,method="nnls")
ant_tree <- as.phylo(hclust(as.dist(cophenetic(ant_tree)),method="average"))
dendextend::cor_cophenetic(ant_tree,org_ant_tree)
ant_tree$node.label <- 1:length(ant_tree$node.label)

ant_tree$tip.label <- gsub("[\"]","",ant_tree$tip.label)
ant_tree

ant_placement <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/ant_genus_placement.csv")
#genus_tree <- make_genus_tree(ant_tree,ant_placement,unique(NMI_temp$sp))
genus_tree <- make_genus_tree(ant_tree,ant_placement,unique(NMI_analysis$species))

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

m <- glmmTMB(num~scaled.water_pca*scaled.niche_water_pca+scaled.temp_pca*scaled.niche_temp_pca+scaled.date+sp.layer+(scaled.water_pca+scaled.temp_pca||species)+(scaled.niche_water_pca+scaled.niche_temp_pca||polygon_name),
             data=NMI_analysis,family="binomial")
summary(m)
car::Anova(m)
r2(m)
write.csv(round(summary(m)$coef$cond,3),"glmm.csv")
#library(blme)
#m1 <- bglmer(num~scaled.water_pca+scaled.niche_water_pca+scaled.temp_pca+scaled.niche_temp_pca+scaled.date+Generalist+Arboreal+(scaled.water_pca+scaled.temp_pca||species)+(scaled.niche_water_pca+scaled.niche_temp_pca||polygon_name),
#             data=NMI_analysis,family="binomial",fixef.prior = normal(cov = diag(9,8)),
#            control=glmerControl(optimizer = "bobyqa"))
#summary(m1)

#residuals(m1)

library(DHARMa)
library(geosphere)
m1 <- glmmTMB(num~scaled.water_pca*scaled.niche_water_pca+scaled.temp_pca*scaled.niche_temp_pca+scaled.date+sp.layer+(scaled.water_pca+scaled.temp_pca||species)+(scaled.niche_water_pca+scaled.niche_temp_pca||polygon_name),
             data=NMI_analysis[order(NMI_analysis$polygon_name),],family="binomial") #re-run model but make sure the factors are ordered alphabetically
summary(m1)

sres <- simulateResiduals(m1)
sres_spatial <- recalculateResiduals(sres,NMI_analysis[order(NMI_analysis$polygon_name),"polygon_name"])
polygon_unique <- NMI_analysis[!duplicated(NMI_analysis$polygon_name),c("polygon_name","x","y")] #make sure the structures are consistent between residual objects and distance M
polygon_unique <- polygon_unique[order(polygon_unique$polygon_name),]
polygon_dist <- distm(polygon_unique[,c("x","y")])
testSpatialAutocorrelation(sres_spatial,distMat=polygon_dist)

m2 <- glmmTMB(num~scaled.water_pca*scaled.niche_water_pca+scaled.temp_pca*scaled.niche_temp_pca+scaled.date+sp.layer+(scaled.water_pca+scaled.temp_pca||species)+(scaled.niche_water_pca+scaled.niche_temp_pca||polygon_name),
              data=NMI_analysis[order(NMI_analysis$species),],family="binomial")
summary(m2) #re-run model but make sure the factors are ordered alphabetically
sres <- simulateResiduals(m2)
phyD <- cophenetic(keep.tip(genus_tree,NMI_analysis$species))
sres_phylo <- recalculateResiduals(sres,NMI_analysis[order(NMI_analysis$species),"species"])
phyD <- phyD[order(rownames(phyD)),order(colnames(phyD))]
testSpatialAutocorrelation(sres_phylo,distMat=phyD)
###################################################Fig 2

library(ggeffects)
v1 <- seq(min(NMI_analysis$scaled.temp_pca),max(NMI_analysis$scaled.temp_pca),length.out=100)
v2 <- seq(min(NMI_analysis$scaled.niche_temp_pca),max(NMI_analysis$scaled.niche_temp_pca),length.out=100)

predict_niche <- ggemmeans(m,terms=c("scaled.temp_pca [v1]","scaled.niche_temp_pca [v2]"),type="fixed",rg.limit=1000*1000)
predict_niche$group <- as.numeric(as.character(predict_niche$group))
NMI_analysis$Status <- ifelse(NMI_analysis$num == 0, "Indoor", "Wild")

theme <- theme(axis.line=element_line(colour="black"),
               axis.text = element_text(size=6),
               axis.title = element_text(size=6),
               legend.position="bottom",
               panel.background=element_rect(colour="white",fill="white"),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               plot.background=element_rect(colour="white",fill="white"),
               legend.key.height = unit(0.25, "cm"),
               legend.key.width = unit(0.3, "cm"),
               legend.title= element_text(size=6),
               legend.text = element_text(size=6),
               legend.margin=margin(t=0.2,unit="cm"))

p2b<- ggplot(predict_niche,aes(y=group,x=x))+
  geom_raster(aes(fill=predicted))+
  geom_point(data=NMI_analysis,aes(y=scaled.niche_temp_pca,x=scaled.temp_pca,colour=Status),size=0.5)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b)",size=3)+
  xlab("\n")+
  ylab("Temperature PCA1 of native range \n (Colder and more seasonal)")+
  scale_fill_continuous(low="#ffffb2",high="#e31a1c",breaks=c(0.1,0.5,0.9))+
  scale_color_manual(values=c("grey75","black"))+
  theme

plot(p2b)

library(ggeffects)
v3 <- seq(min(NMI_analysis$scaled.water_pca),max(NMI_analysis$scaled.water_pca),length.out=100)

predict_water <- ggemmeans(m,terms=c("scaled.temp_pca [v1]","scaled.water_pca[v3]"),type="fixed",rg.limit=100000)
predict_water$group <- as.numeric(as.character(predict_water$group))

prop <- NMI_analysis %>% group_by(polygon_name) %>% summarize(n = n(), sum = sum(num)) %>% mutate(prop = sum/n)
NMI_analysis$prop <-  unlist(prop[match(NMI_analysis$polygon_name,prop$polygon_name),"prop"])

p2a <- ggplot(predict_water,aes(y=group,x=x))+
  geom_raster(aes(fill=predicted))+
  geom_point(data=NMI_analysis,aes(y=scaled.water_pca,x=scaled.temp_pca,colour=prop),size=0.5)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a)",size=3)+
  xlab("Temperature PCA1 of invaded region \n (Colder and more seasonal)")+
  ylab("Water PCA1 of invaded region \n (Wetter and more seasonal)")+
  scale_fill_continuous(low="#ffffb2",high="#e31a1c",name="Predicted probability of establishing in the wild",breaks=c(0.1,0.5,0.9))+
  scale_color_continuous(low="grey75",high="black",name="Proportion of species established in the wild",breaks=c(0.1,0.5,0.9))+
  theme

plot(p2a)

library(ggpubr)
p2 <- ggarrange(p2a,p2b,common.legend=T,legend="bottom")
plot(p2)
ggsave("Indoor_ants/Fig2.tiff",dpi=800,height=8.4,width=16.8,units="cm",compression="lzw",bg="white")
###################################################
library(rnaturalearth)

capacity <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/capacity.csv")
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

not_found <- bentity.shp.df$Final_country[!bentity.shp.df$Final_country %in% capacity$Country_global]
not_found <- sort(unique(not_found))

not_found <- capacity$Country_global[!capacity$Country_global %in% bentity.shp.df$Final_country]
not_found <- sort(unique(not_found))

NMI_analysis$dummy <- 1
######################################################

future_df <- NMI_analysis[,c("species","polygon_name","shp.ID","sp.layer","scaled.date","scaled.niche_temp_pca","scaled.niche_water_pca","dummy")]

source("Indoor_ants/Indoor_ants/rescale_future_climate.R")
NMI_analysis$scaled.temp_pca_2C <- rescale_future_climate(NMI_analysis[,c("MAT.2C","Temp_sd.2C")],
                                                       NMI_analysis[,c("MAT","Temp_sd")],
                                                       NMI_analysis_pca_temp,
                                                       NMI_analysis$temp_pca)
NMI_analysis$scaled.temp_pca_4C <- rescale_future_climate(NMI_analysis[,c("MAT.4C","Temp_sd.4C")],
                                                       NMI_analysis[,c("MAT","Temp_sd")],
                                                       NMI_analysis_pca_temp,
                                                       NMI_analysis$temp_pca)
NMI_analysis$scaled.water_pca_2C <- rescale_future_climate(NMI_analysis[,c("Prec.2C","Prec_sd.2C","SM.2C","SM_sd.2C")],
                                                        NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")],
                                                        NMI_analysis_pca_water,
                                                        NMI_analysis$water_pca)
NMI_analysis$scaled.water_pca_4C<- rescale_future_climate(NMI_analysis[,c("Prec.4C","Prec_sd.4C","SM.4C","SM_sd.4C")],
                                                       NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")],
                                                       NMI_analysis_pca_water,
                                                       NMI_analysis$water_pca)

NMI_analysis$current_status_projection <- predict(m,type="response")

predict_data <- list(NMI_analysis[,c("species","polygon_name","scaled.date","scaled.niche_temp_pca","scaled.niche_water_pca","scaled.temp_pca_2C","scaled.water_pca_2C","sp.layer")],
                     NMI_analysis[,c("species","polygon_name","scaled.date","scaled.niche_temp_pca","scaled.niche_water_pca","scaled.temp_pca_4C","scaled.water_pca_4C","sp.layer")])

predict_data <- lapply(predict_data,setNames,c("species","polygon_name","scaled.date","scaled.niche_temp_pca","scaled.niche_water_pca","scaled.temp_pca","scaled.water_pca","sp.layer"))

future_occ <- do.call(cbind,lapply(predict_data,function(x) predict(m,x,type="response")))
colnames(future_occ) <- c("future_status_2C","future_status_4C")
NMI_analysis <- cbind(NMI_analysis,future_occ)

invasive <- c("Anoplolepis_gracilipes","Azteca_sericeasur","Brachyponera_chinensis","Camponotus_conspicuus_zonatus","Cardiocondyla_wroughtonii","Cardiocondyla_emeryi",
              "Formica_aquilonia","Formica_paralugubris",
              "Lasius_neglectus","Linepithema_humile",
              "Monomorium_floricola","Monomorium_monomorium","Monomorium_pharaonis","Myrmica_rubra","Nylanderia_bourbonica","Nylanderia_fulva",
              "Plagiolepis_alluaudi","Paratrechina_longicornis","Pheidole_megacephala","Pheidole_radoszkowskii","Solenopsis_invicta",
              "Solenopsis_geminata","Solenopsis_papuana","Solenopsis_richteri","Tapinoma_melanocephalum","Technomyrmex_albipes","Technomyrmex difficilis","Tetramorium bicarinatum","Tetramorium simillimum",
              "Trichomyrmex_destructor","Wasmannia_auropunctata")

NMI_analysis$invasive <- ifelse(NMI_analysis$species %in% invasive,"Invasive","Alien")

NMI_analysis$proj_diff_2C <- NMI_analysis$future_status_2C-NMI_analysis$current_status_projection
NMI_analysis$proj_diff_4C <- NMI_analysis$future_status_4C-NMI_analysis$current_status_projection
NMI_analysis$warming_diff <- NMI_analysis$future_status_4C-NMI_analysis$future_status_2C

######################################################
site_summary <- NMI_analysis %>% group_by(polygon_name,shp.ID) %>% summarize(future_sum_2C=sum(future_status_2C),
                                                                          future_sum_4C=sum(future_status_4C),
                                                                          current_sum = sum(num),
                                                                          current_indoor = sum(dummy)-sum(num),
                                                                          current_projection_sum = sum(current_status_projection),
                                                                          proj_diff_2C_net = sum(proj_diff_2C),
                                                                          proj_diff_4C_net = sum(proj_diff_4C),
                                                                          warming_diff_net = sum(warming_diff),
                                                                          n = sum(dummy),
                                                                          future_sum_invasive_2C=sum(future_status_2C[invasive=="Invasive"]),
                                                                          future_sum_invasive_4C=sum(future_status_4C[invasive=="Invasive"]),
                                                                          current_invasive_sum = sum(num[invasive=="Invasive"]),
                                                                          current_projection_invasive_sum = sum(current_status_projection[invasive=="Invasive"]),
                                                                          proj_diff_invasive_2C_net = sum(proj_diff_2C[invasive=="Invasive"]),
                                                                          proj_diff_invasive_4C_net = sum(proj_diff_4C[invasive=="Invasive"]),
                                                                          warming_diff_invasive_net = sum(warming_diff[invasive=="Invasive"]),
                                                                          n_invasive = sum(dummy[invasive=="Invasive"]))

empty_poly <- data.frame(polygon_name=bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary$polygon_name,"BENTITY2_N"],
           shp.ID = rownames(bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary$polygon_name,]))

site_summary <- plyr::rbind.fill(site_summary,empty_poly)
site_summary[is.na(site_summary)] <- 0


######################################################

library(ggplot2)

clim_invasion_df$invasive <- clim_invasion_df$species %in% gsub("_",".",invasive)
indoor.sr <- subset(clim_invasion_df,num==0) %>% group_by(shp.ID) %>% count(num)
indoor.sr.invasive <- subset(clim_invasion_df, num == 0 & invasive == T) %>% group_by(shp.ID) %>% count(num)

bentity.df$indoor.sr <- unlist(indoor.sr[match(bentity.df$geom,indoor.sr$shp.ID),"n"])
bentity.df$indoor.sr.invasive <- unlist(indoor.sr.invasive[match(bentity.df$geom,indoor.sr.invasive$shp.ID),"n"])

outdoor.sr <- subset(clim_invasion_df,num==1) %>% group_by(shp.ID) %>% count(num)
outdoor.sr.invasive <- subset(clim_invasion_df, num == 1 & invasive == T) %>% group_by(shp.ID) %>% count(num)

bentity.df$outdoor.sr <- unlist(outdoor.sr[match(bentity.df$geom,outdoor.sr$shp.ID),"n"])
bentity.df$outdoor.sr.invasive <- unlist(outdoor.sr.invasive[match(bentity.df$geom,outdoor.sr.invasive$shp.ID),"n"])

theme <- theme(axis.line=element_blank(),
               axis.text.x=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               legend.position="bottom",
               panel.background=element_rect(colour="white",fill="white"),
               panel.border=element_blank(),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               plot.background=element_rect(colour="white",fill="white"),
               legend.key.height = unit(0.25, "cm"),
               legend.key.width = unit(0.3, "cm"),
               legend.text = element_text(size=7),
               legend.margin=margin(t=-0.5,unit="cm"),
               legend.background=element_rect(fill="white"))

bentity.df$indoor.sr[bentity.df$indoor.sr == 0] <- NA
bentity.df$outdoor.sr[bentity.df$outdoor.sr == 0] <- NA
bentity.df$indoor.sr.invasive[bentity.df$indoor.sr.invasive == 0] <- NA
bentity.df$outdoor.sr.invasive[bentity.df$outdoor.sr.invasive == 0] <- NA

library(tiff)
library(grid)
indoor_tiff <- readTIFF("Indoor_ants/indoor.tif")
wild_tiff <- readTIFF("Indoor_ants/wild.tif")

g1 <- rasterGrob(indoor_tiff,width = unit(1,"cm"),height=unit(1,"cm"),interpolate=TRUE)
g2 <- rasterGrob(wild_tiff,width = unit(1,"cm"),height=unit(1,"cm"),interpolate=TRUE)

p1a <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=indoor.sr))+
  geom_polygon(colour="black",size=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a) Alien")+
  annotation_custom(g1, xmin=-175, xmax=-175, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0,70),breaks=c(1,35,70))+
  theme

p1b <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=outdoor.sr))+
  geom_polygon(colour="black",size=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b) Alien")+
  annotation_custom(g2, xmin=-175, xmax=-175, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",tran=scales::log1p_trans(),limits=c(0,70),breaks=c(1,35,70))+
  theme

p1c <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=indoor.sr.invasive))+
  geom_polygon(colour="black",size=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(c) Invasive")+
  annotation_custom(g1, xmin=-175, xmax=-175, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0.5,12.5),breaks=c(1,6,12))+
  theme

p1d <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=outdoor.sr.invasive))+
  geom_polygon(colour="black",size=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(d) Invasive")+
  annotation_custom(g2, xmin=-175, xmax=-175, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0.5,12.5),breaks=c(1,6,12))+
  theme

library(ggpubr)
p1_1 <- ggarrange(p1a,p1b,ncol=1,legend="bottom",common.legend=T)
plot(p1_1)
p1_2 <- ggarrange(p1c,p1d,ncol=1,legend="bottom",common.legend=T)
plot(p1_2)

p1 <- ggarrange(p1_1,p1_2)
plot(p1)
ggsave("Indoor_ants/Fig1.tiff",dpi=800,compression="lzw",units="cm",height=12.6,width=16.8,bg="white")
#####################################

site_summary <- as.data.frame(site_summary)

bentity.df$future_gain_2C <- site_summary[match(bentity.df$geom,(site_summary$shp.ID)),"proj_diff_2C_net"]
bentity.df$future_gain_4C <- site_summary[match(bentity.df$geom,site_summary$shp.ID),"proj_diff_4C_net"]
bentity.df$future_gain_invasive_2C <- site_summary[match(bentity.df$geom,(site_summary$shp.ID)),"proj_diff_invasive_2C_net"]
bentity.df$future_gain_invasive_4C <- site_summary[match(bentity.df$geom,site_summary$shp.ID),"proj_diff_invasive_4C_net"]
bentity.df$warming_diff <- site_summary[match(bentity.df$geom,(site_summary$shp.ID)),"warming_diff_net"]
bentity.df$warming_diff_invasive <- site_summary[match(bentity.df$geom,site_summary$shp.ID),"warming_diff_invasive_net"]

theme <- theme(axis.line=element_blank(),
               axis.text.x=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               legend.position="bottom",
               panel.background=element_rect(colour="white",fill="white"),
               panel.border=element_blank(),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               plot.background=element_rect(colour="white",fill="white"),
               legend.key.height = unit(0.25, "cm"),
               legend.key.width = unit(0.3, "cm"),
               legend.text = element_text(size=5),
               legend.margin=margin(t=-0.5,unit="cm"))

library(tiff)
library(grid)
warming_tiff <- readTIFF("Indoor_ants/warming.tif")
g3 <- rasterGrob(warming_tiff,width = unit(0.7,"cm"),height=unit(0.7,"cm"),interpolate=TRUE)

p3a <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_2C))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a) Alien")+
  annotate("text", x = -170, y = -55,label = "2°C",colour="red",size=4,hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,9),breaks = c(0,2,4,6,8))+
  theme
plot(p3a)

p3b <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_4C))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b) Alien")+
  annotate("text", x = -170, y = -55,label = "4°C",colour="red",hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,9),breaks = c(0,2,4,6,8))+
  theme
plot(p3b)

p3c <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_invasive_2C))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(d) Invasive")+
  annotate("text", x = -170, y = -55,label = "2°C",colour="red",hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,3.2),breaks = c(0,1,2,3))+
  theme
plot(p3c)

p3d <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_invasive_4C))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(e) Invasive")+
  annotate("text", x = -170, y = -55,label = "4°C",colour="red",hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,3.2),breaks = c(0,1,2,3))+
  theme
plot(p3d)

p3e <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=warming_diff))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(c) Alien")+
  annotate("text", x = -170, y = -55,label = "4°C vs 2°C",colour="red",hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,6.1),breaks = c(0,2,4,6))+
  theme
plot(p3e)

p3f <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=warming_diff_invasive))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(f) Invasive")+
  annotate("text", x = -170, y = -55,label = "4°C vs 2°C",colour="red",hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,2),breaks = c(0,0.5,1,1.5))+
  theme
plot(p3f)

library(ggpubr)
p3 <- ggarrange(p3a,p3c,p3b,p3d,p3e,p3f,nrow=3,ncol=2)
plot(p3)

ggsave("Indoor_ants/Fig3.tiff",dpi=800,compression="lzw",units="cm",height=16.8,width=16.8,bg="white")


# have to mean them first rather than calculating percent change in each region before averaging - because some regions have inf values (0 invasive currently at outdoor)

site_summary %>% filter(current_indoor >0) %>% reframe(percent_2C = mean(proj_diff_2C_net)/mean(current_projection_sum)*100,
                                                       percent_4C = mean(proj_diff_4C_net)/mean(current_projection_sum)*100,
                                                       percent_2C_invasive = mean(proj_diff_invasive_2C_net)/mean(current_invasive_sum)*100,
                                                       percent_4C_invasive = mean(proj_diff_invasive_4C_net)/mean(current_invasive_sum)*100) 
#################previous clunky code. no longer useful
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
indoor_df <- subset(future_df, current_status==0)
wild_df <- subset(future_df,current_status == 1)
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

invasive <- c("Anoplolepis_gracilipes","Azteca_sericeasur","Brachyponera_chinensis","Camponotus_conspicuus_zonatus","Cardiocondyla_wroughtonii","Cardiocondyla_emeryi",
              "Formica_aquilonia","Formica_paralugubris",
              "Lasius_neglectus","Linepithema_humile",
              "Monomorium_floricola","Monomorium_monomorium","Monomorium_pharaonis","Myrmica_rubra","Nylanderia_bourbonica","Nylanderia_fulva",
              "Plagiolepis_alluaudi","Paratrechina_longicornis","Pheidole_megacephala","Pheidole_radoszkowskii","Solenopsis_invicta",
              "Solenopsis_geminata","Solenopsis_papuana","Solenopsis_richteri","Tapinoma_melanocephalum","Technomyrmex_albipes","Technomyrmex difficilis","Tetramorium bicarinatum","Tetramorium simillimum",
              "Trichomyrmex_destructor","Wasmannia_auropunctata")
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

site_summary <- site_summary[order(site_summary$polygon_name),]
site_summary_invasive <- site_summary_invasive[order(site_summary_invasive$polygon_name),]

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

Country_summary_2C <- site_summary %>% group_by(Final_country) %>% summarize(proj_diff = max(proj_diff),
                                                                             indoor_sr = max(n-current.sum))
Country_summary_2C_invasive <- site_summary_invasive %>% group_by(Final_country) %>% summarize(proj_diff = max(proj_diff),
                                                                                               indoor_sr = max(n-current.sum))

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

invasive <- c("Anoplolepis_gracilipes","Azteca_sericeasur","Brachyponera_chinensis","Camponotus_conspicuus_zonatus","Cardiocondyla_wroughtonii","Cardiocondyla_emeryi",
              "Formica_aquilonia","Formica_paralugubris",
              "Lasius_neglectus","Linepithema_humile",
              "Monomorium_floricola","Monomorium_monomorium","Monomorium_pharaonis","Myrmica_rubra","Nylanderia_bourbonica","Nylanderia_fulva",
              "Plagiolepis_alluaudi","Paratrechina_longicornis","Pheidole_megacephala","Pheidole_radoszkowskii","Solenopsis_invicta",
              "Solenopsis_geminata","Solenopsis_papuana","Solenopsis_richteri","Tapinoma_melanocephalum","Technomyrmex_albipes","Technomyrmex difficilis","Tetramorium bicarinatum","Tetramorium simillimum",
              "Trichomyrmex_destructor","Wasmannia_auropunctata")

#invasive <- c("Acromyrmex_octospinosus","Anoplolepis_gracilipes","Brachyponera_chinensis","Lasius_neglectus","Linepithema_humile",
#"Trichomyrmex_destructor","Monomorium_floricola","Monomorium_pharaonis","Myrmica_rubra","Nylanderia_pubens","Paratrechina_longicornis","Pheidole_megacephala","Solenopsis_invicta",
#"Solenopsis_geminata","Solenopsis_papuana","Solenopsis_richteri","Tapinoma_melanocephalum","Technomyrmex_albipes","Wasmannia_auropunctata")
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
site_summary_4C <- site_summary_4C[order(site_summary_4C$polygon_name),]
site_summary_invasive_4C <- site_summary_invasive_4C[order(site_summary_invasive_4C$polygon_name),]

Country_summary_4C <- site_summary_4C %>% group_by(Final_country) %>% summarize(proj_diff = max(proj_diff),
                                                                                indoor_sr = max(n-current.sum))
Country_summary_4C_invasive <- site_summary_invasive_4C %>% group_by(Final_country) %>% summarize(proj_diff = max(proj_diff),
                                                                                                  indoor_sr = max(n-current.sum))

Country_summary_4C$warming_diff <- Country_summary_4C$proj_diff-Country_summary_2C$proj_diff
Country_summary_4C_invasive$warming_diff <- Country_summary_4C_invasive$proj_diff-Country_summary_2C_invasive$proj_diff
site_summary_4C$warming_diff <- site_summary_4C$proj_diff-site_summary$proj_diff
site_summary_invasive_4C$warming_diff <- site_summary_invasive_4C$proj_diff-site_summary_invasive$proj_diff

mean(site_summary_4C$proj_diff[site_summary_4C$n - site_summary_4C$current.sum > 0])
mean(site_summary$proj_diff[site_summary$n - site_summary$current.sum > 0])
mean(site_summary_invasive_4C$proj_diff[site_summary_invasive_4C$n - site_summary_invasive_4C$current.sum > 0])
mean(site_summary_invasive$proj_diff[site_summary_invasive$n - site_summary_invasive$current.sum > 0])
mean(site_summary_invasive_4C$warming_diff[site_summary_invasive$n - site_summary_invasive$current.sum > 0])
mean(site_summary_4C$warming_diff[site_summary$n - site_summary$current.sum > 0])

mean(site_summary_4C$proj_diff[site_summary_4C$n - site_summary_4C$current.sum == 0])
mean(site_summary$proj_diff[site_summary$n - site_summary$current.sum == 0])
mean(site_summary_invasive_4C$proj_diff[site_summary_invasive_4C$n - site_summary_invasive_4C$current.sum == 0])
mean(site_summary_invasive$proj_diff[site_summary_invasive$n - site_summary_invasive$current.sum == 0])
mean(site_summary_invasive_4C$warming_diff[site_summary_invasive$n - site_summary_invasive$current.sum == 0])
mean(site_summary_4C$warming_diff[site_summary$n - site_summary$current.sum == 0])

###############################
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

theme <- theme(axis.line=element_blank(),
               axis.text.x=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               legend.position="bottom",
               panel.background=element_rect(colour="white",fill="white"),
               panel.border=element_blank(),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               plot.background=element_rect(colour="white",fill="white"),
               legend.key.height = unit(0.25, "cm"),
               legend.key.width = unit(0.3, "cm"),
               legend.text = element_text(size=5),
               legend.margin=margin(t=-0.5,unit="cm"))

warming_tiff <- readTIFF("Indoor_ants/warming.tif")
g3 <- rasterGrob(warming_tiff,width = unit(0.7,"cm"),height=unit(0.7,"cm"),interpolate=TRUE)

p3a <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_2C))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a) Alien")+
  annotate("text", x = -170, y = -55,label = "2°C",colour="red",size=4,hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,9),breaks = c(0,2,4,6,8))+
  theme
plot(p3a)

p3b <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_4C))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b) Alien")+
  annotate("text", x = -170, y = -55,label = "4°C",colour="red",hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,9),breaks = c(0,2,4,6,8))+
  theme
plot(p3b)

p3c <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_invasive_2C))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(d) Invasive")+
  annotate("text", x = -170, y = -55,label = "2°C",colour="red",hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,3.2),breaks = c(0,1,2,3))+
  theme
plot(p3c)

p3d <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=future_gain_invasive_4C))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(e) Invasive")+
  annotate("text", x = -170, y = -55,label = "4°C",colour="red",hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,3.2),breaks = c(0,1,2,3))+
  theme
plot(p3d)

p3e <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=warming_diff))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(c) Alien")+
  annotate("text", x = -170, y = -55,label = "4°C vs 2°C",colour="red",hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,6.1),breaks = c(0,2,4,6))+
  theme
plot(p3e)

p3f <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=warming_diff_invasive))+
  geom_polygon(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(f) Invasive")+
  annotate("text", x = -170, y = -55,label = "4°C vs 2°C",colour="red",hjust=0)+
  annotation_custom(g3, xmin=-180, xmax=-180, ymin=-55, ymax=-55) +
  labs(fill="")+
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,2),breaks = c(0,0.5,1,1.5))+
  theme
plot(p3f)

library(ggpubr)
p3 <- ggarrange(p3a,p3c,p3b,p3d,p3e,p3f,nrow=3,ncol=2)
plot(p3)

ggsave("Indoor_ants/Fig3.tiff",dpi=800,compression="lzw",units="cm",height=16.8,width=16.8,bg="white")
