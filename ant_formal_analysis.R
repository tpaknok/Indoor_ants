setwd("C:/Users/pakno/OneDrive - University of Toronto")
source("Indoor_ants/Indoor_ants/make_genus_tree.R")

library(phytools)
library(ape)
library(stringr)
library(tidyverse)
library(terra)
library(performance)
library(sf)
bentity.shp <- vect("C:/Users/pakno/OneDrive - University of Toronto/Bentity2_shapefile_fullres.shp")
clim_invasion_df <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/clim_invasion_df_22082023.csv")

xy <- as.data.frame(centroids(bentity.shp),geom="XY")
clim_invasion_df <- cbind(clim_invasion_df,xy[match(clim_invasion_df$polygon_name,xy$BENTITY2_N),c("x","y")])
#clim_invasion_df <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/clim_invasion_df_30062023.csv")
#clim_invasion_df_f <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/clim_invasion_df_f_14022023.csv")
#clim_invasion_df <- cbind(clim_invasion_df,clim_invasion_df_f[,6:17])
clim_invasion_df$num <- ifelse(clim_invasion_df$num ==0,0,1)

NMI_analysis <- subset(clim_invasion_df)
#indoor_sp <- unique(NMI_analysis[NMI_analysis$num == 0,"species"])
#NMI_analysis <- NMI_analysis[NMI_analysis$species %in% indoor_sp,]
NMI_analysis <- na.omit(NMI_analysis)
NMI_analysis$species <- gsub(" ",".",NMI_analysis$species)
NMI_analysis_sum <- NMI_analysis %>% group_by(polygon_name) %>% dplyr::summarise(total = sum(num),
                                                                                 count = n(),
                                                                                 indoor = count-total)
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
genus_tree$tip.label <-gsub("_",".",genus_tree$tip.label)

######################################################################
#NMI_analysis$date <- scale(log10(NMI_analysis$date))

#NMI_analysis$species <- gsub("[.]","_",NMI_analysis$species)
#NMI_analysis$species <- gsub(" ","_",NMI_analysis$species)
#NMI_analysis <- na.omit(NMI_analysis)
#NMI_analysis_pca_temp <- prcomp(scale(NMI_analysis[,c("tmin","tmax","tmin_sd","tmax_sd")]))
#pcatest <- PCAtest(scale(NMI_analysis[,c("tmin","tmax","tmin_sd","tmax_sd")]))
#NMI_analysis$temp_pca1 <- apply(scale(NMI_analysis[,c("tmin","tmax","tmin_sd","tmax_sd")]),1, function(x) rowSums(x*t(NMI_analysis_pca_temp$rotation[,1])))
#NMI_analysis$temp_pca2 <- apply(scale(NMI_analysis[,c("tmin","tmax","tmin_sd","tmax_sd")]),1, function(x) rowSums(x*t(NMI_analysis_pca_temp$rotation[,2])))

#NMI_analysis_pca_water <- prcomp(scale(NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")]))
#pcatest <- PCAtest(scale(NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")]))
#NMI_analysis$water_pca1 <- apply(scale(NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")]),1, function(x) rowSums(x*t(NMI_analysis_pca_water$rotation[,1])))
#NMI_analysis$water_pca2 <- apply(scale(NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")]),1, function(x) rowSums(x*t(NMI_analysis_pca_water$rotation[,2])))

#NMI_analysis_niche_pca_temp <- prcomp(scale(NMI_analysis[,c("niche.tmin","niche.tmax","niche.tmin_sd","niche.tmax_sd")]))
#pcatest <- PCAtest(scale(NMI_analysis[,c("niche.tmin","niche.tmax","niche.tmin_sd","niche.tmax_sd")]))
#NMI_analysis$niche_temp_pca1 <- apply(scale(NMI_analysis[,c("niche.tmin","niche.tmax","niche.tmin_sd","niche.tmax_sd")]),1, function(x) rowSums(x*t(NMI_analysis_niche_pca_temp$rotation[,1])))
#NMI_analysis$niche_temp_pca2 <- apply(scale(NMI_analysis[,c("niche.tmin","niche.tmax","niche.tmin_sd","niche.tmax_sd")]),1, function(x) rowSums(x*t(NMI_analysis_niche_pca_temp$rotation[,2])))

#NMI_analysis_niche_pca_water <- prcomp(scale(NMI_analysis[,c("niche.Prec","niche.Prec_sd","niche.SM","niche.SM_sd")]))
#pcatest <- PCAtest(scale(NMI_analysis[,c("niche.Prec","niche.Prec_sd","niche.SM","niche.SM_sd")]))
#NMI_analysis$niche_water_pca1 <- apply(scale(NMI_analysis[,c("niche.Prec","niche.Prec_sd","niche.SM","niche.SM_sd")]),1, function(x) rowSums(x*t(NMI_analysis_niche_pca_water$rotation[,1])))
#NMI_analysis$niche_water_pca2 <- apply(scale(NMI_analysis[,c("niche.Prec","niche.Prec_sd","niche.SM","niche.SM_sd")]),1, function(x) rowSums(x*t(NMI_analysis_niche_pca_water$rotation[,2])))

#NMI_analysis <- cbind(NMI_analysis,scaled=scale(NMI_analysis[,c("water_pca1","water_pca2","niche_water_pca1","niche_water_pca2","temp_pca1","temp_pca2","niche_temp_pca1","niche_temp_pca2","date")]))

#NMI_analysis$num <- ifelse(NMI_analysis$num >0, 1,0)

#####################Characterize environments
options(glmmTMB.cores=6)

library(glmmTMB)
library(lme4)

m1a <- glmmTMB(num~water_pca_log*mean_water_pca_log_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca_log+temp_pca||species)+(mean_water_pca_log_native+mean_temp_pca_native||polygon_name),
             data=NMI_analysis,family="binomial")
summary(m1a)

m1b <- glmmTMB(num~water_pca_sqrt*mean_water_pca_sqrt_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca_sqrt+temp_pca||species)+(mean_water_pca_sqrt_native+mean_temp_pca_native||polygon_name),
             data=NMI_analysis,family="binomial")
summary(m1b)

m1c <- glmmTMB(num~water_pca*mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca+temp_pca||species)+(mean_water_pca_native+mean_temp_pca_native||polygon_name),
             data=NMI_analysis,family="binomial")
summary(m1c)
AIC(m1a,m1b,m1c) #m1C the best

m1c <- glmmTMB(num~water_pca*mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca||species)+(mean_water_pca_native+mean_temp_pca_native||polygon_name),
               data=NMI_analysis,family="binomial")
summary(m1c) #removed one rando meffect as it is too weak and affect the calculation of R2
AIC(m1a,m1b,m1c) 
performance::r2(m1c)

m1c_simplified <- glmmTMB(num~water_pca+mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca||species)+(mean_water_pca_native+mean_temp_pca_native||polygon_name),
               data=NMI_analysis,family="binomial") #remove the weak random effect such that r2 can be calculated
summary(m1c_simplified)
performance::r2(m1c_simplified)
################simplified model
m1e <- glmmTMB(num~water_pca+mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca+temp_pca||species)+(mean_water_pca_native+mean_temp_pca_native||polygon_name),
               data=NMI_analysis,family="binomial")
summary(m1e)
m1e <- glmmTMB(num~water_pca+mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca||species)+(mean_water_pca_native+mean_temp_pca_native||polygon_name),
               data=NMI_analysis,family="binomial") #remove the weak random effect such that r2 can be calculated
summary(m1e)
r2(m1e)

m <- glmmTMB(num~temp_pca*mean_temp_pca_native+(temp_pca||species)+(mean_temp_pca_native||polygon_name),
             data=NMI_analysis,family="binomial")
summary(m)
r2(m)

#see coef, and if there is any sign of quasi-separation...for final model only
m_standardized <- glmmTMB(num~scale(temp_pca)*scale(mean_temp_pca_native)+(scale(temp_pca)||species)+(scale(mean_temp_pca_native)||polygon_name),
                          data=NMI_analysis,family="binomial")
summary(m_standardized)

write.csv(rbind(round(summary(m1c)$coef$cond,3),round(summary(m1c_simplified)$coef$cond,3),round(summary(m)$coef$cond,3)),"Indoor_ants/glmm.csv")

########################################## sensitivity testing (removed weak random effects such that performance::r2 runs - other results won't change)
#here date is scaled to faciliate convergence
region_name <- NMI_analysis_sum[NMI_analysis_sum$indoor > 0,"polygon_name"]
m1d1 <- glmmTMB(num~water_pca*mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca||species)+(mean_water_pca_native+mean_temp_pca_native||polygon_name),
               data=NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,],family="binomial") #subsetted analysis
nrow(NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,])
summary(m1d1)
performance::r2(m1d1)

region_name <- NMI_analysis_sum[NMI_analysis_sum$indoor > 1,"polygon_name"]
m1d2 <- glmmTMB(num~water_pca*mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(temp_pca+water_pca||species)+(0+mean_water_pca_native+mean_temp_pca_native||polygon_name),
               data=NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,],family="binomial") #subsetted analysis
nrow(NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,])
summary(m1d2)
performance::r2(m1d2)

region_name <- NMI_analysis_sum[NMI_analysis_sum$indoor > 2,"polygon_name"]
m1d3 <- glmmTMB(num~water_pca*mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(temp_pca||species)+(0+mean_water_pca_native+mean_temp_pca_native||polygon_name),
               data=NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,],family="binomial") #subsetted analysis
nrow(NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,])
summary(m1d3)
performance::r2(m1d3)

write.csv(rbind(round(summary(m1d1)$coef$cond,3),round(summary(m1d2)$coef$cond,3),round(summary(m1d3)$coef$cond,3)),"Indoor_ants/glmm_additional.csv")

#m1 <- bglmer(num~scaled.water_pca+scaled.niche_water_pca+scaled.temp_pca+scaled.niche_temp_pca+scaled.date+Generalist+Arboreal+(scaled.water_pca+scaled.temp_pca||species)+(scaled.niche_water_pca+scaled.niche_temp_pca||polygon_name),
#             data=NMI_analysis,family="binomial",fixef.prior = normal(cov = diag(9,8)),
#            control=glmerControl(optimizer = "bobyqa"))
#summary(m1)

#residuals(m1)

library(DHARMa)
library(geosphere)
m1 <- glmmTMB(num~temp_pca*mean_temp_pca_native+(temp_pca||species)+(mean_temp_pca_native||polygon_name),
             data=NMI_analysis[order(NMI_analysis$polygon_name),],family="binomial") #re-run model but make sure the factors are ordered alphabetically
summary(m1)

sres <- simulateResiduals(m1)
sres_spatial <- recalculateResiduals(sres,NMI_analysis[order(NMI_analysis$polygon_name),"polygon_name"])
polygon_unique <- NMI_analysis[!duplicated(NMI_analysis$polygon_name),c("polygon_name","x","y")] #make sure the structures are consistent between residual objects and distance M
polygon_unique <- polygon_unique[order(polygon_unique$polygon_name),]
polygon_dist <- distm(polygon_unique[,c("x","y")])
testSpatialAutocorrelation(sres_spatial,distMat=polygon_dist)

m2 <- glmmTMB(num~temp_pca*mean_temp_pca_native+(temp_pca||species)+(mean_temp_pca_native||polygon_name),
              data=NMI_analysis[order(NMI_analysis$species),],family="binomial")
summary(m2) #re-run model but make sure the factors are ordered alphabetically
sres <- simulateResiduals(m2)
phyD <- cophenetic(keep.tip(genus_tree,NMI_analysis$species))
sres_phylo <- recalculateResiduals(sres,NMI_analysis[order(NMI_analysis$species),"species"])
phyD <- phyD[order(rownames(phyD)),order(colnames(phyD))]
testSpatialAutocorrelation(sres_phylo,distMat=phyD)
###################################################Fig 2

library(ggeffects)
v1 <- seq(min(NMI_analysis$temp_pca_4C),max(NMI_analysis$temp_pca),length.out=100)
v2 <- seq(min(NMI_analysis$mean_temp_pca_native),max(NMI_analysis$mean_temp_pca_native),length.out=100)

predict_niche <- ggemmeans(m,terms=c("temp_pca [v1]","mean_temp_pca_native [v2]"),type="fixed",rg.limit=1000*1000)
predict_niche$group <- as.numeric(as.character(predict_niche$group))
NMI_analysis$Status <- ifelse(NMI_analysis$num == 0, "Indoor", "Naturalized")

theme <- theme(axis.line=element_line(colour="black"),
               axis.text = element_text(size=4.5),
               axis.title = element_text(size=4.5),
               legend.position="bottom",
               panel.background=element_rect(colour="white",fill="white"),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               plot.background=element_rect(colour="white",fill="white"),
               legend.key.height = unit(0.25, "cm"),
               legend.key.width = unit(0.3, "cm"),
               legend.title= element_text(size=4.5),
               legend.text = element_text(size=4.5),
               legend.margin=margin(t=0.2,unit="cm"))

p2<- ggplot(predict_niche,aes(y=group,x=x))+
  geom_raster(aes(fill=predicted))+
  geom_point(data=NMI_analysis,aes(y=mean_temp_pca_native,x=temp_pca,colour=Status),size=0.25)+
  #annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b)",size=3)+
  xlab("Temperature PCA1 of invaded region \n (Colder and more seasonal)")+
  ylab("Temperature PCA1 of native range \n (Colder and more seasonal)")+
  scale_fill_continuous(low="#ffffb2",high="#e31a1c",name="Predicted naturalization probability",breaks=c(0.1,0.5,0.9))+
  scale_color_manual(values=c("grey60","black"))+
  theme
plot(p2)

#library(ggeffects)
#v3 <- seq(min(NMI_analysis$water_pca),max(NMI_analysis$water_pca),length.out=100)

#predict_water <- ggemmeans(m,terms=c("temp_pca [v1]","water_pca[v3]"),type="fixed",rg.limit=100000)
#predict_water$group <- as.numeric(as.character(predict_water$group))

#prop <- NMI_analysis %>% group_by(polygon_name) %>% summarize(n = n(), sum = sum(num)) %>% mutate(prop = sum/n)
#NMI_analysis$prop <-  unlist(prop[match(NMI_analysis$polygon_name,prop$polygon_name),"prop"])

#p2a <- ggplot(predict_water,aes(y=group,x=x))+
  #geom_raster(aes(fill=predicted))+
  #geom_point(data=NMI_analysis,aes(y=water_pca,x=temp_pca,colour=prop),size=0.5)+
  #annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a)",size=3)+
  #xlab("Temperature PCA1 of invaded region \n (Colder and more seasonal)")+
  #ylab("Water PCA1 of invaded region \n (Wetter and more seasonal)")+
  #scale_fill_continuous(low="#ffffb2",high="#e31a1c",name="Predicted probability of establishing in the wild",breaks=c(0.1,0.5,0.9))+
  #scale_color_continuous(low="grey75",high="black",name="Proportion of species established in the wild",breaks=c(0.1,0.5,0.9))+
  #theme

#plot(p2a)

#library(ggpubr)
#p2 <- ggarrange(p2,common.legend=T,legend="bottom")
#plot(p2)
plot(p2)
ggsave("Indoor_ants/Fig2.tiff",dpi=800,height=8.4,width=10.4,units="cm",compression="lzw",bg="white")
###################################################
library(rnaturalearth)

capacity <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/capacity.csv")
country <- ne_countries(scale=10)
country <- vect(country)
country_df <- as.data.frame(country)

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

country_df$brk_name[country_df$brk_name == "Côte d'Ivoire"] <- "Cote d'Ivoire"
capacity$total <- capacity$threat+capacity$IAS_list+capacity$Existing_mgmt+capacity$research+capacity$outreach
country_df$total <- capacity[match(country_df$brk_name,capacity$Country_global),"total"]
country$total <- country_df$total
country_sf <- sf::st_as_sf(country)

country_sf <- country_sf[country_sf$brk_name != "Antarctica",]
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
               legend.key.height = unit(0.5, "cm"),
               legend.key.width = unit(2, "cm"),
               legend.text = element_text(size=12),
               legend.title=element_text(size=18),
               legend.margin=margin(t=-0.5,unit="cm"),
               legend.background=element_rect(fill="white"))

country_sf <- st_wrap_dateline(country_sf)
country_sf <- st_transform(country_sf,st_crs("ESRI:54019"))
p_Smap <- ggplot(data=country_sf,aes(fill=total))+
  geom_sf(colour="black",size=0.1)+
  labs(fill="Total response capacity score")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="grey")+
  theme

plot(p_Smap)

ggsave("Indoor_ants/p_Smap.tiff",dpi=800,width=16,height=10,compression="lzw")
######Lesser sunda islands as Indonesia, Aland as Finland, Malyesia & Singapore = Maylesia, 

not_found <- bentity.shp.df$Final_country[!bentity.shp.df$Final_country %in% capacity$Country_global]
not_found <- sort(unique(not_found))

not_found <- capacity$Country_global[!capacity$Country_global %in% bentity.shp.df$Final_country]
not_found <- sort(unique(not_found))

NMI_analysis$dummy <- 1
######################################################

#future_df <- NMI_analysis[,c("species","polygon_name","ID","sp.layer","scaled.date","scaled.niche_temp_pca","scaled.niche_water_pca","dummy")]

#source("Indoor_ants/Indoor_ants/rescale_future_climate.R")
#NMI_analysis$scaled.temp_pca_2C <- rescale_future_climate(NMI_analysis[,c("MAT.2C","Temp_sd.2C")],
                                                       #NMI_analysis[,c("MAT","Temp_sd")],
                                                       #NMI_analysis_pca_temp,
                                                       #NMI_analysis$temp_pca)
#NMI_analysis$scaled.temp_pca_4C <- rescale_future_climate(NMI_analysis[,c("MAT.4C","Temp_sd.4C")],
                                                       #NMI_analysis[,c("MAT","Temp_sd")],
                                                       #NMI_analysis_pca_temp,
                                                       #NMI_analysis$temp_pca)
#NMI_analysis$scaled.water_pca_2C <- rescale_future_climate(NMI_analysis[,c("Prec.2C","Prec_sd.2C","SM.2C","SM_sd.2C")],
                                                        #NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")],
                                                        #NMI_analysis_pca_water,
                                                        #NMI_analysis$water_pca)
#NMI_analysis$scaled.water_pca_4C<- rescale_future_climate(NMI_analysis[,c("Prec.4C","Prec_sd.4C","SM.4C","SM_sd.4C")],
                                                       #NMI_analysis[,c("Prec","Prec_sd","SM","SM_sd")],
                                                       #NMI_analysis_pca_water,
                                                       #NMI_analysis$water_pca)

#NMI_analysis$current_status_projection <- predict(m,type="response")

#predict_data <- list(NMI_analysis[,c("species","polygon_name","scaled.date","scaled.niche_temp_pca","scaled.niche_water_pca","scaled.temp_pca_2C","scaled.water_pca_2C","sp.layer")],
                    #NMI_analysis[,c("species","polygon_name","scaled.date","scaled.niche_temp_pca","scaled.niche_water_pca","scaled.temp_pca_4C","scaled.water_pca_4C","sp.layer")])

#predict_data <- lapply(predict_data,setNames,c("species","polygon_name","scaled.date","scaled.niche_temp_pca","scaled.niche_water_pca","scaled.temp_pca","scaled.water_pca","sp.layer"))

#future_occ <- do.call(cbind,lapply(predict_data,function(x) predict(m,x,type="response")))
#colnames(future_occ) <- c("future_status_2C","future_status_4C")
#NMI_analysis <- cbind(NMI_analysis,future_occ)

NMI_analysis$current_status_projection <- predict(m,type="response")

predict_df <- NMI_analysis[,c("temp_pca_2C","water_pca_2C","mean_water_pca_native","mean_temp_pca_native","sp.layer","date","species","polygon_name")]
colnames(predict_df)[1:2] <- c("temp_pca","water_pca")
NMI_analysis$future_status_2C <- predict(m,newdata=predict_df, type="response")

predict_df <- NMI_analysis[,c("temp_pca_4C","water_pca_4C","mean_water_pca_native","mean_temp_pca_native","sp.layer","date","species","polygon_name")]
colnames(predict_df)[1:2] <- c("temp_pca","water_pca")
NMI_analysis$future_status_4C <- predict(m,newdata=predict_df, type="response")

invasive <- c("Anoplolepis.custodiens","Anoplolepis.gracilipes","Azteca.sericeasur","Brachyponera.chinensis","Camponotus.conspicuus.zonatus","Cardiocondyla.wroughtonii","Cardiocondyla.emeryi",
              "Formica.aquilonia","Formica.paralugubris",
              "Lasius.neglectus","Linepithema.humile",
              "Monomorium.floricola","Monomorium.monomorium","Monomorium.pharaonis","Myrmica.rubra",
              "Nylanderia.bourbonica","Nylanderia.fulva", "Ochetellus.glaber",
              "Plagiolepis.alluaudi","Paratrechina.longicornis","Pheidole.megacephala","Pheidole.radoszkowskii","Solenopsis.invicta",
              "Solenopsis.geminata","Solenopsis.papuana","Solenopsis.richteri",
              "Tapinoma.melanocephalum","Tapinoma.sessile","Technomyrmex.albipes","Technomyrmex.difficilis","Tetramorium.bicarinatum","Tetramorium.simillimum",
              "Trichomyrmex.destructor","Wasmannia.auropunctata") #Lepisolota canescens & Nylanderia pubens excluded due to no impact record

NMI_analysis$invasive <- ifelse(NMI_analysis$species %in% invasive,"Invasive","Alien")

NMI_analysis$proj_diff_2C <- NMI_analysis$future_status_2C-NMI_analysis$current_status_projection
NMI_analysis$proj_diff_4C <- NMI_analysis$future_status_4C-NMI_analysis$current_status_projection
NMI_analysis$warming_diff <- NMI_analysis$future_status_4C-NMI_analysis$future_status_2C

######################################################
site_summary <- NMI_analysis %>% group_by(polygon_name,ID) %>% summarize(current_sum = sum(num),
                                                                         current_indoor = sum(dummy[Status == "Indoor"]),
                                                                         current_invasive = sum(num[invasive == "Invasive"]),
                                                                         current_indoor_invasive = sum(dummy[invasive == "Invasive" & Status == "Indoor"]),
                                                                         current_projection_sum = sum(current_status_projection),
                                                                         future_sum_2C=sum(future_status_2C),
                                                                         future_sum_4C=sum(future_status_4C),
                                                                         future_indoor_gain_2C = sum(future_status_2C[num==0]),
                                                                         future_indoor_gain_4C = sum(future_status_4C[num==0]),
                                                                         proj_diff_2C_net = sum(proj_diff_2C),
                                                                         proj_diff_4C_net = sum(proj_diff_4C),
                                                                         proj_diff_indoor_2C_net = sum(proj_diff_2C[num==0]),
                                                                         proj_diff_indoor_4C_net = sum(proj_diff_4C[num==0]),
                                                                         warming_diff_net = sum(warming_diff),
                                                                         warming_diff_indoor_net = sum(warming_diff[num==0]),
                                                                         current_invasive_sum = sum(num[invasive=="Invasive"]),
                                                                         current_projection_invasive_sum = sum(current_status_projection[invasive=="Invasive"]),
                                                                         future_sum_invasive_2C=sum(future_status_2C[invasive=="Invasive"]),
                                                                         future_sum_invasive_4C=sum(future_status_4C[invasive=="Invasive"]),
                                                                         future_invasive_indoor_gain_2C = sum(future_status_2C[invasive=="Invasive" & num==0]),
                                                                         future_invasive_indoor_gain_4C = sum(future_status_4C[invasive=="Invasive" & num==0]),
                                                                         proj_diff_invasive_2C_net = sum(proj_diff_2C[invasive=="Invasive"]),
                                                                         proj_diff_invasive_4C_net = sum(proj_diff_4C[invasive=="Invasive"]),
                                                                         proj_diff_invasive_indoor_2C_net = sum(proj_diff_2C[invasive=="Invasive" & num==0]),
                                                                         proj_diff_invasive_indoor_4C_net = sum(proj_diff_4C[invasive=="Invasive" & num==0]),
                                                                         warming_diff_invasive_net = sum(warming_diff[invasive=="Invasive"]),
                                                                         warming_diff_invasive_indoor_net = sum(warming_diff[num==0 & invasive=="Invasive"]),
                                                                         Percent_2C = proj_diff_indoor_2C_net/current_sum*100,
                                                                         Percent_4C = proj_diff_indoor_4C_net/current_sum*100,
                                                                         Percent_2C_invasive = proj_diff_invasive_indoor_2C_net/current_sum*100,
                                                                         Percent_4C_invasive = proj_diff_invasive_indoor_4C_net/current_sum*100)
empty_poly <- data.frame(polygon_name=bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary$polygon_name,"BENTITY2_N"],
           ID = rownames(bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary$polygon_name,]))

site_summary <- plyr::rbind.fill(site_summary,empty_poly)
site_summary[is.na(site_summary)] <- 0


######################################################

library(ggplot2)

clim_invasion_df$invasive <- clim_invasion_df$species %in% gsub("_",".",invasive)
indoor.sr <- subset(clim_invasion_df,num==0) %>% group_by(ID,polygon_name) %>% count(num)
indoor.sr.invasive <- subset(clim_invasion_df, num == 0 & invasive == T) %>% group_by(ID,polygon_name) %>% count(num)

#bentity.df$indoor.sr <- unlist(indoor.sr[match(bentity.df$geom,indoor.sr$ID),"n"])
#bentity.df$indoor.sr.invasive <- unlist(indoor.sr.invasive[match(bentity.df$geom,indoor.sr.invasive$ID),"n"])

outdoor.sr <- subset(clim_invasion_df,num==1) %>% group_by(ID,polygon_name) %>% count(num)
outdoor.sr.invasive <- subset(clim_invasion_df, num == 1 & invasive == T) %>% group_by(ID,polygon_name) %>% count(num)

#bentity.df$outdoor.sr <- unlist(outdoor.sr[match(bentity.df$geom,outdoor.sr$ID),"n"])
#bentity.df$outdoor.sr.invasive <- unlist(outdoor.sr.invasive[match(bentity.df$geom,outdoor.sr.invasive$ID),"n"])

#bentity.df$indoor.sr[bentity.df$indoor.sr == 0] <- NA
#bentity.df$outdoor.sr[bentity.df$outdoor.sr == 0] <- NA
#bentity.df$indoor.sr.invasive[bentity.df$indoor.sr.invasive == 0] <- NA
#bentity.df$outdoor.sr.invasive[bentity.df$outdoor.sr.invasive == 0] <- NA

bentity.shp.sf <- sf::st_as_sf(bentity.shp)
bentity.shp.sf$indoor.sr <- unlist(indoor.sr[match(bentity.shp$BENTITY2_N,indoor.sr$polygon_name),"n"])
bentity.shp.sf$indoor.sr.invasive <- unlist(indoor.sr.invasive[match(bentity.shp$BENTITY2_N,indoor.sr.invasive$polygon_name),"n"])
bentity.shp.sf$outdoor.sr <- unlist(outdoor.sr[match(bentity.shp$BENTITY2_N,outdoor.sr$polygon_name),"n"])
bentity.shp.sf$outdoor.sr.invasive <- unlist(outdoor.sr.invasive[match(bentity.shp$BENTITY2_N,outdoor.sr.invasive$polygon_name),"n"])

library(tiff)
library(grid)
indoor_tiff <- readTIFF("Indoor_ants/indoor.tif")
wild_tiff <- readTIFF("Indoor_ants/wild.tif")

g1 <- rasterGrob(indoor_tiff,width = unit(0.75,"cm"),height=unit(0.75,"cm"),interpolate=TRUE)
g2 <- rasterGrob(wild_tiff,width = unit(0.75,"cm"),height=unit(0.75,"cm"),interpolate=TRUE)

#centroid_df <- as.data.frame(centroid,geom="XY")
#centroid_df <- cbind(centroid_df,bentity.shp.sf[,c("indoor.sr","indoor.sr.invasive","outdoor.sr","outdoor.sr.invasive")])
#centroid_df <- centroid_df[,-9]
#centroid_df$area <- expanse(bentity.shp)

#centroid_df$num <- ifelse(str_detect(centroid_df$BENTITY2_N, "Island|Atoll|Archipelago|Isle|"),"Island","Non-island")
########https://rpubs.com/mdavril_gsu/794598

bentity.shp.sf <- st_wrap_dateline(bentity.shp.sf)
bentity.shp.sf <- st_transform(bentity.shp.sf, crs=st_crs("ESRI:54019"))
ratio <- tmaptools::get_asp_ratio(bentity.shp.sf)

xmax <- xmin <- -1.4e07
ymin <- ymax <- -5250000
size <- 2.7
theme <- theme(axis.line=element_blank(),
               axis.text.x=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               panel.background=element_rect(colour="white",fill="white"),
               panel.border=element_blank(),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               plot.background=element_rect(colour="white",fill="white"),
               plot.margin=unit(c(-0.1,0.1,-0.65,-0.1),"cm"),
               legend.key.height = unit(0.15, "cm"),
               legend.key.width = unit(0.25, "cm"),
               legend.text=element_text(size=4),
               legend.spacing.y = unit(0.5, "mm"),
               legend.margin=margin(t=-0.5,unit="cm"),
               legend.background=element_rect(fill="white"),
               legend.direction = "horizontal",
               legend.position = c(0.3, 0),
               legend.justification = c("right", "top"))

p1a <- ggplot(data=bentity.shp.sf,aes(fill=indoor.sr))+
  geom_sf(colour="black",size=0.1)+
  #geom_point(data=centroid_df,aes(x=x,y=y,fill=indoor.sr,shape=num,colour=num),size=2)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a) Alien",size=size)+
  annotation_custom(g1, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)+
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0,70),breaks=c(1,35,70))+
  #scale_colour_manual(values=c("black","transparent"))+
  #scale_shape_manual(values=c(21,1))+
  theme

plot(p1a)

p1b <- ggplot(data=bentity.shp.sf,aes(fill=outdoor.sr))+
  geom_sf(colour="black",size=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b) Alien",size=size)+
  annotation_custom(g2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0,70),breaks=c(1,35,70))+
  theme

p1c <- ggplot(data=bentity.shp.sf,aes(fill=indoor.sr.invasive))+
  geom_sf(colour="black",size=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(c) Harmful",size=size)+
  annotation_custom(g1, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0.5,12.5),breaks=c(1,6,12))+
  theme
plot(p1c)
p1d <- ggplot(data=bentity.shp.sf,aes(fill=outdoor.sr.invasive))+
  geom_sf(colour="black",size=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(d) Harmful",size=size)+
  annotation_custom(g2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0.5,12.5),breaks=c(1,6,12))+
  theme

#p1a <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=indoor.sr))+
  #geom_polygon(colour="black",size=0.1)+
  #annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a) Alien")+
  #annotation_custom(g1, xmin=-175, xmax=-175, ymin=-55, ymax=-55) +
  #labs(fill="")+
  #scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0,70),breaks=c(1,35,70))+
  #theme

#p1b <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=outdoor.sr))+
  #geom_polygon(colour="black",size=0.1)+
  #annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b) Alien")+
  #annotation_custom(g2, xmin=-175, xmax=-175, ymin=-55, ymax=-55) +
  #labs(fill="")+
  #scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",tran=scales::log1p_trans(),limits=c(0,70),breaks=c(1,35,70))+
  #theme

#p1c <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=indoor.sr.invasive))+
  #geom_polygon(colour="black",size=0.1)+
  #annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(c) Invasive")+
  #annotation_custom(g1, xmin=-175, xmax=-175, ymin=-55, ymax=-55) +
  #labs(fill="")+
  #scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0.5,12.5),breaks=c(1,6,12))+
  #theme

#p1d <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=outdoor.sr.invasive))+
  #geom_polygon(colour="black",size=0.1)+
  #annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(d) Invasive")+
  #annotation_custom(g2, xmin=-175, xmax=-175, ymin=-55, ymax=-55) +
  #labs(fill="")+
  #scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0.5,12.5),breaks=c(1,6,12))+
  #theme

library(ggpubr)
#p1_1 <- ggarrange(p1a,p1b,ncol=1,legend="none")
#plot(p1_1)
#p1_2 <- ggarrange(p1c,p1d,ncol=1,legend="none")
#plot(p1_2)

ratio <- tmaptools::get_asp_ratio(bentity.shp.sf)
p1 <- ggarrange(p1a,p1b,p1c,p1d,hjust=0,vjust=0,label.x=0,label.y=0)
plot(p1)
ggsave("Indoor_ants/Fig1.tiff",dpi=800,compression="lzw",units="cm",height=9.5,width=16.8,bg="white")
#####################################

bentity.shp.sf$proj_diff_indoor_2C_net <- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"proj_diff_indoor_2C_net"]
bentity.shp.sf$proj_diff_indoor_4C_net<- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"proj_diff_indoor_4C_net"]
bentity.shp.sf$proj_diff_invasive_indoor_2C_net <- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"proj_diff_invasive_indoor_2C_net"]
bentity.shp.sf$proj_diff_invasive_indoor_4C_net <- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"proj_diff_invasive_indoor_4C_net"]
bentity.shp.sf$warming_diff_indoor_net <- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"warming_diff_indoor_net"]
bentity.shp.sf$warming_diff_invasive_indoor_net<- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"warming_diff_invasive_indoor_net"]
bentity.shp.sf[is.na(bentity.shp.sf$indoor.sr),c("proj_diff_indoor_2C_net","proj_diff_indoor_4C_net","warming_diff_indoor_net")] <- NA
bentity.shp.sf[is.na(bentity.shp.sf$indoor.sr.invasive),c("proj_diff_invasive_indoor_2C_net","proj_diff_invasive_indoor_4C_net","warming_diff_invasive_indoor_net")] <- NA

library(tiff)
library(grid)
warming_tiff <- readTIFF("Indoor_ants/warming.tif")
g3 <- rasterGrob(warming_tiff,width = unit(0.7,"cm"),height=unit(0.7,"cm"),interpolate=TRUE)

xmax <- xmin <- -1.4e07
ymin <- ymax <- -5250000
space <- 900000
size <- 2.7
p3a <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_indoor_2C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a) Alien",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,4.8),breaks = c(0,2,4))+
  theme
plot(p3a)

p3b <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_indoor_4C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b) Alien",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,4.8),breaks = c(0,2,4))+
  theme
plot(p3b)

p3c <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_invasive_indoor_2C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(c) Harmful",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,1.8),breaks = c(0,0.5,1,1.5))+
  theme
plot(p3c)

p3d <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_invasive_indoor_4C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(d) Harmful",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,1.8),breaks = c(0,0.5,1,1.5))+
  theme
plot(p3d)

p3e <- ggplot(data=bentity.shp.sf,aes(fill=warming_diff_indoor_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a) Alien",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C vs 2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,3.5),breaks = c(0,1,2,3,4))+
  theme
plot(p3e)

p3f <- ggplot(data=bentity.shp.sf,aes(fill=warming_diff_invasive_indoor_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b) Harmful",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C vs 2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,1.3),breaks = c(0,0.5,1))+
  theme
plot(p3f)

library(ggpubr)
p3 <- ggarrange(p3a,p3c,p3b,p3d,p3e,p3f,nrow=3,ncol=2)
plot(p3)

ggsave("Indoor_ants/Fig3.tiff",dpi=800,compression="lzw",units="cm",height=13.5,width=16.8,bg="white")


# have to mean them first rather than calculating percent change in each region before averaging - because some regions have inf values (0 invasive currently at outdoor)

site_summary %>% filter(current_indoor >0) %>% summarise(percent_2C = mean(proj_diff_indoor_2C_net)/mean(current_sum)*100,
                                                       percent_4C = mean(proj_diff_indoor_4C_net)/mean(current_sum)*100,
                                                       gain_2C = mean(proj_diff_indoor_2C_net),
                                                       gain_4C = mean(proj_diff_indoor_4C_net))

site_summary %>% filter(current_indoor_invasive >0) %>% summarise(percent_2C_invasive = mean(proj_diff_invasive_indoor_2C_net)/mean(current_invasive)*100,
                                                         percent_4C_invasive = mean(proj_diff_invasive_indoor_4C_net)/mean(current_invasive)*100,
                                                         gain_2C_invasive = mean(proj_diff_invasive_indoor_2C_net),
                                                         gain_4C_invasive = mean(proj_diff_invasive_indoor_4C_net)) 

cor(site_summary$current_sum,site_summary$current_indoor,method="kendall")
cor(site_summary$proj_diff_indoor_2C_net,site_summary$proj_diff_indoor_4C_net,method="kendall")
cor(site_summary$future_sum_invasive_2C,site_summary$future_sum_invasive_4C,method="kendall")
