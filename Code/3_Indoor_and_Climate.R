source("Code/make_genus_tree.R")

library(phytools)
library(ape)
library(stringr)
library(tidyverse)
library(terra)
library(performance)
library(sf)
bentity.shp <- vect("Data/Bentity2_shapefile_fullres.shp")
clim_invasion_df <- read.csv("Data/clim_invasion_df.csv")

xy <- as.data.frame(centroids(bentity.shp),geom="XY")
clim_invasion_df <- cbind(clim_invasion_df,xy[match(clim_invasion_df$polygon_name,xy$BENTITY2_N),c("x","y")])
clim_invasion_df$num <- ifelse(clim_invasion_df$num ==0,0,1)

NMI_analysis <- subset(clim_invasion_df)
NMI_analysis <- na.omit(NMI_analysis)
NMI_analysis$species <- gsub(" ",".",NMI_analysis$species)
NMI_analysis_sum <- NMI_analysis %>% group_by(polygon_name) %>% dplyr::summarise(total = sum(num),
                                                                                 count = n(),
                                                                                 indoor = count-total)
ant_tree <- org_ant_tree <- read.tree("Data/backbone_MLtree_RaxML.tre")
is.ultrametric(ant_tree)
ant_tree <- as.phylo(hclust(as.dist(cophenetic(ant_tree)),method="average"))
dendextend::cor_cophenetic(ant_tree,org_ant_tree)
ant_tree$node.label <- 1:length(ant_tree$node.label)

ant_tree$tip.label <- gsub("[\"]","",ant_tree$tip.label)
ant_tree

ant_placement <- read.csv("Data/ant_genus_placement.csv")
genus_tree <- make_genus_tree(ant_tree,ant_placement,unique(NMI_analysis$species))
genus_tree$tip.label <-gsub("_",".",genus_tree$tip.label)

###GLMM indoor status & climate
options(glmmTMB.cores=6)

library(glmmTMB)

m1 <- glmmTMB(num~water_pca*mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca+temp_pca||species)+(mean_water_pca_native+mean_temp_pca_native||polygon_name),
             data=NMI_analysis,family="binomial") # model does converge. but note the very low random effect variance of temp_pca.
summary(m1)
performance::r2(m1) #singularity!

m1a <- glmmTMB(num~water_pca*mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca||species)+(mean_water_pca_native+mean_temp_pca_native||polygon_name),
               data=NMI_analysis,family="binomial") #converged. Note no changes in everything
car::Anova(m1a)
summary(m1a) #removed one random effect as it is too weak and affect the calculation of R2. 
performance::r2(m1a)

### remove interactions to retest main effect
m1_simplified <- glmmTMB(num~water_pca+mean_water_pca_native+temp_pca*mean_temp_pca_native+log10(date)+sp.layer+(water_pca||species)+(mean_water_pca_native+mean_temp_pca_native||polygon_name),
               data=NMI_analysis,family="binomial") 
summary(m1_simplified)
car::Anova(m1_simplified)
performance::r2(m1_simplified)

#final model
m <- glmmTMB(num~temp_pca*mean_temp_pca_native+(1|species)+(mean_temp_pca_native||polygon_name),
             data=NMI_analysis,family="binomial")
car::Anova(m)
summary(m)
performance::r2(m)

#see coef, and if there is any sign of quasi-separation...for final model only
m_standardized <- glmmTMB(num~scale(temp_pca)*scale(mean_temp_pca_native)+(scale(temp_pca)||species)+(scale(mean_temp_pca_native)||polygon_name),
                          data=NMI_analysis,family="binomial")
summary(m_standardized)

write.csv(rbind(round(car::Anova(m1a),3),round(car::Anova(m1_simplified),3),round(car::Anova(m),3)),"Results/glmm.csv")

### sensitivity testing (removed weak random effects such that performance::r2 runs - other results (e.g. p values) won't change)
### exclude reginos with zero indoor population
region_name <- NMI_analysis_sum[NMI_analysis_sum$indoor > 0,"polygon_name"] 
m1b <- glmmTMB(num~temp_pca*mean_temp_pca_native+(1|species)+(mean_temp_pca_native||polygon_name),
               data=NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,],family="binomial") #subsetted analysis
nrow(NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,])
car::Anova(m1b)
summary(m1b)
performance::r2(m1b)

### exclude regions with <= 1 indoor population
region_name <- NMI_analysis_sum[NMI_analysis_sum$indoor > 1,"polygon_name"]
m1c <- glmmTMB(num~temp_pca*mean_temp_pca_native+(1|species)+(mean_temp_pca_native||polygon_name),
               data=NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,],family="binomial") #subsetted analysis
nrow(NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,])
car::Anova(m1c)
summary(m1c)
performance::r2(m1c)

### exclude regions with <= 2 indoor population
region_name <- NMI_analysis_sum[NMI_analysis_sum$indoor > 2,"polygon_name"]
m1d <- glmmTMB(num~temp_pca*mean_temp_pca_native+(1|species)+(mean_temp_pca_native||polygon_name),
               data=NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,],family="binomial") #subsetted analysis
nrow(NMI_analysis[NMI_analysis$polygon_name %in% region_name$polygon_name,])
summary(m1d)
car::Anova(m1d)

performance::r2(m1d)

write.csv(rbind(round(car::Anova(m1b),3),
                round(car::Anova(m1c),3),
                round(car::Anova(m1d),3)),"Results/glmm_additional.csv")

### Check spatial autocorrelations
library(DHARMa)
library(geosphere)
m_best1 <- glmmTMB(num~temp_pca*mean_temp_pca_native+(1|species)+(mean_temp_pca_native||polygon_name),
             data=NMI_analysis[order(NMI_analysis$polygon_name),],family="binomial") #re-run model but make sure the factors are ordered alphabetically
summary(m_best1)

sres <- simulateResiduals(m_best1)
sres_spatial <- recalculateResiduals(sres,NMI_analysis[order(NMI_analysis$polygon_name),"polygon_name"])
polygon_unique <- NMI_analysis[!duplicated(NMI_analysis$polygon_name),c("polygon_name","x","y")] #make sure the structures are consistent between residual objects and distance M
polygon_unique <- polygon_unique[order(polygon_unique$polygon_name),]
polygon_dist <- distm(polygon_unique[,c("x","y")])
testSpatialAutocorrelation(sres_spatial,distMat=polygon_dist) #results significant, but very low observed values

### Check phylogenetic correlations
m_best2 <- glmmTMB(num~temp_pca*mean_temp_pca_native+(1|species)+(mean_temp_pca_native||polygon_name),
              data=NMI_analysis[order(NMI_analysis$species),],family="binomial")
summary(m_best2) #re-run model but make sure the factors are ordered alphabetically
sres <- simulateResiduals(m_best2)
phyD <- cophenetic(keep.tip(genus_tree,NMI_analysis$species))
sres_phylo <- recalculateResiduals(sres,NMI_analysis[order(NMI_analysis$species),"species"])
phyD <- phyD[order(rownames(phyD)),order(colnames(phyD))]
testSpatialAutocorrelation(sres_phylo,distMat=phyD) #results significant, but very low observed values

### make Fig 1 (Best model predictions)
library(ggeffects)
v1 <- seq(min(NMI_analysis$temp_pca_4C),max(NMI_analysis$temp_pca),length.out=100)
v2 <- seq(min(NMI_analysis$mean_temp_pca_native),max(NMI_analysis$mean_temp_pca_native),length.out=100)

predict_niche <- ggemmeans(m,terms=c("temp_pca [v1]","mean_temp_pca_native [v2]"),type="fixed",rg.limit=1000*1000)
predict_niche$group <- as.numeric(as.character(predict_niche$group))
NMI_analysis$Status <- ifelse(NMI_analysis$num == 0, "Indoor", "Naturalized")

theme <- theme(axis.line=element_line(colour="black"),
               axis.text = element_text(size=7),
               axis.title = element_text(size=7),
               legend.position="bottom",
               panel.background=element_rect(colour="white",fill="white"),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               plot.background=element_rect(colour="white",fill="white"),
               legend.key.height = unit(0.25, "cm"),
               legend.key.width = unit(0.6, "cm"),
               legend.title= element_text(size=7),
               legend.text = element_text(size=7),
               legend.margin=margin(t=0.2,unit="cm"))

p1<- ggplot(predict_niche,aes(y=group,x=x))+
  geom_raster(aes(fill=predicted))+
  geom_point(data=NMI_analysis,aes(y=mean_temp_pca_native,x=temp_pca,colour=Status),size=0.25)+
  #annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b)",size=3)+
  xlab("Temperature PCA1 of invaded region \n (Colder and more seasonal)")+
  ylab("Temperature PCA1 of native range \n (Colder and more seasonal)")+
  scale_fill_continuous(low="#ffffb2",high="#e31a1c",name="Predicted naturalization probability",limits=c(0,1),breaks=c(0,0.25,0.5,0.75,1))+
  scale_color_manual(values=c("grey60","black"),name="")+
  theme+
  guides(fill=guide_colourbar(title.position="top"),
        colour=guide_legend(title.position="top",
                            override.aes = list(size = 5)))
plot(p1)

ggsave("Figures/Model_results.tiff",dpi=800,height=8.4,width=10.4,units="cm",compression="lzw",bg="white")

### Climate change projection
NMI_analysis$dummy <- 1
NMI_analysis$current_status_projection <- predict(m,type="response")
bentity.shp.df <- as.data.frame(bentity.shp)

predict_df <- NMI_analysis[,c("temp_pca_2C","water_pca_2C","mean_water_pca_native","mean_temp_pca_native","sp.layer","date","species","polygon_name")]
colnames(predict_df)[1:2] <- c("temp_pca","water_pca")
NMI_analysis$future_status_2C <- predict(m,newdata=predict_df, type="response")

predict_df <- NMI_analysis[,c("temp_pca_4C","water_pca_4C","mean_water_pca_native","mean_temp_pca_native","sp.layer","date","species","polygon_name")]
colnames(predict_df)[1:2] <- c("temp_pca","water_pca")
NMI_analysis$future_status_4C <- predict(m,newdata=predict_df, type="response")

#############
backup_NMI_analysis <- NMI_analysis
GISS <- read.csv("Data/GISS.csv")
GISS <- GISS[GISS$Health_type != "Indoor (Occurrence)" & GISS$Health_type != "Indoor (Pathogen)" & GISS$Health_type != "Outdoor (Occurrence)" & GISS$Health_type != "Outdoor (Pathogen)",]
library(tidyverse)

GISS[GISS == "NI"] <- NA
GISS_long <- GISS %>% 
  select(Species.name,Plants_score,Animals_score,Competition_score,Ecosystems_score,Diseases_score,Hybridization_score,Crops_score,Animal.production_score,Forestry_score,Infrastructure_score,Health_score,Social_score,
         Plants_Confidence_level,Animals_Confidence_level,Competition_Confidence_level,Ecosystems_Confidence_level,Diseases_Confidence_level,Hybridization_Confidence_level,Crops_Confidence_level,Animal.production_Confidence_level,Forestry_Confidence_level,Infrastructure_Confidence_level,Health_Confidence_level,Social_Confidence_level) %>%
  pivot_longer(!Species.name,names_to=c("Impact",".value"),names_sep="_") %>%
  mutate(score_confidence = replace(score,Confidence < 2, 0)) %>%
  group_by(Species.name,Impact) %>%
  summarize(max.score = max(score_confidence,na.rm=T),
            min.score = min(score_confidence,na.rm=T)) #Inf / -Inf = empty = 0

GISS_long[GISS_long == "-Inf"] <- 0
GISS_long[GISS_long == "Inf"] <- 0

Impact_score <- GISS_long %>% select(-min.score) %>% pivot_wider(names_from=Impact,values_from=max.score)
Impact_score$E.Total <- rowSums(Impact_score[,c("Plants","Animals","Competition","Hybridization","Diseases","Ecosystems")],na.rm=T)
Impact_score$S.Total <- rowSums(Impact_score[,c("Crops","Health","Forestry","Infrastructure","Social","Health")],na.rm=T)

Harmful <- Impact_score %>%
  filter(E.Total > 0 | S.Total > 0) %>%
  select(Species.name,E.Total,S.Total)
###

NMI_analysis$Harmful <- ifelse(gsub("\\."," ",NMI_analysis$species) %in% Harmful$Species.name,"Harmful","Alien")
NMI_analysis$species[NMI_analysis$Harmful == "Harmful" & NMI_analysis$num == 0]
unique(NMI_analysis$species[NMI_analysis$Harmful == "Harmful" & NMI_analysis$num == 0])

NMI_analysis$proj_diff_2C <- NMI_analysis$future_status_2C-NMI_analysis$current_status_projection
NMI_analysis$proj_diff_4C <- NMI_analysis$future_status_4C-NMI_analysis$current_status_projection
NMI_analysis$warming_diff <- NMI_analysis$future_status_4C-NMI_analysis$future_status_2C

NMI_analysis %>% 
  filter(num == 0) %>% 
  summarize(mean(proj_diff_2C),mean(proj_diff_4C),min(proj_diff_2C),max(proj_diff_2C),min(proj_diff_4C),max(proj_diff_4C)) #summary statisticsf for probabilities at population levels
### summary statistics per polygon
site_summary <- NMI_analysis %>% group_by(polygon_name,ID) %>% summarize(current_sum = sum(num),
                                                                         current_indoor = sum(dummy[Status == "Indoor"]),
                                                                         current_Harmful = sum(num[Harmful == "Harmful"]),
                                                                         current_indoor_Harmful = sum(dummy[Harmful == "Harmful" & Status == "Indoor"]),
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
                                                                         current_Harmful_sum = sum(num[Harmful=="Harmful"]),
                                                                         current_projection_Harmful_sum = sum(current_status_projection[Harmful=="Harmful"]),
                                                                         future_sum_Harmful_2C=sum(future_status_2C[Harmful=="Harmful"]),
                                                                         future_sum_Harmful_4C=sum(future_status_4C[Harmful=="Harmful"]),
                                                                         future_Harmful_indoor_gain_2C = sum(future_status_2C[Harmful=="Harmful" & num==0]),
                                                                         future_Harmful_indoor_gain_4C = sum(future_status_4C[Harmful=="Harmful" & num==0]),
                                                                         proj_diff_Harmful_2C_net = sum(proj_diff_2C[Harmful=="Harmful"]),
                                                                         proj_diff_Harmful_4C_net = sum(proj_diff_4C[Harmful=="Harmful"]),
                                                                         proj_diff_Harmful_indoor_2C_net = sum(proj_diff_2C[Harmful=="Harmful" & num==0]),
                                                                         proj_diff_Harmful_indoor_4C_net = sum(proj_diff_4C[Harmful=="Harmful" & num==0]),
                                                                         warming_diff_Harmful_net = sum(warming_diff[Harmful=="Harmful"]),
                                                                         warming_diff_Harmful_indoor_net = sum(warming_diff[num==0 & Harmful=="Harmful"]),
                                                                         Percent_2C = proj_diff_indoor_2C_net/current_sum*100,
                                                                         Percent_4C = proj_diff_indoor_4C_net/current_sum*100,
                                                                         Percent_2C_Harmful = proj_diff_Harmful_indoor_2C_net/current_sum*100,
                                                                         Percent_4C_Harmful = proj_diff_Harmful_indoor_4C_net/current_sum*100)
empty_poly <- data.frame(polygon_name=bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary$polygon_name,"BENTITY2_N"],
           ID = rownames(bentity.shp.df[!bentity.shp.df$BENTITY2_N %in% site_summary$polygon_name,]))

site_summary <- plyr::rbind.fill(site_summary,empty_poly)
site_summary[is.na(site_summary)] <- 0

### Global average gain (%) in alien / harmful species richness
### have to mean them first rather than calculating percent change in each region before averaging - because some regions have inf values (0 harmful species currently at outdoor)

site_summary %>% 
  filter(current_indoor >0) %>% 
  summarise(percent_2C = mean(proj_diff_indoor_2C_net)/mean(current_sum)*100,
            percent_4C = mean(proj_diff_indoor_4C_net)/mean(current_sum)*100,
            gain_2C = mean(proj_diff_indoor_2C_net),
            gain_4C = mean(proj_diff_indoor_4C_net),
            gain_2C_min = min(proj_diff_indoor_2C_net),
            gain_4C_min = min(proj_diff_indoor_4C_net),
            gain_2C_max = max(proj_diff_indoor_2C_net),
            gain_4C_max = max(proj_diff_indoor_4C_net))

site_summary %>% 
  filter(current_indoor_Harmful >0) %>% 
  summarise(percent_2C_Harmful = mean(proj_diff_Harmful_indoor_2C_net)/mean(current_Harmful)*100,
            percent_4C_Harmful = mean(proj_diff_Harmful_indoor_4C_net)/mean(current_Harmful)*100,
            gain_2C_Harmful = mean(proj_diff_Harmful_indoor_2C_net),
            gain_4C_Harmful = mean(proj_diff_Harmful_indoor_4C_net),
            gain_2C_Harmful_min = min(proj_diff_Harmful_indoor_2C_net),
            gain_4C_Harmful_min = min(proj_diff_Harmful_indoor_4C_net),
            gain_2C_Harmful_max = max(proj_diff_Harmful_indoor_2C_net),
            gain_4C_Harmful_max = max(proj_diff_Harmful_indoor_4C_net)) 

cor(site_summary$current_sum,site_summary$current_indoor,method="kendall") #correlation between indoor & outdoor alien richness based on current climate (observed)

###Figure S1a-d (Species richness)
library(ggplot2)

indoor.sr <- subset(NMI_analysis,num==0) %>% group_by(ID,polygon_name) %>% count(num)
indoor.sr.Harmful <- subset(NMI_analysis, num == 0 & Harmful == "Harmful") %>% group_by(ID,polygon_name) %>% count(num)

outdoor.sr <- subset(NMI_analysis,num==1) %>% group_by(ID,polygon_name) %>% count(num)
outdoor.sr.Harmful <- subset(NMI_analysis, num == 1 & Harmful =="Harmful") %>% group_by(ID,polygon_name) %>% count(num)

bentity.shp.sf <- sf::st_as_sf(bentity.shp)
bentity.shp.sf$indoor.sr <- unlist(indoor.sr[match(bentity.shp$BENTITY2_N,indoor.sr$polygon_name),"n"])
bentity.shp.sf$indoor.sr.Harmful <- unlist(indoor.sr.Harmful[match(bentity.shp$BENTITY2_N,indoor.sr.Harmful$polygon_name),"n"])
bentity.shp.sf$outdoor.sr <- unlist(outdoor.sr[match(bentity.shp$BENTITY2_N,outdoor.sr$polygon_name),"n"])
bentity.shp.sf$outdoor.sr.Harmful <- unlist(outdoor.sr.Harmful[match(bentity.shp$BENTITY2_N,outdoor.sr.Harmful$polygon_name),"n"])

library(tiff)
library(grid)
indoor_tiff <- readTIFF("Figures/indoor.tif")
wild_tiff <- readTIFF("Figures/wild.tif")

g1 <- rasterGrob(indoor_tiff,width = unit(0.75,"cm"),height=unit(0.75,"cm"),interpolate=TRUE)
g2 <- rasterGrob(wild_tiff,width = unit(0.75,"cm"),height=unit(0.75,"cm"),interpolate=TRUE)

bentity.shp.sf <- st_wrap_dateline(bentity.shp.sf)
bentity.shp.sf <- st_transform(bentity.shp.sf, crs=st_crs("ESRI:54019"))
ratio <- tmaptools::get_asp_ratio(bentity.shp.sf)

xmax <- xmin <- -1.4e07
ymin <- ymax <- -5250000
size <- 2.85
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
               plot.margin=unit(c(-0.1,0.1,-0.5,-0.1),"cm"),
               legend.key.height = unit(0.15, "cm"),
               legend.key.width = unit(0.25, "cm"),
               legend.text=element_text(size=7),
               legend.spacing.y = unit(0.5, "mm"),
               legend.margin=margin(t=-0.5,unit="cm"),
               legend.background=element_rect(fill="white"),
               legend.direction = "horizontal",
               legend.position = c(0.3, 0),
               legend.justification = c("right", "top"))

p1a <- ggplot(data=bentity.shp.sf,aes(fill=indoor.sr))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a) Non-native",size=size)+
  annotation_custom(g1, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)+
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(1,70),breaks=c(1,24,47,70))+
  theme

plot(p1a)

p1b <- ggplot(data=bentity.shp.sf,aes(fill=outdoor.sr))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b) Non-native",size=size)+
  annotation_custom(g2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(1,70),breaks=c(1,24,47,70))+
  theme

p1c <- ggplot(data=bentity.shp.sf,aes(fill=indoor.sr.Harmful))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(c) Harmful",size=size)+
  annotation_custom(g1, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(1,8),breaks=c(1,3,5,7))+
  theme
plot(p1c)

p1d <- ggplot(data=bentity.shp.sf,aes(fill=outdoor.sr.Harmful))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(d) Harmful",size=size)+
  annotation_custom(g2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(1,16),breaks=c(1,6,11,16))+
  theme

###Fig2a-d
bentity.shp.sf$proj_diff_indoor_2C_net <- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"proj_diff_indoor_2C_net"]
bentity.shp.sf$proj_diff_indoor_4C_net<- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"proj_diff_indoor_4C_net"]
bentity.shp.sf$proj_diff_Harmful_indoor_2C_net <- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"proj_diff_Harmful_indoor_2C_net"]
bentity.shp.sf$proj_diff_Harmful_indoor_4C_net <- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"proj_diff_Harmful_indoor_4C_net"]
bentity.shp.sf$warming_diff_indoor_net <- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"warming_diff_indoor_net"]
bentity.shp.sf$warming_diff_Harmful_indoor_net<- site_summary[match(bentity.shp.sf$BENTITY2_N,site_summary$polygon_name),"warming_diff_Harmful_indoor_net"]
bentity.shp.sf[is.na(bentity.shp.sf$indoor.sr),c("proj_diff_indoor_2C_net","proj_diff_indoor_4C_net","warming_diff_indoor_net")] <- NA
bentity.shp.sf[is.na(bentity.shp.sf$indoor.sr.Harmful),c("proj_diff_Harmful_indoor_2C_net","proj_diff_Harmful_indoor_4C_net","warming_diff_Harmful_indoor_net")] <- NA

library(tiff)
library(grid)
warming_tiff <- readTIFF("Figures/warming.tif")
g3 <- rasterGrob(warming_tiff,width = unit(0.7,"cm"),height=unit(0.7,"cm"),interpolate=TRUE)

xmax <- xmin <- -1.4e07
ymin <- ymax <- -5250000
space <- 900000
size <- 2.85
p2a <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_indoor_2C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a) Non-native",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0,5.3),breaks = c(0,2.5,5))+
  theme
plot(p2a)

p2b <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_indoor_4C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b) Non-native",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0,5.3),breaks = c(0,2.5,5))+
  theme
plot(p2b)

p2c <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_Harmful_indoor_2C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(c) Harmful",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0,1.5),breaks = c(0,0.5,1,1.5))+
  theme
plot(p2c)

p2d <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_Harmful_indoor_4C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(d) Harmful",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0,1.5),breaks = c(0,0.5,1,1.5))+
  theme
plot(p2d)

###Fig S4a-b
pS4a <- ggplot(data=bentity.shp.sf,aes(fill=warming_diff_indoor_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(a) Non-native",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C vs 2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0,4),breaks = c(0,2,4))+
  theme
plot(pS4a)

pS4b <- ggplot(data=bentity.shp.sf,aes(fill=warming_diff_Harmful_indoor_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(b) Harmful",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C vs 2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(0,1),breaks = c(0,0.5,1))+
  theme
plot(pS4b)

