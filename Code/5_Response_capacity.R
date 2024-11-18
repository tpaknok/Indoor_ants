### Response capacity (country-level) analyses
### Homogenizing the two datasets (country vs region)
library(rnaturalearth)

capacity <- read.csv("Data/capacity.csv") #only differences between the source and our file is that we added a column "Country_global_modified"
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

#Aggregate bentity systems according to the capacity data
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

#slight change here such that different naming between system matches
country_df$brk_name[country_df$brk_name == "Côte d'Ivoire"] <- "Cote d'Ivoire"

#Make Fig S6, and prepare data for other graphs
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
               legend.key.height = unit(0.25, "cm"),
               legend.key.width = unit(1, "cm"),
               legend.text = element_text(size=12),
               legend.title=element_text(size=12),
               legend.margin=margin(t=0,unit="cm"),
               legend.background=element_rect(fill="white"))

country_sf <- st_wrap_dateline(country_sf)
country_sf <- st_transform(country_sf,st_crs("ESRI:54019"))

pS6 <- ggplot(data=country_sf,aes(fill=total))+
  geom_sf(colour="black",size=0.1)+
  labs(fill="Total response capacity score")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white")+
  theme

plot(pS6)

ggsave("Figures/Response_capacity.tiff",dpi=800,width=12,height=8,compression="lzw",units="cm")

not_found <- bentity.shp.df$Final_country[!bentity.shp.df$Final_country %in% capacity$Country_global]
not_found <- sort(unique(not_found)) #Bentity polygon not considered

not_found <- capacity$Country_global[!capacity$Country_global %in% bentity.shp.df$Final_country]
not_found <- sort(unique(not_found)) #capacity data not considered (Note that many of these polygons have no data anyway)

### Country-level analyses
NMI_analysis$Final_country <- bentity.shp.df$Final_country[match(NMI_analysis$polygon_name,bentity.shp.df$BENTITY2_N)]
NMI_analysis <- cbind(NMI_analysis,capacity[match(NMI_analysis$Final_country,capacity$Country_global),c("threat","IAS_list","Existing_mgmt","research","outreach")])

NMI_analysis$total_score <- apply(NMI_analysis[,c("threat","IAS_list","Existing_mgmt","research","outreach")],1,sum)

length(which(is.na(NMI_analysis$total_score)))

sp_country <- NMI_analysis %>% group_by(species,Final_country,Harmful,total_score) %>% filter(num == 0) %>% summarise(proj_diff_2C = max(proj_diff_2C),
                                                                                                                       proj_diff_4C = max(proj_diff_4C),
                                                                                                                       Impact_2C.E.Total = max(Impact_2C.E.Total),
                                                                                                                       Impact_4C.E.Total = max(Impact_4C.E.Total),
                                                                                                                       Impact_2C.S.Total = max(Impact_2C.S.Total),
                                                                                                                       Impact_4C.S.Total = max(Impact_4C.S.Total))

### Distribution of climate change effects along different response capacities
record_score <- sp_country %>% group_by(total_score) %>% summarise("2°C Non-native ant richness" = sum(proj_diff_2C),
                                                                   "4°C Non-native ant richness" = sum(proj_diff_4C),
                                                                   "2°C Harmful ant richness" = sum(proj_diff_2C[Harmful == "Harmful"]),
                                                                   "4°C Harmful ant richness" = sum(proj_diff_4C[Harmful=="Harmful"]),
                                                                   "2°C Environmental impact" = sum(Impact_2C.E.Total),
                                                                   "4°C Environmental impact" = sum(Impact_4C.E.Total),
                                                                   "2°C Socioeconomic impact" = sum(Impact_2C.S.Total),
                                                                   "4°C Socioeconomic impact" = sum(Impact_4C.S.Total))

record_score <- record_score %>% pivot_longer(!total_score)
na.omit(record_score) %>% group_by(name) %>% summarize(percent = value[total_score==5]/sum(value)*100) # impacts in countries with tl 5 scores

record_score$name <- factor(record_score$name,levels=c("2°C Non-native ant richness",
                                                       "4°C Non-native ant richness", 
                                                       "2°C Harmful ant richness",
                                                       "4°C Harmful ant richness",
                                                       "2°C Environmental impact",
                                                       "4°C Environmental impact",
                                                       "2°C Socioeconomic impact",
                                                       "4°C Socioeconomic impact"))

### Make Fig3
p3 <- ggplot(data=record_score)+
  geom_bar(aes(x=total_score,weight=value))+
  xlab("Total response capcity score (0-5)")+
  ylab("Increases driven by climate change")+
  facet_wrap(~name,scales="free",nrow=4,ncol=2)+
  theme_classic()+
  scale_x_continuous(breaks=c(0,1,2,3,4,5))+
  theme(legend.position="bw")

plot(p3)
ggsave("Figures/Score_and_impact.tiff",dpi=800,width=16.8,height=16.8,units="cm",compression="lzw",bg="white")

### Make FigS5 (Country level variations - given the same response capacities)
country_record_score <- sp_country %>% group_by(total_score,Final_country) %>% summarise("2°C Non-native ant richness" = sum(proj_diff_2C),
                                                                           "4°C Non-native ant richness" = sum(proj_diff_4C),
                                                                           "2°C Harmful ant richness" = sum(proj_diff_2C[Harmful == "Harmful"]),
                                                                           "4°C Harmful ant richness" = sum(proj_diff_4C[Harmful=="Harmful"]),
                                                                           "2°C Environmental impact" = sum(Impact_2C.E.Total),
                                                                           "4°C Environmental impact" = sum(Impact_4C.E.Total),
                                                                           "2°C Socioeconomic impact" = sum(Impact_2C.S.Total),
                                                                           "4°C Socioeconomic impact" = sum(Impact_4C.S.Total))

sum_stat <- country_record_score %>% group_by(total_score) %>% summarise(n= n(),
                                                                   mean_2C_alien = mean(`2°C Non-native`),
                                                                   sd_2C_alien = sd(`2°C Non-native`),
                                                                   mean_4C_alien = mean(`4°C Non-native`),
                                                                   sd_4C_alien = sd(`4°C Non-native`),
                                                                   mean_2C_harmful = mean(`2°C Harmful`),
                                                                   sd_2C_harmful = sd(`2°C Harmful`),
                                                                   mean_4C_harmful = mean(`4°C Harmful`),
                                                                   sd_4C_harmful = sd(`4°C Harmful`),
                                                                   mean_2C_Environmental = mean(`2°C Environmental`),
                                                                   sd_2C_Environmental = sd(`2°C Environmental`),
                                                                   mean_2C_Socioeconomic = mean(`2°C Socioeconomic`),
                                                                   sd_2C_Socioeconomic = sd(`2°C Socioeconomic`),
                                                                   mean_4C_Environmental = mean(`4°C Environmental`),
                                                                   sd_4C_Environmental = sd(`4°C Environmental`),
                                                                   mean_4C_Socioeconomic = mean(`4°C Socioeconomic`),
                                                                   sd_4C_Socioeconomic = sd(`4°C Socioeconomic`))
                                                                   


country_record_score_long <- na.omit(country_record_score) %>% pivot_longer(cols=-c(Final_country,total_score),names_to="scenario")
country_record_score_long$scenario <- factor(country_record_score_long$scenario,levels=c("2°C Non-native ant richness",
                                                                                         "4°C Non-native ant richness", 
                                                                                         "2°C Harmful ant richness",
                                                                                         "4°C Harmful ant richness",
                                                                                         "2°C Environmental impact",
                                                                                         "4°C Environmental impact",
                                                                                         "2°C Socioeconomic impact",
                                                                                         "4°C Socioeconomic impact"))

pS5 <- ggplot(data=country_record_score_long)+
  geom_boxplot(aes(x=total_score,y=value,group=total_score))+
  #geom_jitter(aes(x=total_score,y=value),width=0.1)+
  xlab("Total response capacity score (0-5)")+
  ylab("Increases driven by climate change")+
  facet_wrap(~scenario,scales="free",nrow=4,ncol=2)+
  theme_classic()+
  theme(legend.position="bw")

plot(pS5)
ggsave("Figures/Score_and_impact_boxplot.tiff",dpi=800,width=16.8,height=16.8,units="cm",compression="lzw",bg="white")

### figure 4
country_score <- sp_country %>% group_by(Final_country) %>% summarise("2°C Non-native" = sum(proj_diff_2C),
                                                                      "4°C Non-native" = sum(proj_diff_4C),
                                                                      "2°C Harmful" = sum(proj_diff_2C[Harmful == "Harmful"]),
                                                                      "4°C Harmful" = sum(proj_diff_4C[Harmful=="Harmful"]),
                                                                      "2°C Environmental" = sum(Impact_2C.E.Total),
                                                                      "4°C Environmental" = sum(Impact_4C.E.Total),
                                                                      "2°C Socioeconomic" = sum(Impact_2C.S.Total),
                                                                      "4°C Socioeconomic" = sum(Impact_4C.S.Total))

country_score <- cbind(country_score,capacity[match(country_score$Final_country,capacity$Country_global),c("threat","IAS_list","Existing_mgmt","research","outreach")])
var = c("2°C Non-native","4°C Non-native","2°C Harmful","4°C Harmful","2°C Environmental","4°C Environmental","2°C Socioeconomic","4°C Socioeconomic")
top_n <- 20
country_score <- na.omit(country_score)
country_score$Final_country[country_score$Final_country == "Dem. Rep. Korea"] <- "North Korea" #change name
country_score$Final_country[country_score$Final_country == "Republic of Korea"] <- "South Korea" #change name

prop_df <- NULL
for (num in 1:8){
  top_20_pos <- order(country_score[,var[[num]]],decreasing=T)[1:top_n]
  plot_df <- country_score[top_20_pos,]
  
  top_10 <- c(var=var[[num]],top="top 10",colSums(ifelse(plot_df[1:10,c("threat","IAS_list","Existing_mgmt","research","outreach")]==1,1,0)/10*100))
  top_20 <- c(var=var[[num]],top="top 20",colSums(ifelse(plot_df[1:20,c("threat","IAS_list","Existing_mgmt","research","outreach")]==1,1,0)/20*100))
  prop_df <- rbind(prop_df,top_10,top_20)
  plot_df$Final_country <- factor(plot_df$Final_country,plot_df$Final_country[order(plot_df[,var[[num]]])])
  
  plot_df <- plot_df %>% select(Final_country,threat,IAS_list,Existing_mgmt,research,outreach) %>% pivot_longer(!Final_country)
  plot_df$value <- as.factor(plot_df$value)
  plot_df$name <- as.factor(plot_df$name)

  plot_df$name <- recode_factor(plot_df$name, Existing_mgmt  = "Management", 
                               IAS_list = "List",
                               outreach = "Monitoring",
                               research = "Research",
                               threat = "Threat")
  
  assign(paste0("p_cap_",num),ggplot(plot_df,aes(y=Final_country,x=name,fill=value))+
                                       geom_tile(colour="black")+
                                       ylab("")+
                                       xlab("")+
                                       labs(fill="Score")+
                                       scale_fill_manual(values=c("#bd0026","#fd8d3c","#ffffb2"))+
                                       ggtitle("")+
                                       theme_classic()+
                                       theme(axis.text.x = if(num==1){element_text(angle = 10,size=7)}else{element_text(angle=10,colour = "white",size=8)},
                                             axis.text.y = element_text(size=5.5)))
  }

p4 <- ggarrange(p_cap_1,p_cap_2,p_cap_3,p_cap_4,p_cap_5,p_cap_6,p_cap_7,p_cap_8,
                labels = c("2°C Non-native ant richness",
                           "4°C Non-native ant richness",
                           "2°C Harmful ant richness",
                           "4°C Harmful ant richness",
                           "2°C Environmental impact",
                           "4°C Environmental impact",
                           "2°C Socioeconomic impact",
                           "4°C Socioeconomic impact"),
                hjust=-0.2,
                vjust=0,
                label.y=0.9,
                font.label=list(size=12),
                nrow=4,ncol=2,common.legend=T,legend="bottom")
plot(p4)

ggsave("Figures/Fig4.tiff",dpi=800,width=16.8,height=24,units="cm",compression="lzw",bg="white")
write.csv(prop_df,"Results/prop.csv")

### Summary statistics Table S5 (% of climate change effects on the top-10 / 20 countries)
sum(country_score$`2°C Non-native`[order(country_score$`2°C Non-native`,decreasing=T)][1:10])/sum(country_score$`2°C Non-native`)
sum(country_score$`2°C Non-native`[order(country_score$`2°C Non-native`,decreasing=T)][1:20])/sum(country_score$`2°C Non-native`)
sum(country_score$`4°C Non-native`[order(country_score$`4°C Non-native`,decreasing=T)][1:10])/sum(country_score$`4°C Non-native`)
sum(country_score$`4°C Non-native`[order(country_score$`4°C Non-native`,decreasing=T)][1:20])/sum(country_score$`4°C Non-native`)
sum(country_score$`2°C Harmful`[order(country_score$`2°C Harmful`,decreasing=T)][1:10])/sum(country_score$`2°C Harmful`)
sum(country_score$`2°C Harmful`[order(country_score$`2°C Harmful`,decreasing=T)][1:20])/sum(country_score$`2°C Harmful`)
sum(country_score$`4°C Harmful`[order(country_score$`4°C Harmful`,decreasing=T)][1:10])/sum(country_score$`4°C Harmful`)
sum(country_score$`4°C Harmful`[order(country_score$`4°C Harmful`,decreasing=T)][1:20])/sum(country_score$`4°C Harmful`)
sum(country_score$`2°C Environmental`[order(country_score$`2°C Environmental`,decreasing=T)][1:10])/sum(country_score$`2°C Environmental`)
sum(country_score$`2°C Environmental`[order(country_score$`2°C Environmental`,decreasing=T)][1:20])/sum(country_score$`2°C Environmental`)
sum(country_score$`4°C Environmental`[order(country_score$`4°C Environmental`,decreasing=T)][1:10])/sum(country_score$`4°C Environmental`)
sum(country_score$`4°C Environmental`[order(country_score$`4°C Environmental`,decreasing=T)][1:20])/sum(country_score$`4°C Environmental`)
sum(country_score$`2°C Socioeconomic`[order(country_score$`2°C Socioeconomic`,decreasing=T)][1:10])/sum(country_score$`2°C Socioeconomic`)
sum(country_score$`2°C Socioeconomic`[order(country_score$`2°C Socioeconomic`,decreasing=T)][1:20])/sum(country_score$`2°C Socioeconomic`)
sum(country_score$`4°C Socioeconomic`[order(country_score$`4°C Socioeconomic`,decreasing=T)][1:10])/sum(country_score$`4°C Socioeconomic`)
sum(country_score$`4°C Socioeconomic`[order(country_score$`4°C Socioeconomic`,decreasing=T)][1:20])/sum(country_score$`4°C Socioeconomic`)

### for shiny app
NMI_analysis <- cbind(NMI_analysis,Score_analysis)
write.csv(NMI_analysis, "ant_indoor_analysis.csv")