### Impact projections under climate change
Score_analysis <- Impact_score[match(NMI_analysis$species,gsub(" ",".",Impact_score$Species.name)),-1]
Score_analysis[is.na(Score_analysis)] <- 0

NMI_analysis$Impact_4C <- NMI_analysis$proj_diff_4C*Score_analysis #although here used proj_diff_4C, later analyses will only include indoor populations. Therefore the change will be 100% contributed by indoor population, not outdoor population experiencing altered climatic suitability.
NMI_analysis$Impact_2C <- NMI_analysis$proj_diff_2C*Score_analysis
NMI_analysis$Impact_current <- NMI_analysis$num*Score_analysis # num = 1 for naturalized species
NMI_analysis <- do.call(cbind,NMI_analysis)

site_impact_summary <- NMI_analysis %>% 
  group_by(polygon_name,ID) %>% 
  summarize(impact_total_current_E = sum(Impact_current.E.Total,na.rm=T)[num=1],
            impact_total_2C_E = sum(Impact_2C.E.Total[num == 0]),
            impact_total_4C_E = sum(Impact_4C.E.Total[num == 0]),
            impact_increase_2C_E = impact_total_2C_E/impact_total_current_E*100,
            impact_increase_4C_E = impact_total_4C_E/impact_total_current_E*100,
            impact_total_current_S = sum(Impact_current.S.Total)[num=1],
            impact_total_2C_S = sum(Impact_2C.S.Total[num == 0]),
            impact_total_4C_S = sum(Impact_4C.S.Total[num == 0]),
            impact_increase_2C_S = impact_total_2C_S/impact_total_current_S*100,
            impact_increase_4C_S = impact_total_4C_S/impact_total_current_S*100,
            current_Harmful_indoor = sum(dummy[Harmful=="Harmful" & num == 0])) %>% filter(current_Harmful_indoor > 0)

mean(site_impact_summary$impact_total_2C_E)/mean(site_impact_summary$impact_total_current_E)*100 #percentage gain in E_score in 2C scenario
mean(site_impact_summary$impact_total_4C_E)/mean(site_impact_summary$impact_total_current_E)*100 #percentage gain in E_score in 4C scenario
mean(site_impact_summary$impact_total_2C_E)
mean(site_impact_summary$impact_total_4C_E)
min(site_impact_summary$impact_total_2C_E)
min(site_impact_summary$impact_total_4C_E)
max(site_impact_summary$impact_total_2C_E)
max(site_impact_summary$impact_total_4C_E)

mean(site_impact_summary$impact_total_2C_S)/mean(site_impact_summary$impact_total_current_S)*100 #percentage gain in S_score in 2C scenario
mean(site_impact_summary$impact_total_4C_S)/mean(site_impact_summary$impact_total_current_S)*100 #percentage gain in S_score in 4C scenario
mean(site_impact_summary$impact_total_2C_S)
mean(site_impact_summary$impact_total_4C_S)
min(site_impact_summary$impact_total_2C_S)
min(site_impact_summary$impact_total_4C_S)
max(site_impact_summary$impact_total_2C_S)
max(site_impact_summary$impact_total_4C_S)

### Figure S2
Global_impact <- as_tibble(data.frame(polygon_name=NMI_analysis$polygon_name,num=NMI_analysis$num,Harmful=NMI_analysis$Harmful,warm.4C=NMI_analysis$proj_diff_4C*Score_analysis,warm.2C=NMI_analysis$proj_diff_2C*Score_analysis)) %>%
                 group_by(polygon_name) %>%
                 filter(num==0 & Harmful == "Harmful") %>%
                 summarize(across(warm.4C.Animal.production:warm.2C.S.Total, function(x) sum(x,na.rm=T))) %>%
                 ungroup() %>%
                 summarize(across(warm.4C.Animal.production:warm.2C.S.Total,list(sum=~sum(.x,na.rm=T),
                                                                                 n=~length(.x[.x>0]))))

Global_impact <- data.frame(Score=t(Global_impact[seq(1,56,by=2)]),
                            n=t(Global_impact[seq(2,56,by=2)]),
                            Warming = c(rep("4°C",14),rep("2°C",14)),
                            Impact = rep(names(NMI_analysis$proj_diff_4C*Score_analysis),2))

Global_impact$Type <- ifelse(Global_impact$Impact == "Animals" | Global_impact$Impact == "Competition" | Global_impact$Impact == "Diseases" | Global_impact$Impact == "Ecosystems" | Global_impact$Impact == "Hybridization" | Global_impact$Impact == "Plants", "Environmental","Socioeconomic")

Global_impact <- subset(Global_impact, Impact != "E.Total" & Impact != "S.Total")

Global_impact$Impact <- recode(Global_impact$Impact,Animal.production = "Animal Production")

Global_impact$Impact <- factor(Global_impact$Impact,levels=c("Plants","Hybridization","Ecosystems","Diseases","Competition","Animals",
                                                             "Social","Infrastructure","Health","Forestry","Crops","Animal Production"))

Global_impact$Warming <- factor(Global_impact$Warming,levels=c("4°C","2°C")) 

pS2 <- ggplot(Global_impact,aes(label=n,y=Impact,x=Score,fill=Warming))+
  geom_bar(stat="identity",position = 'dodge')+
  geom_text(aes(label=n), position=position_dodge(width=0.9), hjust=-0.25,size=2)+
  xlab("Cumulative score increases across regions under warming")+
  ylab("")+
  scale_fill_manual(values=c("#bd0026","#feb24c"))+
  theme_classic()+
  facet_wrap(~Type,scales="free_y")+
  xlim(0,112)+
  theme_bw()+
  theme(legend.position="bottom",
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),
        legend.text=element_text(size=6),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.title = element_text(size=6))

plot(pS2)

ggsave("Figures/FigS2.tiff",width=16,height=8,units="cm",dpi=600,compression="lzw")
### Sub-figures on impacts, including Figure 2, S1, S4
bentity.shp.sf$proj_diff_impact_E_2C_net <- unlist(site_impact_summary[match(bentity.shp.sf$BENTITY2_N,site_impact_summary$polygon_name),"impact_total_2C_E"])
bentity.shp.sf$proj_diff_impact_E_4C_net<- unlist(site_impact_summary[match(bentity.shp.sf$BENTITY2_N,site_impact_summary$polygon_name),"impact_total_4C_E"])
bentity.shp.sf$proj_diff_impact_S_2C_net <- unlist(site_impact_summary[match(bentity.shp.sf$BENTITY2_N,site_impact_summary$polygon_name),"impact_total_2C_S"])
bentity.shp.sf$proj_diff_impact_S_4C_net<- unlist(site_impact_summary[match(bentity.shp.sf$BENTITY2_N,site_impact_summary$polygon_name),"impact_total_4C_S"])
bentity.shp.sf$warming_diff_impact_E <- unlist(site_impact_summary[match(bentity.shp.sf$BENTITY2_N,site_impact_summary$polygon_name),"impact_total_4C_E"]-site_impact_summary[match(bentity.shp.sf$BENTITY2_N,site_impact_summary$polygon_name),"impact_total_2C_E"]) #4C - 2C E
bentity.shp.sf$warming_diff_impact_S <- unlist(site_impact_summary[match(bentity.shp.sf$BENTITY2_N,site_impact_summary$polygon_name),"impact_total_4C_S"]-site_impact_summary[match(bentity.shp.sf$BENTITY2_N,site_impact_summary$polygon_name),"impact_total_2C_S"]) #4C - 2C S

bentity.shp.sf$current_impact_E<- unlist(site_impact_summary[match(bentity.shp.sf$BENTITY2_N,site_impact_summary$polygon_name),"impact_total_current_E"])
bentity.shp.sf$current_impact_S<- unlist(site_impact_summary[match(bentity.shp.sf$BENTITY2_N,site_impact_summary$polygon_name),"impact_total_current_S"])

library(tiff)
library(grid)
warming_tiff <- readTIFF("Figures/warming.tif")
g3 <- rasterGrob(warming_tiff,width = unit(0.7,"cm"),height=unit(0.7,"cm"),interpolate=TRUE)

xmax <- xmin <- -1.4e07
ymin <- ymax <- -5250000
space <- 900000
size <- 3

p1e <- ggplot(data=bentity.shp.sf,aes(fill=current_impact_E))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(e) Environmental",size=size)+
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,max(bentity.shp.sf$current_impact_E,na.rm=T)))+
  theme
plot(p1e)

p1f <- ggplot(data=bentity.shp.sf,aes(fill=current_impact_S))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(f) Socioeconomic",size=size-0.2)+
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,max(bentity.shp.sf$current_impact_S,na.rm=T)))+
  theme
plot(p1f)

library(ggpubr)
p1 <- ggarrange(p1a,p1b,p1c,p1d,p1e,p1f,hjust=0,vjust=0,label.x=0,label.y=0,nrow=3,ncol=2)
plot(p1)
ggsave("Figures/FigS1.tiff",dpi=800,compression="lzw",units="cm",height=9.5*1.5,width=16.8,bg="white")

###

p2e <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_impact_E_2C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(e) Environmental",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,12.2))+
  theme

p2f <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_impact_E_4C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(f) Environmental",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,12.2))+
  theme

p2g <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_impact_S_2C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(g) Socioeconomic",size=size-0.2)+
  annotate("text", x = xmin+space, y = ymin+space,label = "2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,7.2))+
  theme

p2h  <- ggplot(data=bentity.shp.sf,aes(fill=proj_diff_impact_S_4C_net))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(h) Socioeconomic",size=size-0.2)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.1,7.2))+
  theme

pS4c <- ggplot(data=bentity.shp.sf,aes(fill=warming_diff_impact_E))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(c) Environmental",size=size)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C vs 2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.01,8.4))+
  theme

pS4d <- ggplot(data=bentity.shp.sf,aes(fill=warming_diff_impact_S))+
  geom_sf(colour="black",linewidth=0.1)+
  annotate("text", x = -Inf, y = Inf, hjust=0,vjust=1,label = "(d) Socioeconomic",size=size-0.2)+
  annotate("text", x = xmin+space, y = ymin+space,label = "4°C vs 2°C",colour="red",size=size,hjust=0)+
  annotation_custom(g3, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
  labs(fill="")+
  xlim(-14800000,14800000)+
  ylim(-6500000,9000000)+
  scale_fill_continuous(low="#ffffb2",high="#bd0026",na.value="white",limits=c(-0.01,4.5))+
  theme

p2 <- ggarrange(p2a,p2b,p2c,p2d,p2e,p2f,p2g,p2h,nrow=4,ncol=2)
plot(p2)
ggsave("Figures/Fig2.tiff",dpi=800,compression="lzw",units="cm",height=13.5/3*4.5,width=16.8,bg="white")

###Species level summary on impactsm and Fig S2
Species_summary <- NMI_analysis %>% group_by(species) %>% filter(Harmful == "Harmful" & num == 0)  %>% summarize(increase_2C = sum(proj_diff_2C),
                                                                                                                   increase_4C = sum(proj_diff_4C),
                                                                                                                   Impact_E_2C = sum(Impact_2C.E.Total),
                                                                                                                   Impact_E_4C = sum(Impact_4C.E.Total),
                                                                                                                   Impact_S_2C = sum(Impact_2C.S.Total),
                                                                                                                   Impact_S_4C = sum(Impact_4C.S.Total),
                                                                                                                   indoor_region = n())
Species_summary_long <- Species_summary %>% pivot_longer(!species)

Species_summary_long$Warming <- ifelse(grepl("2C",Species_summary_long$name),"2°C","4°C")
Species_summary_long$Type <- ifelse(grepl("_E_",Species_summary_long$name),"Environmental",
ifelse(grepl("_S",Species_summary_long$name),"Socioeconomic","Naturalization probablity"))

Species_summary_long$Warming <- factor(Species_summary_long$Warming,levels=c("4°C","2°C")) 
Species_summary_long$Type <- factor(Species_summary_long$Type,levels=c("Naturalization probablity","Environmental","Socioeconomic")) 
Species_summary_long$species <- factor(Species_summary_long$species,levels= sort(unique(Species_summary_long$species),decreasing=T))

pS3 <- ggplot(data=Species_summary_long,aes(y=species,x=value,fill=Warming))+
  geom_bar(stat="Identity",position="dodge")+
  facet_wrap(~Type,scale="free_x")+
  xlab("Cumulative increases across indoor populations under warming")+
  ylab("Species (Number of indoor populations)")+
  scale_fill_manual(values=c("#bd0026","#feb24c"))+
  scale_y_discrete(labels = paste0("<i>",gsub("\\."," ",levels(Species_summary_long$species)),"</i>"," (",rev(Species_summary$indoor_region),")"))+
  theme_bw()+
  theme(legend.position="bottom",
        axis.text.x=element_text(size=6),
        axis.text.y=ggtext::element_markdown(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),
        legend.text=element_text(size=6),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.title = element_text(size=6))
plot(pS3)
ggsave("Figures/FigS3.tiff",dpi=800,compression="lzw",units="cm",height=13.5,width=16.8,bg="white")

###

pS4<- ggarrange(pS4a,pS4b,pS4c,pS4d)
plot(pS4)
ggsave("Figures/FigS4.tiff",dpi=800,compression="lzw",units="cm",height=8.4*1.25,width=16.8,bg="white")

#pairwise correlations between gain in alien, harmful, E and S impacts between the two scenarios. Only consider polygon with >= 1 populations.
cor(subset(site_summary,current_indoor != 0)$proj_diff_indoor_2C_net,subset(site_summary,current_indoor != 0)$proj_diff_indoor_4C_net,method="kendall") 
cor(subset(site_summary,current_indoor_Harmful != 0)$proj_diff_Harmful_indoor_2C_net,subset(site_summary,current_indoor_Harmful != 0)$proj_diff_Harmful_indoor_4C_net,method="kendall")
cor(site_impact_summary$impact_total_2C_E,site_impact_summary$impact_total_4C_E,method="kendall")
cor(site_impact_summary$impact_total_2C_S,site_impact_summary$impact_total_4C_S,method="kendall")
