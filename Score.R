NMI_analysis$Final_country <- bentity.shp.df$Final_country[match(NMI_analysis$polygon_name,bentity.shp.df$BENTITY2_N)]
NMI_analysis <- cbind(NMI_analysis,capacity[match(NMI_analysis$Final_country,capacity$Country_global),c("threat","IAS_list","Existing_mgmt","research","outreach")])

NMI_analysis$total_score <- apply(NMI_analysis[,c("threat","IAS_list","Existing_mgmt","research","outreach")],1,sum)

length(which(is.na(NMI_analysis$total_score)))

sp_country <- NMI_analysis %>% group_by(species,Final_country,invasive,total_score) %>% filter(num == 0) %>% summarise(proj_diff_2C = max(proj_diff_2C),
                                                                                                                       proj_diff_4C = max(proj_diff_4C))

record_score <- sp_country %>% group_by(total_score) %>% summarise("2°C Alien" = sum(proj_diff_2C),
                                                                   "4°C Alien" = sum(proj_diff_4C),
                                                                   "2°C Invasive" = sum(proj_diff_2C[invasive == "Invasive"]),
                                                                   "4°C Invasive" = sum(proj_diff_4C[invasive=="Invasive"]))

record_score <- record_score %>% pivot_longer(!total_score)
record_score %>% group_by(name) %>% summarize(percent = value[total_score==5]/sum(value)*100)

p6 <- ggplot(data=record_score)+
  geom_point(aes(x=total_score,y=value))+
  geom_line(aes(x=total_score,y=value))+
  xlab("Total score (0-5)")+
  ylab("Effects of climate change")+
  facet_wrap(~name,scales="free")+
  theme_classic()+
  theme(legend.position="bw")

plot(p6)
ggsave("Fig6.tiff",dpi=800,width=16.8,height=16.8,units="cm",compression="lzw",bg="white")

country_record_score <- sp_country %>% group_by(total_score,Final_country) %>% summarise("2°C Alien" = sum(proj_diff_2C),
                                                                           "4°C Alien" = sum(proj_diff_4C),
                                                                           "2°C Invasive" = sum(proj_diff_2C[invasive == "Invasive"]),
                                                                           "4°C Invasive" = sum(proj_diff_4C[invasive=="Invasive"]))

sum_stat <- country_record_score %>% group_by(total_score) %>% summarise(n= n(),
                                                                   mean_2C_alien = mean(`2°C Alien`),
                                                                   sd_2C_alien = sd(`2°C Alien`),
                                                                   mean_4C_alien = mean(`4°C Alien`),
                                                                   sd_4C_alien = sd(`4°C Alien`),
                                                                   mean_2C_invasive = mean(`2°C Invasive`),
                                                                   sd_2C_invasive = sd(`2°C Invasive`),
                                                                   mean_4C_invasive = mean(`4°C Invasive`),
                                                                   sd_4C_invasive = sd(`4°C Invasive`),)


country_record_score_long <- na.omit(country_record_score) %>% pivot_longer(cols=-c(Final_country,total_score),names_to="scenario")

pS1 <- ggplot(data=country_record_score_long)+
  geom_boxplot(aes(x=total_score,y=value,group=total_score))+
  #geom_jitter(aes(x=total_score,y=value),width=0.1)+
  xlab("Total score (0-5)")+
  ylab("Effects of climate change")+
  facet_wrap(~scenario,scales="free")+
  theme_classic()+
  theme(legend.position="bw")

plot(pS1)
ggsave("FigS1.tiff",dpi=800,width=16.8,height=16.8,units="cm",compression="lzw",bg="white")

library(glmmTMB)
library(performance)
r2(glmmTMB(`2°C Alien`~scale(total_score)+I(scale(total_score)^2),disp=~scale(total_score)+I(scale(total_score)^2),data=country_record_score))
summary(glmmTMB(`2°C Alien`~scale(total_score)+I(scale(total_score)^2),disp=~scale(total_score)+I(scale(total_score)^2),data=country_record_score))
summary(glmmTMB(`2°C Invasive`~scale(total_score)+I(scale(total_score)^2),disp=~scale(total_score)+I(scale(total_score)^2),data=country_record_score))
summary(glmmTMB(`4°C Alien`~scale(total_score)+I(scale(total_score)^2),disp=~scale(total_score)+I(scale(total_score)^2),data=country_record_score))
summary(glmmTMB(`4°C Invasive`~scale(total_score)+I(scale(total_score)^2),disp=~scale(total_score)+I(scale(total_score)^2),data=country_record_score))



#####################################################################################################################
country_score <- sp_country %>% group_by(Final_country) %>% summarise("2°C Alien" = sum(proj_diff_2C),
                                                                      "4°C Alien" = sum(proj_diff_4C),
                                                                      "2°C Invasive" = sum(proj_diff_2C[invasive == "Invasive"]),
                                                                      "4°C Invasive" = sum(proj_diff_4C[invasive=="Invasive"]))
country_score <- na.omit(country_score)

country_score <- cbind(country_score,capacity[match(country_score$Final_country,capacity$Country_global),c("threat","IAS_list","Existing_mgmt","research","outreach")])
var = c("2°C Alien","4°C Alien","2°C Invasive","4°C Invasive" )
top_n <- 20
country_score <- na.omit(country_score)

for (num in 1:4){
  top_20_pos <- order(country_score[,var[[num]]],decreasing=T)[1:top_n]
  plot_df <- country_score[top_20_pos,]
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
                                       scale_fill_manual(values=c("#ffffb2","#fd8d3c","#bd0026"))+
                                       ggtitle(var[[num]])+
                                       theme_classic()+
                                       theme(axis.text.x = if(num==1){element_text(angle = 10,size=6)}else{element_text(angle=10,colour = "white",size=6)}))
  }

p5 <- ggarrange(p_cap_1,p_cap_3,p_cap_2,p_cap_4,common.legend=T,legend="bottom")
ggsave("Fig5.tiff",dpi=800,width=16.8,height=16.8,units="cm",compression="lzw",bg="white")
