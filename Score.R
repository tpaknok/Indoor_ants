NMI_analysis$Final_country <- bentity.shp.df$Final_country[match(NMI_analysis$polygon_name,bentity.shp.df$BENTITY2_N)]
NMI_analysis <- cbind(NMI_analysis,capacity[match(NMI_analysis$Final_country,capacity$Country_global),c("threat","IAS_list","Existing_mgmt","research","outreach")])

NMI_analysis$total_score <- apply(NMI_analysis[,c("threat","IAS_list","Existing_mgmt","research","outreach")],1,sum)

length(which(is.na(NMI_analysis$total_score)))

sp_country <- NMI_analysis %>% group_by(species,Final_country,invasive,total_score) %>% filter(num == 0) %>% summarise(proj_diff_2C = max(proj_diff_2C),
                                                                                                                       proj_diff_4C = max(proj_diff_4C),
                                                                                                                       Impact_2C.E.Total = max(Impact_2C.E.Total),
                                                                                                                       Impact_4C.E.Total = max(Impact_4C.E.Total),
                                                                                                                       Impact_2C.S.Total = max(Impact_2C.S.Total),
                                                                                                                       Impact_4C.S.Total = max(Impact_4C.S.Total))

record_score <- sp_country %>% group_by(total_score) %>% summarise("2°C Alien" = sum(proj_diff_2C),
                                                                   "4°C Alien" = sum(proj_diff_4C),
                                                                   "2°C Harmful" = sum(proj_diff_2C[invasive == "Invasive"]),
                                                                   "4°C Harmful" = sum(proj_diff_4C[invasive=="Invasive"]),
                                                                   "2°C Environmental" = sum(Impact_2C.E.Total),
                                                                   "4°C Environmental" = sum(Impact_4C.E.Total),
                                                                   "2°C Socioeconomic" = sum(Impact_2C.S.Total),
                                                                   "4°C Socioeconomic" = sum(Impact_4C.S.Total))

record_score <- record_score %>% pivot_longer(!total_score)
na.omit(record_score) %>% group_by(name) %>% summarize(percent = value[total_score==5]/sum(value)*100)

record_score$name <- factor(record_score$name,levels=c("2°C Alien","4°C Alien", "2°C Harmful","4°C Harmful","2°C Environmental","4°C Environmental","2°C Socioeconomic","4°C Socioeconomic"))
p6 <- ggplot(data=record_score)+
  geom_bar(aes(x=total_score,weight=value))+
  #geom_line(aes(x=total_score,y=value))+
  xlab("Total response capcity score (0-5)")+
  ylab("Increases driven by climate change")+
  facet_wrap(~name,scales="free",nrow=4,ncol=2)+
  theme_classic()+
  scale_x_continuous(breaks=c(0,1,2,3,4,5))+
  theme(legend.position="bw")

plot(p6)
ggsave("Indoor_ants/Fig6.tiff",dpi=800,width=16.8,height=16.8,units="cm",compression="lzw",bg="white")

country_record_score <- sp_country %>% group_by(total_score,Final_country) %>% summarise("2°C Alien" = sum(proj_diff_2C),
                                                                           "4°C Alien" = sum(proj_diff_4C),
                                                                           "2°C Harmful" = sum(proj_diff_2C[invasive == "Invasive"]),
                                                                           "4°C Harmful" = sum(proj_diff_4C[invasive=="Invasive"]),
                                                                           "2°C Environmental" = sum(Impact_2C.E.Total),
                                                                           "4°C Environmental" = sum(Impact_4C.E.Total),
                                                                           "2°C Socioeconomic" = sum(Impact_2C.S.Total),
                                                                           "4°C Socioeconomic" = sum(Impact_4C.S.Total))

sum_stat <- country_record_score %>% group_by(total_score) %>% summarise(n= n(),
                                                                   mean_2C_alien = mean(`2°C Alien`),
                                                                   sd_2C_alien = sd(`2°C Alien`),
                                                                   mean_4C_alien = mean(`4°C Alien`),
                                                                   sd_4C_alien = sd(`4°C Alien`),
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
country_record_score_long$scenario <- factor(country_record_score_long$scenario,levels=c("2°C Alien","4°C Alien", "2°C Harmful","4°C Harmful","2°C Environmental","4°C Environmental","2°C Socioeconomic","4°C Socioeconomic"))

pS1 <- ggplot(data=country_record_score_long)+
  geom_boxplot(aes(x=total_score,y=value,group=total_score))+
  #geom_jitter(aes(x=total_score,y=value),width=0.1)+
  xlab("Total response capacity score (0-5)")+
  ylab("Increases driven by climate change")+
  facet_wrap(~scenario,scales="free",nrow=4,ncol=2)+
  theme_classic()+
  theme(legend.position="bw")

plot(pS1)
ggsave("Indoor_ants/FigS1.tiff",dpi=800,width=16.8,height=16.8,units="cm",compression="lzw",bg="white")

##########regression not needed
#library(glmmTMB)
#library(performance)
#r2(glmmTMB(`2°C Alien`~scale(total_score)+I(scale(total_score)^2),disp=~scale(total_score)+I(scale(total_score)^2),data=country_record_score))
#summary(glmmTMB(`2°C Alien`~scale(total_score)+I(scale(total_score)^2),disp=~scale(total_score)+I(scale(total_score)^2),data=country_record_score))
#summary(glmmTMB(`2°C Invasive`~scale(total_score)+I(scale(total_score)^2),disp=~scale(total_score)+I(scale(total_score)^2),data=country_record_score))
#summary(glmmTMB(`4°C Alien`~scale(total_score)+I(scale(total_score)^2),disp=~scale(total_score)+I(scale(total_score)^2),data=country_record_score))
#summary(glmmTMB(`4°C Invasive`~scale(total_score)+I(scale(total_score)^2),disp=~scale(total_score)+I(scale(total_score)^2),data=country_record_score))



#####################################################################################################################
country_score <- sp_country %>% group_by(Final_country) %>% summarise("2°C Alien" = sum(proj_diff_2C),
                                                                      "4°C Alien" = sum(proj_diff_4C),
                                                                      "2°C Harmful" = sum(proj_diff_2C[invasive == "Invasive"]),
                                                                      "4°C Harmful" = sum(proj_diff_4C[invasive=="Invasive"]),
                                                                      "2°C Environmental" = sum(Impact_2C.E.Total),
                                                                      "4°C Environmental" = sum(Impact_4C.E.Total),
                                                                      "2°C Socioeconomic" = sum(Impact_2C.S.Total),
                                                                      "4°C Socioeconomic" = sum(Impact_4C.S.Total))

country_score <- cbind(country_score,capacity[match(country_score$Final_country,capacity$Country_global),c("threat","IAS_list","Existing_mgmt","research","outreach")])
var = c("2°C Alien","4°C Alien","2°C Harmful","4°C Harmful","2°C Environmental","4°C Environmental","2°C Socioeconomic","4°C Socioeconomic")
top_n <- 20
country_score <- na.omit(country_score)
country_score$Final_country[country_score$Final_country == "Dem. Rep. Korea"] <- "North Korea"
country_score$Final_country[country_score$Final_country == "Republic of Korea"] <- "South Korea"

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
                                       scale_fill_manual(values=c("#ffffb2","#fd8d3c","#bd0026"))+
                                       ggtitle(var[[num]])+
                                       theme_classic()+
                                       theme(axis.text.x = if(num==1){element_text(angle = 10,size=7)}else{element_text(angle=10,colour = "white",size=8)},
                                             axis.text.y = element_text(size=5.5)))
  }

p5 <- ggarrange(p_cap_1,p_cap_2,p_cap_3,p_cap_4,p_cap_5,p_cap_6,p_cap_7,p_cap_8,nrow=4,ncol=2,common.legend=T,legend="bottom")
plot(p5)

ggsave("Indoor_ants/Fig5.tiff",dpi=800,width=16.8,height=24,units="cm",compression="lzw",bg="white")
write.csv(prop_df,"Indoor_ants/prop.csv")
##########################
sum(country_score$`2°C Alien`[order(country_score$`2°C Alien`,decreasing=T)][1:10])/sum(country_score$`2°C Alien`)
sum(country_score$`2°C Alien`[order(country_score$`2°C Alien`,decreasing=T)][1:20])/sum(country_score$`2°C Alien`)
sum(country_score$`4°C Alien`[order(country_score$`4°C Alien`,decreasing=T)][1:10])/sum(country_score$`4°C Alien`)
sum(country_score$`4°C Alien`[order(country_score$`4°C Alien`,decreasing=T)][1:20])/sum(country_score$`4°C Alien`)
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
