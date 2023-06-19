library(remotes)
install_github('vincentarelbundock/WDI')
library(WDI)

search_result <- WDIsearch('gdp.*constant')

GDP_dat = WDI(indicator='NY.GDP.MKTP.KD', start=2021,end=2022)
GDP_dat$country[GDP_dat$country == "Bahamas, The"] <- "Bahamas"
GDP_dat$country[GDP_dat$country == "Bosnia and Herzegovina"] <- "Bosnia and Herz."
GDP_dat$country[GDP_dat$country == "Cayman Islands"] <- "Cayman Is."
GDP_dat$country[GDP_dat$country == "Central African Republic"] <- "Central African Rep."
GDP_dat$country[GDP_dat$country == "Korea, Dem. People's Rep."] <- "Dem. Rep. Korea"
GDP_dat$country[GDP_dat$country == "Congo, Dem. Rep."] <- "Democratic Republic of the Congo"
GDP_dat$country[GDP_dat$country == "Egypt, Arab Rep."] <- "Egypt"
GDP_dat$country[GDP_dat$country == "Equatorial Guinea"] <- "Eq. Guinea"
GDP_dat$country[GDP_dat$country == "French Polynesia"] <- "Fr. Polynesia"
GDP_dat$country[GDP_dat$country == "Gambia, The"] <- "Gambia"
GDP_dat$country[GDP_dat$country == "Hong Kong SAR, China"] <- "Hong Kong"
GDP_dat$country[GDP_dat$country == "Iran, Islamic Rep."] <- "Iran"
GDP_dat$country[GDP_dat$country == "Channel Islands"] <- "Jersey"
GDP_dat$country[GDP_dat$country == "Macao SAR, China"] <- "Macao"
GDP_dat$country[GDP_dat$country == "Marshall Islands"] <- "Marshall Is."
GDP_dat$country[GDP_dat$country == "Micronesia, Fed. Sts."] <- "Micronesia"
GDP_dat$country[GDP_dat$country == "Korea, Rep."] <- "Republic of Korea"
GDP_dat$country[GDP_dat$country == "Congo, Rep."] <- "Republic of the Congo"
GDP_dat$country[GDP_dat$country == "Russia"] <- "Russian Federation"
GDP_dat$country[GDP_dat$country == "Slovak Republic"] <- "Slovakia"
GDP_dat$country[GDP_dat$country == "Solomon Islands"] <- "Solomon Is."
GDP_dat$country[GDP_dat$country == "Syrian Arab Republic"] <- "Syria"
GDP_dat$country[GDP_dat$country == "Solomon Islands"] <- "Solomon Is."
GDP_dat$country[GDP_dat$country == "Kyrgyz Republic"] <- "Taiwan"
GDP_dat$country[GDP_dat$country == "Lao PDR"] <- "Laos"
GDP_dat$country[GDP_dat$country == "Turkiye"] <- "Turkey"
GDP_dat$country[GDP_dat$country == "Turks and Caicos Islands"] <- "Turks and Caicos Is."
GDP_dat$country[GDP_dat$country == "Venezuela, RB"] <- "Venezuela"
GDP_dat$country[GDP_dat$country == "Yemen, Rep."] <- "Yemen"

site_summary$indoor.sr <- site_summary$n - site_summary$current.sum

gain_df <- data.frame(shp = site_summary$polygon_name,
                      shp.ID = site_summary$shp.ID,
                      total.sr = site_summary$n,
                      indoor.sr=site_summary$indoor.sr,
                      indoor.sr.invasive=site_summary_invasive$n-site_summary_invasive$current.sum,
                      alien.gain = site_summary$proj_diff,
                      alien.gain.4C = site_summary_4C$proj_diff)

gain_df <- data.frame(gain_df,
                 invasive.gain = site_summary_invasive[match(gain_df$shp,site_summary_invasive$polygon_name),"proj_diff"],
                 invasive.gain.4C = site_summary_invasive_4C[match(gain_df$shp,site_summary_invasive_4C$polygon_name),"proj_diff"])

colnames(gain_df) <- c("shp","shp.ID","sr","indoor.sr","indoor.sr.invasive","alien.gain","alien.gain.4C","invasive.gain","invasive.gain.4C")

gain_df$Country <- bentity.shp.df[match(gain_df$shp,bentity.shp.df$BENTITY2_N),"Final_country"]

gain_df$cap.R <- capacity[match(gain_df$Country,capacity$Country_global_modified),"reactive"]
gain_df$cap.P <- capacity[match(gain_df$Country,capacity$Country_global_modified),"proactive"]

gain_df$cap <- gain_df$cap.R+gain_df$cap.P

gain_df$temp_pca_change <- future_df[match(gain_df$shp,future_df$polygon_name),"scaled.temp_pca"]-future_df[match(gain_df$shp,future_df$polygon_name),"scaled.temp_pca_current"]

cor(gain_df$alien.gain,gain_df$invasive.gain,method="pearson")
cor(gain_df$alien.gain.4C,gain_df$invasive.gain.4C,method="pearson")

##########
gain_df$normal.cap <- (gain_df$cap - min(gain_df$cap,na.rm=T))/(max(gain_df$cap,na.rm=T)-min(gain_df$cap,na.rm=T))
gain_df$normal.invasive.4C.gain <- (gain_df$invasive.gain.4C - min(gain_df$invasive.gain.4C,na.rm=T))/(max(gain_df$invasive.gain.4C,na.rm=T)-min(gain_df$invasive.gain.4C,na.rm=T))
gain_df$normal.invasive.gain <- (gain_df$invasive.gain - min(gain_df$invasive.gain,na.rm=T))/(max(gain_df$invasive.gain,na.rm=T)-min(gain_df$invasive.gain,na.rm=T))

gain_df$shortfall_4C <- gain_df$normal.invasive.4C.gain-gain_df$normal.cap

bentity.df$shortfall_4C <- gain_df[match(bentity.df$geom,gain_df$shp.ID),"shortfall_4C"]
bentity.df$normal.cap <- gain_df[match(bentity.df$geom,gain_df$shp.ID),"normal.cap"]

#p4b <- ggplot(data=bentity.df,aes(x=x,y=y,group=id,fill=shortfall_4C))+
  #geom_polygon(colour="black")+
  #scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-1,0.35))+
  #labs(fill="Shortfall")+
  #theme+
  #coord_equal()
#plot(p4b)

###########

avg_country <- aggregate(cbind(invasive.gain,invasive.gain.4C,alien.gain,alien.gain.4C)~Country,data=gain_df,FUN=max,na.rm=T)
Country_summary_2C_invasive <- cbind(Country_summary_2C_invasive,avg_country[match(Country_summary_2C_invasive$Final_country,avg_country$Country),c("invasive.gain","alien.gain")])
Country_summary_2C_invasive$cap.R <- capacity[match(Country_summary_2C_invasive$Final_country,capacity$Country_global),"reactive"]
Country_summary_2C_invasive$cap.P <- capacity[match(Country_summary_2C_invasive$Final_country,capacity$Country_global),"proactive"]
Country_summary_2C_invasive$cap <- Country_summary_2C_invasive$cap.P+Country_summary_2C_invasive$cap.R
Country_summary_2C_invasive <- cbind(Country_summary_2C_invasive,capacity[match(Country_summary_2C_invasive$Final_country,capacity$Country_global),7:12])
Country_summary_2C_invasive$GDP <- GDP_dat[match(Country_summary_2C_invasive$Final_country,GDP_dat$country),"NY.GDP.MKTP.KD"]
shortfall_country_2C <- Country_summary_2C_invasive[!is.na(Country_summary_2C_invasive$cap),]
shortfall_country_2C$warming <- "2C"


########################################

Country_summary_4C_invasive <- cbind(Country_summary_4C_invasive,avg_country[match(Country_summary_4C_invasive$Final_country,avg_country$Country),c("invasive.gain.4C","alien.gain.4C")])
Country_summary_4C_invasive$cap.R <- capacity[match(Country_summary_4C_invasive$Final_country,capacity$Country_global),"reactive"]
Country_summary_4C_invasive$cap.P <- capacity[match(Country_summary_4C_invasive$Final_country,capacity$Country_global),"proactive"]
Country_summary_4C_invasive$cap <- Country_summary_4C_invasive$cap.P+Country_summary_4C_invasive$cap.R
Country_summary_4C_invasive <- cbind(Country_summary_4C_invasive,capacity[match(Country_summary_4C_invasive$Final_country,capacity$Country_global),7:12])
Country_summary_4C_invasive$GDP <- GDP_dat[match(Country_summary_4C_invasive$Final_country,GDP_dat$country),"NY.GDP.MKTP.KD"]
shortfall_country_4C <- Country_summary_4C_invasive[!is.na(Country_summary_4C_invasive$cap),]
shortfall_country_4C$warming <- "4C"

var <- c("outreach","research","Existing_mgmt","IAS_list","threat")
reg_df <- list()
reg_csv <- list()
list_num <- 1

for (response in 1:2) {
  for (num in 1:5) {
    reg_data <- if (response == 1) {get("Country_summary_2C_invasive")} else {get("Country_summary_4C_invasive")}
    reg_data <- subset(reg_data,indoor_sr > 0)
    formula <- as.formula(paste0(ifelse(response==1,"invasive.gain~","invasive.gain.4C~"),var[[num]]))
    m_reg <- glmmTMB(formula,disp=as.formula(paste0("~",var[[num]])),data=reg_data)
    p <- summary(m_reg)$coefficients$cond[2,4]
    slope <- summary(m_reg)$coefficients$cond[2,1]
    reg_df[[list_num]] <- data.frame(response=ifelse(response==1,"gain.2C","gain.4C"),predictor=var[[num]],p=p,warming=ifelse(response == 1,"2°C","4°C"))
    reg_csv[[list_num]] <- data.frame(summary(m_reg)$coefficients$cond[,c(1:2,4)],r2(m_reg)$R2_marginal)
    list_num <- list_num+1
  }
}

reg_df <- do.call(rbind,reg_df)
reg_df$id <- paste0("Gain_IAS","_",reg_df$predictor,"_",reg_df$warming)
write.csv(round(do.call(rbind,reg_csv),3),"reg_csv_alien.csv")

long_df <- subset(Country_summary_2C_invasive,indoor_sr > 0) %>% select(IAS_list,threat,Existing_mgmt,research,outreach,"invasive.gain") %>% na.omit() %>% pivot_longer(!invasive.gain) %>% mutate(warming = "2°C")
long_df_4C <- subset(Country_summary_4C_invasive,indoor_sr > 0) %>% select(IAS_list,threat,Existing_mgmt,research,outreach,"invasive.gain.4C") %>% na.omit() %>% pivot_longer(!invasive.gain.4C) %>% mutate(warming = "4°C")
colnames(long_df_4C)[[1]] <- colnames(long_df)[[1]] <- "Gain_IAS"

long_df <- rbind(long_df,long_df_4C)
long_df$id <- paste0("Gain_IAS","_",long_df$name,"_",long_df$warming)
long_df$p <- reg_df[match(long_df$id,reg_df$id),"p"]
long_df$sig <- ifelse(long_df$p < 0.05, "sig","insig")
long_df$name <- as.factor(long_df$name)
long_df$name <- recode_factor(long_df$name,Existing_mgmt = "Control",IAS_list = "List",outreach = "Monitoring",research = "Research",threat = "Threat")

p5 <- ggplot(data=long_df,aes(x=value,y=Gain_IAS))+
  facet_wrap(name~warming,nrow=5,ncol=2,scales="free_y")+
  geom_jitter(aes(group=value),width=0.1)+
  geom_smooth(aes(linetype=sig),method="lm",se=F)+ #basically the same estimate as glmmTMB. use lm here for simplicity
  ylab("Projected gain in invasive species in wild environments")+
  xlab("Management capacity score")+
  scale_x_continuous(breaks=c(0,0.5,1))+
  scale_linetype_manual(values=c("dashed","solid"),guide = "none")+
  scale_color_manual(values=c("blue","red"),name="Tau",labels=c("Median","90th percentile"))+
  theme_bw()+
  theme(legend.position="bottom")

plot(p5)

ggsave("Fig4.tiff",width=16.8,height=20,dpi=800,units="cm",compression="lzw",bg="white")
#############################
inv_df <- list(Country_summary_2C_invasive,Country_summary_4C_invasive)
topn = c(10,20,30)
for (scen in 1:2) {
  temp_df <- inv_df[[scen]]
  for (i in 1:3) {
    long_df <- temp_df %>% select(Final_country,IAS_list,threat,Existing_mgmt,research,outreach,ifelse(scen==1,"max.gain.2C","max.gain.4C")) %>% na.omit() %>% top_n(n=topn[[i]]) %>% pivot_longer(!Final_country)
    temp_gain <- long_df[long_df$name == ifelse(scen==1,"max.gain.2C","max.gain.4C"),]
    temp_gain <- temp_gain$Final_country[order(temp_gain$value,decreasing=F)]
    long_df$Final_country <- factor(long_df$Final_country,levels=temp_gain)
    
    long_df <- long_df[long_df$name != ifelse(scen==1,"max.gain.2C","max.gain.4C"),]
    long_df$value <- as.factor(long_df$value)
    title <- ifelse(i != 2,"",ifelse(scen==1,"2°C","4°C"))
    assign(paste0("p",topn[[i]],"_",ifelse(scen==1,"2C","4C")),ggplot(data=long_df,aes(x=name,y=Final_country,fill=value))+
                  geom_tile(color="black")+
                  scale_fill_manual(values=c("blue","purple","red"),name="Score")+
                  ggtitle(title)+
                  ylab("")+
                  xlab("")+
                  scale_x_discrete(labels = c('Control','List','Monitoring','Research','Threat'))+
                  theme_classic()+
                  theme(axis.text.x = element_text(angle = 15,size=8),
                        axis.text.y = element_text(size=8)))
    
  }
}

ggarrange(p20_2C,p20_4C,ncol=2,nrow=1,common.legend=T,legend="bottom")
ggsave("Fig5.tiff",width=16.8,height=12.6,units="cm",compression="lzw",bg="white")

###################################################
country_geom <- as.data.frame(geom(makeValid(country)))
country_df$brk_name[country_df$brk]
country_geom$country <- country_df$brk_name[country_geom$geom]
country_geom$gain_2C <- unlist(Country_summary_2C_invasive[match(country_geom$country,Country_summary_2C_invasive$Final_country),"proj_diff"])
country_geom$id <- paste0(country_geom$geom,"_",country_geom$part)
country_geom <- subset(country_geom,hole==0)
country_geom[country_geom$geom == 22 & (country_geom$y < 40 | country_geom$y > 55),"gain_2C"] <- NA

p4a <- ggplot(country_geom,aes(x=x,y=y,group=id,fill=gain_2C))+
  geom_polygon(colour="black")+  
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,2))+
  labs(fill="Gain in invasive ant species under 2°C warming")+
  theme+
  coord_equal()
plot(p4a)

country_geom$gain_4C <- unlist(Country_summary_4C_invasive[match(country_geom$country,Country_summary_4C_invasive$Final_country),"proj_diff"])
country_geom[country_geom$geom == 22 & (country_geom$y < 40 | country_geom$y > 55),"gain_4C"] <- NA
p4b <- ggplot(country_geom,aes(x=x,y=y,group=id,fill=gain_4C))+
  geom_polygon(colour="black")+  
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-0.1,2))+
  labs(fill="Gain in invasive ant species under 4°C warming")+
  theme+
  coord_equal()
plot(p4b)

#####################################
temp_shortfall <- subset(shortfall_df,warming=="2C")
country_geom$shortfall_2C <- unlist(temp_shortfall[match(country_geom$country,temp_shortfall$Final_country),"shortfall"])
country_geom[country_geom$geom == 22 & (country_geom$y < 40 | country_geom$y > 55),"shortfall_2C"] <- NA

p4c <- ggplot(country_geom,aes(x=x,y=y,group=id,fill=shortfall_2C))+
  geom_polygon(colour="black")+  
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-1,0.34))+
  labs(fill="Shortfall under 2°C warming")+
  theme+
  coord_equal()
plot(p4c)

temp_shortfall <- subset(shortfall_df,warming=="4C")
country_geom$shortfall_4C <- unlist(temp_shortfall[match(country_geom$country,temp_shortfall$Final_country),"shortfall"])
country_geom[country_geom$geom == 22 & (country_geom$y < 40 | country_geom$y > 55),"shortfall_4C"] <- NA

p4d <- ggplot(country_geom,aes(x=x,y=y,group=id,fill=shortfall_4C))+
  geom_polygon(colour="black")+  
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(-1,0.34))+
  labs(fill="Shortfall under 4°C warming")+
  theme+
  coord_equal()
plot(p4d)

p4 <- ggarrange(p4a,p4b,p4c,p4d)
plot(p4)

#########
country_geom$cap <- unlist(Country_summary_4C_invasive[match(country_geom$country,Country_summary_4C_invasive$Final_country),"cap"])
country_geom[country_geom$geom == 22 & (country_geom$y < 40 | country_geom$y > 55),"cap"] <- NA

p4c <- ggplot(country_geom,aes(x=x,y=y,group=id,fill=cap))+
  geom_polygon(colour="black")+  
  scale_fill_continuous(low="blue",high="red",na.value="white",limits=c(0,6))+
  labs(fill="Total capacity")+
  theme+
  coord_equal()
plot(p4c)



########
library(ggeffects)

predict_df <- ggemmeans(m_4C_country,terms="cap")
p4d <- ggplot(data=shortfall_country_4C,aes(x=cap,y=proj_diff))+
  geom_line(data=predict_df,aes(x=x,y=predicted))+
  geom_point()+
  ylab("Gain in alien ant species under 4°C warming")+
  xlab("Total capacity")+
  theme_bw()

plot(p4d)  
