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
                      alien.gain = site_summary$proj_diff,
                      alien.gain.4C = site_summary_4C$proj_diff)

gain_df <- data.frame(gain_df,
                 invasive.gain = site_summary_invasive[match(gain_df$shp,site_summary_invasive$polygon_name),"proj_diff"],
                 invasive.gain.4C = site_summary_invasive_4C[match(gain_df$shp,site_summary_invasive_4C$polygon_name),"proj_diff"])

colnames(gain_df) <- c("shp","shp.ID","sr","indoor.sr","alien.gain","alien.gain.4C","invasive.gain","invasive.gain.4C")

gain_df$Country <- bentity.shp.df[match(gain_df$shp,bentity.shp.df$BENTITY2_N),"Final_country"]

gain_df$cap.R <- capacity[match(gain_df$Country,capacity$Country_global_modified),"reactive"]
gain_df$cap.P <- capacity[match(gain_df$Country,capacity$Country_global_modified),"proactive"]

gain_df$cap <- gain_df$cap.R+gain_df$cap.P

gain_df$temp_pca_change <- future_df[match(gain_df$shp,future_df$polygon_name),"scaled.temp_pca"]-future_df[match(gain_df$shp,future_df$polygon_name),"scaled.temp_pca_current"]

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
library(quantreg)
min_max <- function (x) {
  scaled.var <- (x-min(x))/(max(x)-min(x))
  return(scaled.var)
}

avg_country <- aggregate(cbind(invasive.gain,invasive.gain.4C)~Country,data=gain_df,FUN=max,na.rm=T)
Country_summary_2C_invasive$max.gain.2C <- avg_country[match(Country_summary_2C_invasive$Final_country,avg_country$Country),"invasive.gain"]
Country_summary_2C_invasive$cap.R <- capacity[match(Country_summary_2C_invasive$Final_country,capacity$Country_global),"reactive"]
Country_summary_2C_invasive$cap.P <- capacity[match(Country_summary_2C_invasive$Final_country,capacity$Country_global),"proactive"]
Country_summary_2C_invasive$cap <- Country_summary_2C_invasive$cap.P+Country_summary_2C_invasive$cap.R
Country_summary_2C_invasive <- cbind(Country_summary_2C_invasive,capacity[match(Country_summary_2C_invasive$Final_country,capacity$Country_global),7:12])
Country_summary_2C_invasive$GDP <- GDP_dat[match(Country_summary_2C_invasive$Final_country,GDP_dat$country),"NY.GDP.MKTP.KD"]
shortfall_country_2C <- Country_summary_2C_invasive[!is.na(Country_summary_2C_invasive$cap),]
shortfall_country_2C$warming <- "2C"

m_2C_rq <- rq(max.gain.2C~log(GDP),data=Country_summary_2C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_2C_rq,se="boot",R=1000)

m_2C_rq <- rq(max.gain.2C~outreach,data=Country_summary_2C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_2C_rq,se="boot",R=1000)

m_2C_rq <- rq(max.gain.2C~research,data=Country_summary_2C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_2C_rq,se="boot",R=1000)

m_2C_rq <- rq(max.gain.2C~Existing_mgmt,data=Country_summary_2C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_2C_rq,se="boot",R=1000)

m_2C_rq <- rq(max.gain.2C~IAS_list,data=Country_summary_2C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_2C_rq,se="boot",R=1000)

m_2C_rq <- rq(max.gain.2C~threat,data=Country_summary_2C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_2C_rq,se="boot",R=1000)

Country_summary_4C_invasive$max.gain.4C <- avg_country[match(Country_summary_4C_invasive$Final_country,avg_country$Country),"invasive.gain.4C"]
Country_summary_4C_invasive$cap.R <- capacity[match(Country_summary_4C_invasive$Final_country,capacity$Country_global),"reactive"]
Country_summary_4C_invasive$cap.P <- capacity[match(Country_summary_4C_invasive$Final_country,capacity$Country_global),"proactive"]
Country_summary_4C_invasive$cap <- Country_summary_4C_invasive$cap.P+Country_summary_4C_invasive$cap.R
Country_summary_4C_invasive <- cbind(Country_summary_4C_invasive,capacity[match(Country_summary_4C_invasive$Final_country,capacity$Country_global),7:12])
Country_summary_4C_invasive$GDP <- GDP_dat[match(Country_summary_4C_invasive$Final_country,GDP_dat$country),"NY.GDP.MKTP.KD"]
shortfall_country_4C <- Country_summary_4C_invasive[!is.na(Country_summary_4C_invasive$cap),]
shortfall_country_4C$warming <- "4C"

m_4C_rq <- rq(max.gain.4C~log(GDP),data=Country_summary_4C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_4C_rq,se="boot",R=1000)

m_4C_rq <- rq(max.gain.4C~outreach,data=Country_summary_4C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_4C_rq,se="boot",R=1000)

m_4C_rq <- rq(max.gain.4C~research,data=Country_summary_4C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_4C_rq,se="boot",R=1000)

m_4C_rq <- rq(max.gain.4C~Existing_mgmt,data=Country_summary_4C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_4C_rq,se="boot",R=1000)

m_4C_rq <- rq(max.gain.4C~IAS_list,data=Country_summary_4C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_4C_rq,se="boot",R=1000)

m_4C_rq <- rq(max.gain.4C~threat,data=Country_summary_4C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_4C_rq,se="boot",R=1000)

colSums(na.omit(Country_summary_4C_invasive)[order(na.omit(Country_summary_4C_invasive)$max.gain.4C,decreasing=T)[1:30],9:14])/30
colSums(na.omit(Country_summary_2C_invasive)[order(na.omit(Country_summary_2C_invasive)$max.gain.2C,decreasing=T)[1:30],9:14])/30

m_4C_rq <- rq(Country_summary_4C_invasive$max.gain.4C-Country_summary_2C_invasive$max.gain.2C~log(GDP),data=Country_summary_4C_invasive,tau=c(0.5,0.7,0.9),method="fn")
summary(m_4C_rq,se="boot",R=1000)

#############################
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
