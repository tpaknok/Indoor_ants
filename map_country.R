##########################

### map_country level

GBR_post <- subset(GBR, period == "1951-2014")

country <- ne_countries(scale=10,type="map_units")
country_df <- as.data.frame(country)

### YAR + YMD = YEM;
### WDEU + EDEU = DEU;
### ROM = ROU;
### USSR = LTU+GEO+EST+LVA+UKR+MDA+KGZ+UZB+TJK+ARM+AZE+TKM+RUS+BLR+KAZ;
### YUG = HRV+SVN+MKD+BIH+SRB+MNE
### ANT = ABW+CUW+BES+SXM; (only use data till 2010, as BES SXM no data. ARU became independent in 1986, merge with ANT)
### CZSK = CZE+SVK

GBR_post$revised_iso_o <- case_match(GBR_post$iso_o, 
                                     c("CZE","SVK")~"CZSK",
                                     c("HRV","SVN","MKD","BIH","SRB","MNE") ~ "YUG",
                                     c("ROM") ~ "ROU",
                                     c("WDEU","EDEU") ~ "DEU",
                                     c("YAR","YMD") ~ "YEM",
                                     c("LTU","GEO","EST","LVA","UKR","MDA","KGZ","UZB","TJK","ARM","AZE","TKM","RUS","BLR","KAZ") ~ "USSR",
                                     c("ABW","CUW","BES","SXM") ~ "ANT",
                                     .default=GBR_post$iso_o)

country$revised_iso_a3_eh <- case_match(country$iso_a3_eh, 
                                        c("CZE","SVK")~"CZSK",
                                        c("HRV","SVN","MKD","BIH","SRB","MNE") ~ "YUG",
                                        c("ROM") ~ "ROU",
                                        c("WDEU","EDEU") ~ "DEU",
                                        c("YAR","YMD") ~ "YEM",
                                        c("LTU","GEO","EST","LVA","UKR","MDA","KGZ","UZB","TJK","ARM","AZE","TKM","RUS","BLR","KAZ") ~ "USSR",
                                        c("-99") ~ NA,
                                        c("ABW","CUW","BES","SXM") ~ "ANT",
                                        .default=country$iso_a3_eh)


library(sf)
country.sf <- st_as_sf(country)
country.sf <- st_wrap_dateline(country.sf)
country.sf <- st_transform(country.sf, crs=st_crs("ESRI:54019"))
country.sf <- country.sf[country.sf$brk_name != "Antarctica",]

country.merged.sf <- st_as_sf(country.sf) %>% group_by(revised_iso_a3_eh) %>% summarize(n= n())

prelim_map <- ggplot(data=country.merged.sf,aes(fill=n))+
  geom_sf(colour="black",size=0.1)+
  labs(fill="")
plot(prelim_map)

remove(prelim_map)

GBR_post_flow <- GBR_post %>% 
  group_by(revised_iso_o,year) %>% 
  filter(revised_iso_o != "ANT" | (revised_iso_o == "ANT" & year <=2010)) %>%
  summarize(sum_flow = sum(FLOW,na.rm=T)) %>%
  ungroup %>%
  group_by(revised_iso_o) %>%
  summarize(average_flow = mean(sum_flow,na.rm=T))

country.merged.sf <- cbind(country.merged.sf,c(GBR_post_flow[match(country.merged.sf$revised_iso_a3_eh,GBR_post_flow$revised_iso_o),"average_flow"]))

prelim_map <- ggplot(data=country.merged.sf,aes(fill=log10(average_flow)))+
  geom_sf(colour="black",size=0.1)+
  labs(fill="")+
  scale_fill_gradient(low="#ffffb2",high="#bd0026",na.value="grey")
plot(prelim_map)

################ convert distributino from regional to country

list_df <- split(distribution.records %>% filter(Date < UK_first & UK_first > 1950 & UK_first != "not-given"),distribution.records %>% filter(Date < UK_first & UK_first > 1950 & UK_first != "not-given") %>% select("valid_species_name"))

list_shp <- lapply(list_df,function(x) bentity.shp.sf[bentity.shp.sf$BENTITY2_N %in% x[,"bentity2_name"],])
list_shp <- lapply(list_shp, st_centroid)

list_shp_country <- lapply(list_shp,function(x) st_intersection(x, country.merged.sf))
country_total <- c(table(unlist(lapply(list_shp_country, function(x) unique(x$revised_iso_a3_eh)))))

###
GBR_post_flow$distribution <- country_total[match(GBR_post_flow$revised_iso_o,names(country_total))]
GBR_post_flow$distribution[is.na(GBR_post_flow$distribution)] <- 0
plot(GBR_post_flow$average_flow,GBR_post_flow$distribution)

centroid<- as.data.frame(centroids(vect(st_transform(country.merged.sf, crs=st_crs("EPSG:4326")))),geom="XY") #st_centroid is probelmatic
GBR_post_flow <- cbind(GBR_post_flow,centroid[match(GBR_post_flow$revised_iso_o,centroid$revised_iso_a3_eh),c("x","y")])

library(spaMM)

m <- fitme(distribution~average_flow+Matern(1|x+y),data=GBR_post_flow,family="poisson")
summary(m)
