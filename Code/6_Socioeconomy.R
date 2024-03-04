library(tidyverse)
library(ggplot2)
worldbank <- read.csv("Data/socioeconomy.csv")

worldbank_data <- worldbank %>% na.omit() 

cor_df <- worldbank_data %>% 
          group_by(Series.Name) %>% 
          summarize(r=cor(YR2015,`YR2022.YR2021`),n=n())

cor_df$label <- paste0("N=",cor_df$n,"\nPearson's r=",round(cor_df$r,3))

p7 <- ggplot(worldbank_data,aes(y=YR2022.YR2021,x=YR2015))+
      geom_point()+
      geom_abline(slope=1)+
      geom_text(data=cor_df,aes(label=label),x=+Inf,y=-Inf,vjust=-0.25,hjust=1,size=3)+
      labs(y="Year 2021/2022",x = "Year 2015")+
      geom_smooth(method="lm")+
      facet_wrap(~Series.Name,scales="free")+
      theme_classic()
plot(p7)

ggsave("Figures/socioeconomy.tiff",dpi=600,compression="lzw",width=16,height = 16,units="cm")
