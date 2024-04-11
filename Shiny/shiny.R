library(shiny)
library(tidyverse)
library(shinyWidgets)
library(DT)
library(leaflet)
library(rmapshaper)
library(shinydashboard)
library(leaflegend)
library(sf)
alldata <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/R/shiny/indoor_ants_analysis.csv")
indoor_ants <- alldata %>% filter(num == 0)
indoor_ants$num <- 1

bentity.shp.sf <- st_read("C:/Users/pakno/OneDrive - University of Toronto/Indoor_ants/R/shiny/bentity_shp_sf.shp")
bentity.shp.sf$polygon_name <- bentity.shp.sf$BENTITY
bentity.shp.sf <- ms_simplify(bentity.shp.sf,keep=0.1)
bentity.shp.sf$id <- 1:nrow(bentity.shp.sf)
bentity.shp.sf <- bentity.shp.sf %>% 
                    left_join(
                      distinct(alldata[,c("polygon_name",
                                           "Final_country",
                                           "threat",
                                           "IAS_list",
                                           "Existing_mgmt",
                                           "research",
                                           "outreach",
                                           "total_score")]),
                      by="polygon_name") 
  
indoor_ants <- indoor_ants %>% left_join(bentity.shp.sf[,c("polygon_name","id")],by="polygon_name")

impact_names <- rep(c("Animal Production",
                       "Animals",
                       "Competition",
                       "Crops",
                       "Diseases",
                       "Ecosystems",
                       "Forestry",
                       "Health",
                       "Hybridization",
                       "Infrastructure",
                       "Plants",
                       "Social",
                       "E.Total",
                       "S.Total"),2)

impact_display <- c("E.Total", "S.Total")
# theme <- theme(axis.line=element_blank(),
#                axis.text.x=element_blank(),
#                axis.text.y=element_blank(),
#                axis.ticks=element_blank(),
#                axis.title.x=element_blank(),
#                axis.title.y=element_blank(),
#                panel.background=element_rect(colour="white",fill="white"),
#                panel.border=element_blank(),
#                panel.grid.major=element_blank(),
#                panel.grid.minor=element_blank(),
#                plot.background=element_rect(colour="white",fill="white"),
#                legend.position="bottom")

Climate <- c("Current","2°C","4°C")

css <- '.nav-tabs>li>a {
  font-family: "Lucida Sans", sans-serif;
  color: black;
}'

# ##
ui <-
  dashboardPage(
  dashboardHeader(title= "Indoor Ant Invasions under Climate Change"),
  dashboardSidebar(
    virtualSelectInput("Species",
                           "Species",
                           choices=list("Harmful" = unique(indoor_ants$species[indoor_ants$Harmful == "Harmful"]),
                                        "Non-harmful" = unique(indoor_ants$species[indoor_ants$Harmful != "Harmful"])),
                           #options = list(`actions-box` = TRUE),
                           selected = unique(indoor_ants$species),
                           multiple = T),
    radioButtons("Climate",
                     "Climatic conditinos",
                     Climate,
                     selected="Current")
  ),
  dashboardBody(
    leafletOutput("climate_map"),
    fluidRow(
      tabsetPanel(
        tabPanel("Data",
                 tags$head(tags$style(HTML(css))),
                 DTOutput("Dynamic")),
        tabPanel("Detailed Impact",
                  plotOutput("impact_region_plot")),
        tabPanel("Species impact",
                 plotOutput("Species_impact_plot")
                 ),
        tabPanel("Response Capacity",
                 DTOutput("RC")),
        type="pills"
      )
    )
  )
)

##
# ui <- 
#     dashboardPage(
#     dashboardHeader(title= "Indoor Ant Invasions under Climate Change"),
#     dashboardSidebar(
#       virtualSelectInput("Species",
#                              "Species",
#                              choices=list("Harmful" = unique(indoor_ants$species[indoor_ants$Harmful == "Harmful"]),
#                                           "Non-harmful" = unique(indoor_ants$species[indoor_ants$Harmful != "Harmful"])),
#                              #options = list(`actions-box` = TRUE),
#                              selected = unique(indoor_ants$species),
#                              multiple = T),
#       radioButtons("Climate",
#                        "Climatic conditinos",
#                        Climate,
#                        selected="Current")
#     ),
#     dashboardBody(
#       tags$style(type = "text/css", "#climate_map {height: calc(100vh - 80px) !important;}"),
#       leafletOutput("climate_map",width = "100%", height = "100%")
#     )
#   )
  
server <- function(input, output, session) {
  df <- reactive(indoor_ants %>% 
                   filter(species %in% input$Species))
  
  df_map <- reactive(
    df() %>%
    group_by(polygon_name) %>% 
    summarize(Value = ifelse(input$Climate == "Current", 
                             sum(num), 
                             ifelse(input$Climate == "2°C", 
                                    round(sum(proj_diff_2C),3),
                                    round(sum(proj_diff_4C),3)))) %>%
    right_join(bentity.shp.sf[,c("polygon_name","BENTITY","id"),by="polygon_name"]) %>%
    st_as_sf() 
    )
  
  pal <- reactive(
    if (length(na.omit(unique(df_map()$Value))) > 1) colorNumeric("YlOrRd", domain = df_map()$Value)
    else  colorFactor("#bd0026", domain = df_map()$Value)
    )
  
  legend_title <- reactive(
    if (length(input$Species) == 1 & input$Climate=="Current") {
      "Indoor populations"
    } else {
    if (length(input$Species) > 1 & input$Climate=="Current") {
      "Number of indoor non-native ant species"
    } else {
    if (length(input$Species) == 1 & input$Climate!="Current") {
      "Projected gain in naturalization probability"
    }
    else {"Projected gain in outdoor non-native ant species"
      }
    }
    }
  )
  
  output$climate_map <- renderLeaflet({
     map <- df_map() %>%
     leaflet() %>% 
     addTiles(options = providerTileOptions(minZoom = 1)) %>%
     addPolygons(fillColor = ~pal()(Value),
                  color = "black",
                  fillOpacity = 0.7,
                  weight=0.2,
                  popup = paste0(df_map()$polygon_name,":",df_map()$Value),
                  layerId = ~id)

        if (length(unique(na.omit(df_map()$Value))) > 1) {
         map %>%  addLegendNumeric(pal=pal(),
                 values=~Value,
                 position = "bottomleft",
                 decreasing=F,
                 naLabel = "No indoor record",
                 orientation = "horizontal",
                 width=100,
                 height=20,
                 title= legend_title())
        }
       else {
         map %>% addLegendFactor(pal=pal(),
                              values=~Value,
                              position = "bottomleft",
                              naLabel="No indoor record",
                              width=100,
                              height=20,
                              title= legend_title())
       }
   })
   
  rv <- reactiveVal()
  
  observeEvent(input$climate_map_shape_click, {
    if(!is.null(rv()) && rv() == input$climate_map_shape_click$id)
      rv(NULL)     # Reset filter
    else
      rv(input$climate_map_shape_click$id)
    })
  
  region_df <-reactive({
      if(is.null(rv())) { 
        df()
      } else { 
        df() %>% filter(id.y == rv())
      }
  })
  
  output$Dynamic <- renderDataTable({
    datatable(region_df() %>% select("species","polygon_name","date","sp.layer"))
  }
 )
  
  impact_region <- reactive({
    Global_impact <- region_df() %>%  
      group_by(polygon_name) %>%
      summarize(across(Impact_4C.Animal.production:Impact_2C.S.Total, function(x) sum(x,na.rm=T))) %>%
      ungroup() %>%
      summarize(across(Impact_4C.Animal.production:Impact_2C.S.Total,list(sum=~sum(.x,na.rm=T),
                                                                      n=~length(.x[.x>0]))))
    
    Global_impact <- data.frame(Score=t(Global_impact[seq(1,56,by=2)]),
                                n=t(Global_impact[seq(2,56,by=2)]),
                                Warming = c(rep("4°C",14),rep("2°C",14)),
                                Impact = rep(impact_names,2))
    
    Global_impact$Type <- ifelse(Global_impact$Impact == "Animals" | Global_impact$Impact == "Competition" | Global_impact$Impact == "Diseases" | Global_impact$Impact == "Ecosystems" | Global_impact$Impact == "Hybridization" | Global_impact$Impact == "Plants", "Environmental","Socioeconomic")
    
    Global_impact <- subset(Global_impact, Impact != "E.Total" & Impact != "S.Total")
    
    Global_impact$Impact <- recode(Global_impact$Impact,Animal.production = "Animal Production")
    
    Global_impact$Impact <- factor(Global_impact$Impact,levels=c("Plants","Hybridization","Ecosystems","Diseases","Competition","Animals",
                                                                 "Social","Infrastructure","Health","Forestry","Crops","Animal Production"))
    
    Global_impact$Warming <- factor(Global_impact$Warming,levels=c("4°C","2°C")) 
    
    return(Global_impact)
  })
  
  output$impact_region_df <- renderDataTable(impact_region())
  
  output$impact_region_plot <- renderPlot({
      impact_region() %>% 
      ggplot(aes(label=n,y=Impact,x=Score,fill=Warming)) +
      geom_bar(stat="identity",position = 'dodge') +
      scale_fill_manual(values=c("#bd0026","#feb24c"))+
      xlab("Cumulative score increases under warming") +
      ylab("") +
      facet_wrap(~Type,scales="free_y") +
      theme_bw() +
      theme(legend.position="bottom",
            axis.text.x=element_text(size=12),
            axis.text.y=element_text(size=12),
            axis.title.x=element_text(size=16),
            axis.title.y=element_text(size=16),
            legend.text=element_text(size=12),
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.5, "cm"),
            legend.title = element_text(size=12),
            strip.text = element_text(size = 16))
    }
  )
  
  Species_impact <- reactive({
    Species_impact <- region_df() %>%  
      filter(Harmful == "Harmful") %>%
      group_by(species) %>%
      summarize(across(Impact_4C.Animal.production:Impact_2C.S.Total, function(x) sum(x,na.rm=T)),
                n = n()) %>%
      pivot_longer(-c(species,n),values_to = "Score",names_to="Impact")
    
    Species_impact$Warming <- case_match(str_match(Species_impact$Impact, "Impact_(.*?)\\.")[,2],
                                    "4C" ~ "4°C",
                                    "2C" ~ "2°C",
                                    "Current" ~ "Current")
    
    Species_impact$Impact <- rep(impact_names,nrow(Species_impact)/length(impact_names))
    Species_impact$Type <- ifelse(Species_impact$Impact == "Animals" | 
                                    Species_impact$Impact == "Competition" | 
                                    Species_impact$Impact == "Diseases" | 
                                    Species_impact$Impact == "Ecosystems" | 
                                    Species_impact$Impact == "Hybridization" | 
                                    Species_impact$Impact == "Plants" | 
                                    Species_impact$Impact == "E.Total", 
                                  "Environmental",
                                  "Socioeconomic")
    
    Species_impact$Impact <- recode(Species_impact$Impact,Animal.production = "Animal Production")
    
    Species_impact$Impact <- factor(Species_impact$Impact,levels=c("E.Total","Plants","Hybridization","Ecosystems","Diseases","Competition","Animals",
                                                                 "S.Total","Social","Infrastructure","Health","Forestry","Crops","Animal Production"))
    
    Species_impact$Warming <- factor(Species_impact$Warming,levels=c("4°C","2°C")) 
    
    Species_impact$species <- as.factor(Species_impact$species)
    return(Species_impact)
  })
  
  output$Species_impact_plot <- renderPlot({
    Species_plot_df <- Species_impact() %>% 
      filter(Impact %in% impact_display)
      
    if (nrow(Species_plot_df) == 0) {
      ggplot()+
        annotate(geom="text",x=0,y=0,label = "No harmful species restricted to indoor environments \nbased on the filters",size=12)+
        theme_void()
    } else {
      ggplot(data=Species_plot_df, aes(y=species,x=Score,fill=Warming))+
      geom_bar(stat="Identity",position="dodge")+
      facet_wrap(~Type,scale="free_x")+
      xlab("Cumulative increases across indoor populations under warming")+
      ylab("Species (Number of indoor populations)")+
      scale_fill_manual(values=c("#bd0026","#feb24c"))+
      scale_y_discrete(labels = paste0(gsub("\\."," ",levels(Species_plot_df$species))," (",Species_plot_df$n[!duplicated(Species_plot_df$species)],")"))+
      theme_bw()+
      theme(legend.position="bottom",
              axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.x=element_text(size=16),
              axis.title.y=element_text(size=16),
              legend.text=element_text(size=12),
              legend.key.height = unit(0.5, "cm"),
              legend.key.width = unit(0.5, "cm"),
              legend.title = element_text(size=12),
              strip.text = element_text(size = 16))
    }
  })
  
  response_cap <- reactive(
    if(is.null(rv())) {
      bentity.shp.sf %>%
        select("polygon_name",Final_country:total_score) %>%
        as_tibble() %>%
        select(!geometry)
      } else {
          bentity.shp.sf %>%
          filter(id == rv()) %>%
          select("polygon_name",Final_country:total_score) %>%
          as_tibble() %>%
          select(!geometry)
      }
  )
  
  output$RC <- renderDataTable(response_cap())
}
shinyApp(ui, server)
