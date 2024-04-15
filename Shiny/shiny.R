library(shiny)
library(tidyverse)
library(shinyWidgets)
library(DT)
library(leaflet)
library(rmapshaper)
library(shinydashboard)
library(leaflegend)
library(sf)
library(shinythemes)
library(htmlTable)
alldata <- read.csv("indoor_ant_analysis.csv")
alldata$label <- paste0(alldata$species,"_",alldata$polygon_name)
see <- alldata %>% group_by(label) %>% summarize(n=n()) %>% filter(n>1)
alldata$Status <- as.factor(alldata$Status)
#bentity.shp.sf <- st_read("bentity_shp_sf.shp")
#bentity.shp.sf <- ms_simplify(bentity.shp.sf,keep=0.05)

bentity.shp.sf <- st_read("simplified_bentity.shp")
bentity.shp.sf$polygon_name <- bentity.shp.sf$BENTITY
bentity.shp.sf$id <- 1:nrow(bentity.shp.sf)
bentity.shp.sf <- bentity.shp.sf %>% 
                    left_join(
                      distinct(alldata[alldata$Status != "Native",
                                       c("polygon_name",
                                         "Final_country",
                                         "threat",
                                         "IAS_list",
                                         "Existing_mgmt",
                                         "research",
                                         "outreach",
                                         "total_score")]),
                      by="polygon_name") 
  
alldata <- alldata %>% left_join(bentity.shp.sf[,c("polygon_name","id")],by="polygon_name")

all_harmful <- unique(alldata$species[alldata$Harmful == "Harmful"])
alldata[alldata$species %in% all_harmful,"Harmful"] <- "Harmful"

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

metadata <- data.frame(Variable=c("Species",
                                  "Region",
                                  "Date",
                                  "Strata",
                                  "Status",
                                  "E.Total",
                                  "S.Total",
                                  "Outdoor.Prob.Current",
                                  "Gain.Outdoor.Prob.2C",
                                  "Gain.Outdoor.Prob.4C"),
                       Definition = c("Species name",
                                    "Region name",
                                    "The year of the earliest occurrence record in the region",
                                    "Vertical habitat strata (six categories)",
                                    "Indoor (Only indoor records)/Naturalized (at least one outdoor record)",
                                    "Environmental Impact score of the species (the sum of six impacts)",
                                    "Socioeconomic impact score of the species (the sum of six impacts)",
                                    "Naturalization probability under current climate",
                                    "Gain in naturalization probablity under 2°C warming compared with current climate",
                                    "Gain in naturalization probablity under 4°C warming compared with current climate")
)

# ##
ui <-
  navbarPage("Indoor Ant Invasions", 
             theme = shinytheme("cosmo"),
             tabPanel("Map",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          radioButtons("Climate",
                                       "Climatic conditions",
                                       Climate,
                                       selected="Current"),
                          sliderInput("E.Impact",
                                      "Environmental Impact",
                                      min=0,
                                      max=14,
                                      value=c(0,14),
                                      step= 1),
                          sliderInput("S.Impact",
                                      "Socioeconomic Impact",
                                      min=0,
                                      max=16,
                                      value=c(0,16),
                                      step= 1),
                          virtualSelectInput("Species",
                                             "Species (Leave search empty and tick the box to select/deselect all)",
                                             choices=list("Harmful" = unique(alldata$species[alldata$Harmful == "Harmful"]),
                                                          "Non-harmful" = unique(alldata$species[alldata$Harmful != "Harmful"])),
                                             hideClearButton = F,
                                             selected = unique(alldata$species),
                                             search = T,
                                             multiple = T,
                                             optionsCount=5,
                                             updateOn = "close"),
                          # uiOutput("Species_UI"),
                          virtualSelectInput("Region",
                                             "Region (Leave search empty and tick the box to select/deselect all)",
                                             choices=sort(unique(alldata$polygon_name)),
                                             selected = sort(unique(alldata$polygon_name)),
                                             multiple = T,
                                             search = T,
                                             optionsCount=5,
                                             updateOn = "close"),
                          virtualSelectInput("Metric",
                                             "Metric",
                                             choices = c("Number of species","Environmental impact","Socioeconomic impact"),
                                             multiple = F,
                                             selected = "Number of species"),
                          virtualSelectInput("Environment",
                                             "Record",
                                             choices=c("Indoor records only","Outdoor records only","All records"),
                                             multiple = F,
                                             selected = "Indoor records only"),
                          actionButton("reset_input","Remove map-click filter")
                        ),
                        mainPanel(
                          width=9,
                          h6("Click polygon to view the data of specific region. 
                             Re-click the region or click the 'remove map-click filter' button to remove the filter"),
                          leafletOutput("climate_map"),
                          fluidRow(
                            tabsetPanel(
                              tabPanel("Data",
                                       tags$head(tags$style(HTML(css))),
                                       DTOutput("Dynamic")),
                              tabPanel("Current Impacts",
                                       plotOutput("Current_impact_plot")),
                              tabPanel("Warming Impacts (Region)",
                                       plotOutput("impact_region_plot")),
                              tabPanel("Warming Impacts (Species)",
                                       plotOutput("Species_impact_plot")
                              ),
                              tabPanel("Response Capacity",
                                       DTOutput("RC")),
                              type="pills"
                            )
                          )
                        )
                      )
             ),
             tabPanel("About",
                      h2("What is this?"),
                      h4("This is a shiny app visualizing ecological forecasts of how climate change facilitate indoor non-native ants spreading outdoors (i.e., naturalization). A relevant manuscript has  been submitted"),
                      h2("What do different columns and terms mean?"),
                      tableOutput("metadata"),
                      h2("What are the impacts?"),
                      h4("We included six environmental and six socioeconomic impacts as defined in the Gruber et al. (2022) paper (References provided below).
                         For example, environmental impacts include competition with other species, while socioeconomic impacts include crop loss.
                         Note that we are forecasting ant impacts in outdoor environments. Thus we did not consider records indicating records in indoor environments.
                         This includes occurrences of non-native ants in residential buildings as a nuisance, or contaminating food in storage.
                         Thus some species have impact records in the Gruber et al. database, but they are not listed as harmful here."),
                      h2("What are response capacities?"),
                      h4("This is extracted from a study evaluating CBD reports by different countries/adminstrative areas in 2015.
                         High response capacities indicate better invasive species management.
                         Note that the study evaluated the resposne capacity of different countries / adminstrative arears but not regions.
                         Thus, for example, all the states in the US have the same response capcities. 
                         However, the US outlying islands are considered as different unit due to very different management.
                         See the references listed below for details"),
                      h2("How reliable are the climate change forecasts?"),
                      h4("Our model has >90% accucracy predicting the status of each population under current climate.
                      Thus we assume it does a decent job predicting the future."),
                      h2("Why do you have climate change forecast for different types of records?"),
                      h4("In the manuscript, we only included the forecasts for indoor populations only. 
                      The cumulative gains in naturalization probability under warming are then inferred as gains in outdoor non-native species.
                      However, it's also possible to forecast the naturalization probability of outdoor populations. 
                      Imagine a naturalized population that have 50% naturalization probability under current climate, and warming increase the probability to 70%.
                      It's more tricky to interpret the cumulative gains when outdoor populations are included (since they have already naturalized, 
                      even though naturalization probability is not 100%). Thus the default settings are using indoor populations only
                      But we provide the options here in case you want to see it."),
                      h2("How updated are the data?"),
                      h4("The study integrated data from four databases, all obtained in 2023. Full references are listed below. Click the blue names to visit the websites hosting the datasets.
                         Note that the response capacity data was based on CBD report in 2015."),
                      h3(a("Ant occurrences", href="https://antmaps.org/")),
                      h4("Guénard, B., Weiser, M. D., Gomez, K., Narula, N., & Economo, E. P. (2017). The Global Ant Biodiversity Informatics (GABI) database: synthesizing data on the geographic distribution of ant species (Hymenoptera: Formicidae). Myrmecological News, 24, 83–89."),
                      h4("Wong, M. K., Economo, E. P., & Guénard, B. (2023). The global spread and invasion capacities of alien ants. Current Biology, 33(3), 566-571."),
                      h3(a("Vertical habitat strata", href="https://doi.org/10.6084/m9.figshare.21666191")),
                      h4("Wong, M. K., Economo, E. P., & Guénard, B. (2023). The global spread and invasion capacities of alien ants. Current Biology, 33(3), 566-571."),
                      h3(a("Climate", href="https://www.climatologylab.org/terraclimate.html")),
                      h4("Abatzoglou, J. T., Dobrowski, S. Z., Parks, S. A., & Hegewisch, K. C. (2018). TerraClimate, a high-resolution global dataset of monthly climate and climatic water balance from 1958–2015. Scientific Data, 5(1), 1–12."),
                      h3(a("Impact scores", href="https://datadryad.org/stash/dataset/doi:10.5061/dryad.2280gb5t2")),
                      h4("Gruber, M. A. M., Santoro, D., Cooling, M., Lester, P. J., Hoffmann, B. D., Boser, C., & Lach, L. (2022). A global review of socioeconomic and environmental impacts of ants reveals new insights for risk assessment. Ecological Applications, 32(4), e2577"),
                      h3(a("Response capacity"), href="https://www.fabiogeography.com/data"),
                      h4("Early, R., Bradley, B. A., Dukes, J. S., Lawler, J. J., Olden, J. D., Blumenthal, D. M., … Miller, L. P. (2016). Global threats from invasive alien species in the twenty-first century and national response capacities. Nature Communications, 7(1), 12485."),
                      h2("How to report a bug? Have questions? Want more features?"),
                      h4("Send an email to Toby P.N. Tsang (tpaknok@gmail.com)"),
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
  output$metadata <- renderTable(metadata)
  
  df <- reactive({
    record_df <- alldata %>% 
      filter(species %in% input$Species) %>%
      filter(polygon_name %in% input$Region) %>%
      filter((E.Total >= input$E.Impact[1] & E.Total <= input$E.Impact[2]) | Status == "Native") %>%
      filter((S.Total >= input$S.Impact[1] & S.Total <= input$S.Impact[2]) | Status == "Native")

    # if(!is.null(input$Species)) {
    #   record_df <- record_df %>% filter(species %in% input$Species)
    # } else {
    #   record_df <- record_df
    # }
    
    record_df <- switch(input$Environment,
                         "Indoor records only" = record_df %>% filter(num==0 | Status == "Native"),
                         "Outdoor records only" = record_df %>% filter(num==1 | Status == "Native"),
                         record_df)
  })
  
  # output$Species_UI <- renderUI({
  #   virtualSelectInput("Species",
  #                      "Species",
  #                      choices=list("Harmful" = unique(df()$species[df()$Harmful == "Harmful" & df()$Status != "Native"]),
  #                                   "Non-harmful" = unique(df()$species[df()$Harmful != "Harmful" & df()$Status != "Native"])),
  #                      selected = unique(df()$species[df()$Status != "Native"]),
  #                      search = T,
  #                      multiple = T)
  # })

  
  df_map <- reactive({
    var_name <- switch(
      paste0(input$Metric,"_",input$Climate),
      "Number of species_Current" = "dummy",
      "Number of species_2°C" = "proj_diff_2C",
      "Number of species_4°C" = "proj_diff_4C",
      "Environmental impact_Current" = "E.Total",
      "Environmental impact_2°C" = "Impact_2C.E.Total",
      "Environmental impact_4°C" = "Impact_4C.E.Total",
      "Socioeconomic impact_Current" = "S.Total",
      "Socioeconomic impact_2°C" = "Impact_2C.S.Total",
      "Socioeconomic impact_4°C" = "Impact_4C.S.Total")
      
      if (length(input$Species) == 1 & input$Climate == "Current") {
        df() %>% 
          mutate(Value = Status) %>%
          right_join(bentity.shp.sf[,c("polygon_name","BENTITY","id"),by="polygon_name"]) %>%
          st_as_sf() #Singe species mapping
      } else {
        df() %>% 
          group_by(polygon_name) %>% 
          summarize(Value = round(sum(get(var_name),na.rm=T),3)) %>%
          mutate(Value = na_if(Value,0)) %>%
          right_join(bentity.shp.sf[,c("polygon_name","BENTITY","id"),by="polygon_name"]) %>%
          st_as_sf()
      }
  })
  
  pal <- reactive(
    if (is.factor(df_map()$Value)) {
      colorFactor(c("#ff829b","#bd0026","Darkgreen"), levels = c("Indoor","Naturalized","Native"))
    } else {
      if (is.numeric(df_map()$Value) & length(na.omit(unique(df_map()$Value))) > 1) {
      colorNumeric("YlOrRd", domain = df_map()$Value)
    } else {
      colorFactor("#bd0026",domain = df_map()$Value)
      }
    }
  )
  
  legend_title <- reactive({
    if (length(input$Species) == 1) {
      switch(
        paste0(input$Metric,"_",input$Climate),
        "Number of species_Current" = "Populations",
        "Number of species_2°C" = "Projected gain in naturalization probability",
        "Number of species_4°C" = "Projected gain in naturalization probability",
        "Environmental impact_Current" = "Environmental impacts",
        "Environmental impact_2°C" = "Projected gain in environmental impacts",
        "Environmental impact_4°C" = "Projected gain in environmental impacts",
        "Socioeconomic impact_Current" = "Socioeconomic impacts",
        "Socioeconomic impact_2°C" = "Projected gain in socioeconomic impacts",
        "Socioeconomic impact_4°C" = "Projected gain in socioeconomic impacts")
    } else {
      switch(
        paste0(input$Metric,"_",input$Climate),
        "Number of species_Current" = "Number of non-native ant species",
        "Number of species_2°C" = "Projected gain in outdoor non-native ant species",
        "Number of species_4°C" = "Projected gain in outdoor non-native ant species",
        "Environmental impact_Current" = "Environmental impacts",
        "Environmental impact_2°C" = "Projected gain in environmental impacts",
        "Environmental impact_4°C" = "Projected gain in environmental impacts",
        "Socioeconomic impact_Current" = "Socioeconomic impacts",
        "Socioeconomic impact_2°C" = "Projected gain in socioeconomic impacts",
        "Socioeconomic impact_4°C" = "Projected gain in socioeconomic impacts")
    }
  })
  
  # legend_title <- reactive(
  #   if (length(input$Species) == 1 & input$Climate=="Current") {
  #     "Indoor populations"
  #   } else {
  #   if (length(input$Species) > 1 & input$Climate=="Current") {
  #     "Number of indoor non-native ant species"
  #   } else {
  #   if (length(input$Species) == 1 & input$Climate!="Current") {
  #     "Projected gain in naturalization probability"
  #   }
  #   else {"Projected gain in outdoor non-native ant species"
  #     }
  #   }
  #   }
  # )
  
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

        if (is.numeric(df_map()$Value) & length(na.omit(unique(df_map()$Value))) > 1) {
         map %>%  addLegendNumeric(pal=pal(),
                 values=~Value,
                 position = "bottomleft",
                 decreasing=F,
                 naLabel = "No record",
                 orientation = "horizontal",
                 width=100,
                 height=20,
                 title= legend_title())
        }
       else {
         map %>% addLegendFactor(pal=pal(),
                              values=~Value,
                              position = "bottomleft",
                              naLabel="No record",
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
  
   observe({
     req(input$reset_input)
     rv(NULL)
   })
  
  region_df <-reactive({
      if(is.null(rv())) { 
        df()
      } else { 
        df() %>% filter(id.y == rv())
      }
  })
  
  output$Dynamic <- renderDataTable({
    
    DataTable <- region_df() %>% 
      select("species",
             "polygon_name",
             "date",
             "sp.layer", 
             "Status",
             "E.Total",
             "S.Total",
             "current_status_projection",
             "proj_diff_2C",
             "proj_diff_4C") %>% 
      mutate(across(current_status_projection:proj_diff_4C, function(x) (round(x,3)))) %>%
      rename(Species = species,
             Region = polygon_name,
             Date = date,
             Strata = sp.layer,
             Outdoor.Prob.Current = current_status_projection,
             Gain.Outdoor.Prob.2C = proj_diff_2C ,
             Gain.Outdoor.Prob.4C = proj_diff_4C)
                                                 
     datatable(DataTable,
      extensions = "Buttons",
      options = list(
        paging=TRUE,
        dom = 'ftpB',
        buttons = list('copy', 'csv', 'excel')
      ),
      class = "display"
    )
  })
  
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
  
  impact_current <- reactive({
    Current_impact <- region_df() %>%  
      summarize(across(Animal.production:S.Total, function(x) sum(x,na.rm=T))) %>% 
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var="Impact")
    
    Current_impact$Type <- ifelse(Current_impact$Impact == "Animals" | Current_impact$Impact == "Competition" | Current_impact$Impact == "Diseases" | Current_impact$Impact == "Ecosystems" | Current_impact$Impact == "Hybridization" | Current_impact$Impact == "Plants", "Environmental","Socioeconomic")
    
    Current_impact <- subset(Current_impact, Impact != "E.Total" & Impact != "S.Total")
    
    Current_impact$Impact <- recode(Current_impact$Impact,Animal.production = "Animal Production")
    
    Current_impact$Impact <- factor(Current_impact$Impact,levels=c("Plants","Hybridization","Ecosystems","Diseases","Competition","Animals",
                                                                   "Social","Infrastructure","Health","Forestry","Crops","Animal Production"))
    
    return(Current_impact)
  })
  
  output$Current_impact_plot <- renderPlot({
    impact_current() %>% 
      ggplot(aes(y=Impact,x=V1)) +
      geom_bar(stat="identity",position = 'dodge',fill="#bd0026") +
      xlab("Cumulative scores under current climates") +
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
  })
  
  Species_impact <- reactive({
    Species_impact <- region_df() %>%  
      filter(Harmful == "Harmful") %>%
      group_by(species) %>%
      summarise(across(c("Impact_4C.E.Total","Impact_4C.S.Total","Impact_2C.E.Total","Impact_2C.S.Total","proj_diff_2C","proj_diff_4C"), ~mean(.x,na.rm=T)),
                n = n()) %>%
      pivot_longer(-c(species,n),values_to = "Score",names_to="Impact")
    
    Species_impact$Warming <- "4°C"
    Species_impact$Warming[grep("2C",Species_impact$Impact)] <- "2°C"
    
    Species_impact$Type <- rep(c("Environmental","Socioeconomic","Environmental","Socioeconomic","Naturalization probability","Naturalization probability"),
                                 nrow(Species_impact)/6)
    
    Species_impact$Warming <- factor(Species_impact$Warming,levels=c("4°C","2°C")) 
    Species_impact$Type <- factor(Species_impact$Type,levels=c("Naturalization probability","Environmental","Socioeconomic")) 
    
    Species_impact$species <- as.factor(Species_impact$species)
    return(Species_impact)
  })
  
  output$Species_impact_plot <- renderPlot({
    Species_plot_df <- Species_impact()
      
    if (nrow(Species_plot_df) == 0) {
      ggplot()+
        annotate(geom="text",x=0,y=0,label = "No harmful species \nbased on the filters",size=12)+
        theme_void()
    } else {
      ggplot(data=Species_plot_df, aes(y=species,x=Score,fill=Warming))+
      geom_bar(stat="Identity",position="dodge")+
      facet_wrap(~Type,scale="free_x")+
      xlab("Cumulative increases across populations under warming")+
      ylab("Species (Number of populations)")+
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
        select(!geometry & !Final_country)
      } else {
          bentity.shp.sf %>%
          filter(id == rv()) %>%
          select("polygon_name",Final_country:total_score) %>%
          as_tibble() %>%
          select(!geometry & !Final_country)
      }
  )
  
  output$RC <- renderDataTable(response_cap())
}
shinyApp(ui, server)
