### excel file

### regional sheet
S6a <- bentity.shp.sf %>%   
  as_tibble %>%
  as.data.frame() %>%
  select(Region_name=BENTITY2_N,
         Observed_indoor_richness=indoor.sr,
         Observed_naturalized_richness=outdoor.sr,
         Projected_naturalized_richness_indoor_current=proj_indoor_Current_total,
         Projected_naturalized_richness_indoor_2C=proj_indoor_2C_total,
         Projected_naturalized_richness_indoor_4C=proj_indoor_4C_total,
         Projected_naturalized_richness_difference_2C=proj_diff_indoor_2C_net,
         Projected_naturalized_richness_difference_4C=proj_diff_indoor_4C_net,
         Projected_naturalized_richness_difference_warm_scen = warming_diff_indoor_net,
         Observed_harmful_indoor_richness=indoor.sr.Harmful,
         Observed_harmful_naturalized_richness=outdoor.sr.Harmful,
         Projected_harmful_naturalized_richness_indoor_current=proj_Harmful_indoor_Current_total,
         Projected_harmful_naturalized_richness_indoor_2C=proj_Harmful_indoor_2C_total,
         Projected_harmful_naturalized_richness_indoor_4C= proj_Harmful_indoor_4C_total,
         Projected_harmful_naturalized_richness_difference_2C=proj_diff_Harmful_indoor_2C_net,
         Projected_harmful_naturalized_richness_difference_4C=proj_diff_Harmful_indoor_4C_net,
         Projected_harmful_naturalized_richness_warm_scen=warming_diff_Harmful_indoor_net,
         Observed_env_impact_indoor=current_impact_E_indoor,
         Observed_env_impact_naturalized=current_impact_E_outdoor,
         Cumulative_env_SWIS_indoor_current=proj_impact_E_Current_total,
         Cumulative_env_SWIS_indoor_2C=proj_impact_E_2C_total,
         Cumulative_env_SWIS_indoor_4C=proj_impact_E_4C_total,
         Projected_env_SWIS_difference_2C=proj_diff_impact_E_2C_net,
         Projected_env_SWIS_difference_4C=proj_diff_impact_E_4C_net,
         Projected_env_SWIS_difference_warm_scen=warming_diff_impact_E,
         Observed_soc_impact_indoor=current_impact_S_indoor,
         Observed_soc_impact_naturalized=current_impact_S_outdoor,
         Cumulative_soc_SWIS_indoor_current=proj_impact_S_Current_total,
         Cumulative_soc_SWIS_indoor_2C=proj_impact_S_2C_total,
         Cumulative_soc_SWIS_indoor_4C=proj_impact_S_4C_total,
         Projected_soc_SWIS_difference_2C=proj_diff_impact_S_2C_net,
         Projected_soc_SWIS_difference_4C=proj_diff_impact_S_4C_net,
         Projected_soc_SWIS_difference_warm_scen=warming_diff_impact_S,
  )

pos <- which(is.na(S6a$Observed_indoor_richness))
S6a[pos,c("Projected_naturalized_richness_indoor_current","Projected_naturalized_richness_indoor_2C","Projected_naturalized_richness_indoor_4C")] <- NA

pos <- which(is.na(S6a$Observed_harmful_indoor_richness))
S6a[pos,c("Projected_harmful_naturalized_richness_indoor_current","Projected_harmful_naturalized_richness_indoor_2C","Projected_harmful_naturalized_richness_indoor_4C")] <- NA

### country

sp_country_table <- NMI_analysis %>%  
  group_by(species,Final_country,Harmful,total_score,num) %>%  
  summarise(current_status_projection = max(current_status_projection),
            future_status_2C = max(future_status_2C), #note that max of proj_diff_2C might not be the same as future_status_2C - current_status_projection (if different populations were selected)
            future_status_4C = max(future_status_4C), #note that max of proj_diff_4C might not be the same as future_status_4C - current_status_projection (if different populations were selected)
            proj_diff_2C = max(proj_diff_2C),
            proj_diff_4C = max(proj_diff_4C),
            observed_Impact_current.E.Total = max(observed_Impact_current.E.Total), 
            projected_Impact_current.E.Total=max(projected_Impact_current.E.Total), 
            projected_Impact_2C.E.Total = max(projected_Impact_2C.E.Total),
            projected_Impact_4C.E.Total = max(projected_Impact_4C.E.Total),
            proj_impact_diff_2C.E.Total = max(proj_impact_diff_2C.E.Total), #might not euqal max(2C)-max(current) because of different indoor-restricted populations in the countries
            proj_impact_diff_4C.E.Total = max(proj_impact_diff_4C.E.Total),
            observed_Impact_current.S.Total = max(observed_Impact_current.S.Total),
            projected_Impact_current.S.Total=max(projected_Impact_current.S.Total),
            projected_Impact_2C.S.Total = max(projected_Impact_2C.S.Total),
            projected_Impact_4C.S.Total = max(projected_Impact_4C.S.Total),
            proj_impact_diff_2C.S.Total = max(proj_impact_diff_2C.S.Total),
            proj_impact_diff_4C.S.Total = max(proj_impact_diff_4C.S.Total)) %>%
  mutate(max.num = max(num)) #indoor restricted at country levels?

S6b <- sp_country_table %>%   
  group_by(Final_country) %>%
  summarize(Observed_indoor_richness=n_distinct(species[max.num==0]), #only consider indoor-restricted species at country levels
            Observed_naturalized_richness=n_distinct(species[max.num==1]),
            Projected_naturalized_richness_indoor_current=sum(current_status_projection[max.num==0]), 
            Projected_naturalized_richness_indoor_2C=sum(future_status_2C[max.num==0]),
            Projected_naturalized_richness_indoor_4C=sum(future_status_4C[max.num==0]),
            Projected_naturalized_richness_difference_2C=sum(proj_diff_2C[num==0]), # use num not max num.This will also consider indoor-restricted species projected to expand distribution, even though they may have already naturalized in other parts of the country
            Projected_naturalized_richness_difference_4C=sum(proj_diff_4C[num==0]),
            Projected_naturalized_richness_difference_warm_scen = Projected_naturalized_richness_difference_4C-Projected_naturalized_richness_difference_2C,
            Observed_harmful_indoor_richness=n_distinct(species[max.num==0 & Harmful == "Harmful"]),
            Observed_harmful_naturalized_richness=n_distinct(species[max.num==1 & Harmful == "Harmful"]),
            Projected_harmful_naturalized_richness_indoor_current=sum(current_status_projection[max.num==0 & Harmful == "Harmful"]),
            Projected_harmful_naturalized_richness_indoor_2C=sum(future_status_2C[max.num==0 & Harmful == "Harmful"]),
            Projected_harmful_naturalized_richness_indoor_4C=sum(future_status_4C[max.num==0 & Harmful == "Harmful"]),
            Projected_harmful_naturalized_richness_difference_2C=sum(proj_diff_2C[num==0 & Harmful == "Harmful"]),
            Projected_harmful_naturalized_richness_difference_4C=sum(proj_diff_2C[num==0 & Harmful == "Harmful"]),
            Projected_harmful_naturalized_richness_warm_scen=Projected_harmful_naturalized_richness_difference_4C-Projected_harmful_naturalized_richness_difference_2C,
            Observed_env_impact_indoor=sum(observed_Impact_current.E.Total[max.num==0]),
            Observed_env_impact_naturalized=sum(observed_Impact_current.E.Total[max.num==1]),
            Cumulative_env_SWIS_indoor_current=sum(projected_Impact_current.E.Total[max.num==0]),
            Cumulative_env_SWIS_indoor_2C=sum(projected_Impact_2C.E.Total[max.num==0]),
            Cumulative_env_SWIS_indoor_4C=sum(projected_Impact_4C.E.Total[max.num==0]),
            Projected_env_SWIS_difference_2C=sum(proj_impact_diff_2C.E.Total[num==0]),
            Projected_env_SWIS_difference_4C=sum(proj_impact_diff_4C.E.Total[num==0]),
            Projected_env_SWIS_difference_warm_scen=Projected_env_SWIS_difference_4C- Projected_env_SWIS_difference_2C,
            Observed_soc_impact_indoor=sum(observed_Impact_current.S.Total[max.num==0]),
            Observed_soc_impact_naturalized=sum(observed_Impact_current.S.Total[max.num==1]),
            Cumulative_soc_SWIS_indoor_current=sum(projected_Impact_current.S.Total[max.num==0]),
            Cumulative_soc_SWIS_indoor_2C=sum(projected_Impact_2C.S.Total[max.num==0]),
            Cumulative_soc_SWIS_indoor_4C=sum(projected_Impact_4C.S.Total[max.num==0]),
            Projected_soc_SWIS_difference_2C=sum(proj_impact_diff_2C.S.Total[num==0]),
            Projected_soc_SWIS_difference_4C=sum(proj_impact_diff_4C.S.Total[num==0]),
            Projected_soc_SWIS_difference_warm_scen=Projected_soc_SWIS_difference_4C- Projected_soc_SWIS_difference_2C,
  ) %>%
  left_join(capacity[,c("Country_global_modified","threat","IAS_list","Existing_mgmt","research","outreach")],by=join_by(Final_country == Country_global_modified)) %>%
  filter(!is.na(Final_country))

write.csv(S6a,"Results/S6a.csv")
write.csv(S6b,"Results/S6b.csv")

###
 