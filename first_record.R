set <- exotic.records %>% 
  filter(Improved.Classification == "Exotic" | Improved.Classification == "Indoor Introduced") %>% 
  group_by(valid_species_name,bentity2_name,Improved.Classification) %>% 
  arrange(Date) %>% 
  slice(1L) %>%
  ungroup() %>%
  group_by(pop) %>%
  summarize(diff = Date[Improved.Classification == "Indoor Introduced"] - Date[Improved.Classification == "Exotic"])
