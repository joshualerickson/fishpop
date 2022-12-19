library(tidyverse)
library(sf)

pat_fish <- readxl::read_xlsx("C:/Users/joshualerickson/Box/EUR/2620FishPlng/data by waterbody/District fish hab data compilation.xlsx")

#split data into small and large based on mean width
#this helps stratify the streams by size

pat_fish <- pat_fish %>% mutate(rowID = str_remove_all(`Stream ID`,"-"),
                                    rowID = str_sub(rowID, end = -5))

pat_fish <- pat_fish %>% mutate(size = ifelse(`Mean Width (ft)` >= 15, paste0("Large"), paste0("Small")))

pat_fish <- pat_fish %>% mutate(`INFS LWD/mile` = str_replace(`INFS LWD/mile`, '-', NA_character_),
                                `INFS LWD/mile` = as.numeric(`INFS LWD/mile`)) %>% mutate(across(is.numeric, round, 2),across(is.character, factor)) 

pat_fish <- pat_fish %>%  
  mutate(`Fish Species` = ifelse(`Fish Species` == "EBT, RBT", paste('EBT,RBT'), paste(`Fish Species`)))

pat_fish <- pat_fish %>% janitor::remove_empty(which = 'cols')

pat_fish <- pat_fish %>% mutate(stream_segment = str_remove_all(`Stream ID`, '-[^-]+$'))

pat_fish <- pat_fish %>% janitor::clean_names()

pat_fish_sf <- sf::read_sf('data/eda_fish_spatial.gpkg', layer = 'eda_fish_spatial') %>% 
               filter(strm_sg != 'Fortine 2')

pat_fish_df <- pat_fish_sf %>% st_drop_geometry() %>% 
  select(starts_with('crop'), 'stream_segment' = strm_sg) %>% 
  right_join(pat_fish, by = 'stream_segment')

pat_fish_df <- pat_fish_df %>% relocate(starts_with('crop'), .after = everything())

#for all the data ts
pat_fish_ts <- pat_fish_df %>% janitor::clean_names() %>% 
  mutate(
    area_acres = crop_fac_taudem_17all_int*30*30*0.000247105,
    silt = as.numeric(crop_silt_100_cpg_all),
    silt = as.numeric(crop_sand_100_cpg_all),
    permeability = as.numeric(crop_permeable_100_cpg_all),
    precip = crop_us_precip_1981_2010_cpg_all,
    temp = crop_us_tmax_1981_2010_int_cpg_all) %>% 
  select(stream_segment:reach,year, mean_width_ft:fish_species, area_acres:temp)

write_csv(pat_fish_ts, 'data/fish_data_all.csv')

#for mean aggregation
pat_fish_df <- pat_fish_df %>% 
  group_by(stream_segment) %>% 
  mutate(across(mean_width_ft:psd150mm, ~mean(.x, na.rm = TRUE), .names = 'mean_{.col}')) %>% ungroup()

pat_fish_df <- pat_fish_df %>% group_by(stream_segment) %>% slice(n=1) %>% ungroup()
fish_stuff2 <- pat_fish_df %>% janitor::clean_names() %>% 
  mutate(
         area_acres = crop_fac_taudem_17all_int*30*30*0.000247105,
         silt = as.numeric(crop_silt_100_cpg_all),
         silt = as.numeric(crop_sand_100_cpg_all),
         permeability = as.numeric(crop_permeable_100_cpg_all),
         precip = crop_us_precip_1981_2010_cpg_all,
         temp = crop_us_tmax_1981_2010_int_cpg_all) %>% 
  select(stream_segment:reach, mean_mean_width_ft:temp)
write_csv(fish_stuff2, 'data/fish_data_mean.csv')
