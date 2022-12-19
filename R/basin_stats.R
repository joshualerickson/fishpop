#### getting the zonal stats per year eda

library(sf)
library(terra)
library(rgee)
library(tidyrgee)
library(exploreRGEE)
ee_Initialize()

library(tidyverse)
# read in spatial data we want to look at

fish_spat <- read_sf('data/eda_fish_spatial.gpkg', 'eda_fish_spatial')%>% 
  filter(strm_sg != 'Fortine 2')
fish_all <- read_csv('data/fish_data_all.csv')

fish_spat <- fish_spat %>% mapedit::editFeatures()
# get basins 
basins <- fish_spat %>% 
          split(.$strm_sg) %>% 
          map(~nhdplusTools::get_split_catchment(.))

basins2 <- fish_spat %>% 
  filter(strm_sg %in% c('Deep-1', 'Dodge-2',
                        'Fortine-2',
                        'Pinkham-4', 'Pinkham-5')) %>% 
  split(.$strm_sg) %>% 
  map(~nhdplusTools::get_split_catchment(.))

names_basins <- names(basins)

basins_merge <- map2(basins, names_basins, ~.x %>% mutate(stream_segment = .y))




basins_merge <- bind_rows(basins_merge) %>% 
                filter(id %in% c('mergedCatchment', 'splitCatchment'))
basins_merge %>%  
  mapview::mapview(alpha.regions = 0)  + 
  mapview::mapview(fish_spat,
                    cex = 1,
                    color = 'red') 

ic <- ee$ImageCollection("GRIDMET/DROUGHT")

ic_tidy <- as_tidyee(ic) 

ic_tidy <- ic_tidy %>% 
           filter(month == 8)

ic_tidy <- ic_tidy %>% 
          filter(doy == 235)

ee_Initialize(gcs = T)  
basins_merge <- basins_merge %>% st_cast('MULTIPOLYGON')

bas_fc <- sf_as_ee(basins_merge)

ext_basins <- ee_extract_tidy(ic_tidy, y = bas_fc, 
                scale  = 4000,
                via = 'gcs',
                container = 'jle_rasters')  
# 
# write_csv(ext_basins, 'data/ext_basins.csv')
# ext_basins <- read_csv('data/ext_basins.csv')
  
ext_basins <- ext_basins %>% mutate(year = lubridate::year(date))


ext_basins <- ext_basins %>% left_join(fish_all, by = c('stream_segment', 'year'))

mod_df <- ext_basins %>%   
  filter(str_detect(parameter, 'spi2y|spi1y|spi5y|spei2y|spei1y|spei5y|eddi2y|eddi1y|eddi5y')) %>% 
  pivot_wider(names_from = parameter, values_from = value) %>% 
  filter(!is.na(density_fish_m2)) %>% 
  mutate(year = abs(min(year)-year))

mod_pca <- ext_basins %>%   
  filter(str_detect(parameter, 'spi2y|spi1y|spi5y|spei2y|spei1y|spei5y|eddi2y|eddi1y|eddi5y')) %>% 
  pivot_wider(names_from = parameter, values_from = value) %>% 
  filter(!is.na(density_fish_m2)) %>% 
  mutate(year = abs(min(year)-year)) %>% 
  select(eddi1y:spi5y, stream_segment) %>% 
  group_by(stream_segment) %>% 
  nest() %>% 
  mutate(pca1 = map(data, ~pcaMethods::pca(., nPcs = 1) %>% 
                      .@sDev)) %>% 
  select(stream_segment, pca1) %>% 
  ungroup() %>% 
  mutate(pca1 = as.numeric(pca1))

mod_pca <- mod_df %>% left_join(mod_pca, by = 'stream_segment')

library(lme4)
lmm_year <- lmer(density_fish_m2 ~ area_acres +
                   pca1 + 
                   year + 
                   (year|stream_segment), 
                 mod_pca,
                 control = lmerControl(optimizer = 'Nelder_Mead'))
jtools::summ(lmm_year)

library(lme4)
lmm_year <- lmer(density_fish_m2 ~ area_acres + 
                   silt + 
                      year + 
                      (year|stream_segment) , 
                    mod_df ,
                    control = lmerControl(optimizer = 'Nelder_Mead'))
jtools::summ(lmm_year)

lmm_precip <- lmer(density_fish_m2 ~ area_acres + 
                     silt + 
                     spei5y +
                   (spei5y|stream_segment), 
                 mod_df ,
                 control = lmerControl(optimizer = 'Nelder_Mead'))

anova(lmm_year, lmm_precip)
jtools::summ(lmm_precip)
summary(lmm_precip)
m <- lm(density_fish_m2~area_acres+silt+spei5y, data = mod_df)
summary(m)
mod_mean <- read_csv('data/fish_data_mean.csv')
modm1 <- lm(mean_density_fish_m2~area_acres+pca1,
            data = mod_mean %>% left_join(mod_pca %>% 
                                            select(pca1, stream_segment), by = 'stream_segment'))
summary(modm1)
mod_mean %>% left_join(mod_pca %>% 
                         select(pca1, stream_segment), by = 'stream_segment') %>% 
  mutate(area_cut = cut_interval(area_acres, 3)) %>% 
  ggplot(aes(area_acres, mean_density_fish_m2)) + 
  geom_point(aes(color = pca1), size = 3) +
  geom_smooth(method = 'lm', se = F) +
  stat_smooth(geom = 'line',
              aes(group = area_cut), 
              method = 'lm') + 
  resourceviz::custom_theme()

anova(lmm_year, lmm_precip)
jtools::summ(lmm_precip)
summary(lmm_precip)
ext_basins %>% 
  filter(str_detect(parameter, 'spi2y|spi1y|spi5y|spei2y|spei1y|spei5y|eddi2y|eddi1y|eddi5y')) %>% 
  pivot_wider(names_from = parameter, values_from = value) %>% 
  select(spi2y, spi1y, spi5y,spei2y,spei1y,spei5y,eddi2y,eddi1y,eddi5y, year,stream_segment) %>% view()
  group_by(stream_segment) %>% nest() %>% view() 
  mutate(pca1 = map(data, ~pcaMethods::pca(., nPcs = 1) %>% 
                      .@completeObs)) %>% 
  select(-data) %>% 
  view()
  ggplot(aes(value, density_fish_m2, color = parameter)) + 
  geom_point() +
  geom_smooth(aes(group = parameter), method = 'lm', se = F) + 
  facet_wrap(~stream_segment) + 
  resourceviz::custom_theme() 


mod_df %>% 
  #filter(str_detect(parameter, 'eddi2y|eddi1y|eddi5y|spei5y')) %>% 
  ggplot(aes(spei5y, density_fish_m2)) + 
  geom_point() +
  geom_smooth(method = 'lm', se = F) + 
  facet_wrap(~stream_segment, scales = 'free') + 
  resourceviz::custom_theme() 

ext_basins %>% 
  filter(str_detect(parameter, 'spei2y|spei5y')) %>% 
  ggplot(aes(value, density_fish_m2, color = parameter)) + 
  geom_point() +
  geom_smooth(aes(group = parameter), method = 'lm', se = F) + 
  facet_wrap(~stream_segment) + 
  resourceviz::custom_theme() 

ext_basins %>% 
  filter(str_detect(parameter, 'spei2y|spei5y')) %>% 
  ggplot(aes(value, density_fish_m2, color = parameter)) + 
  geom_point() +
  geom_smooth(aes(group = parameter), method = 'lm', se = F) + 
  facet_wrap(~stream_segment) + 
  resourceviz::custom_theme() 

ext_basins %>% 
  filter(str_detect(parameter, 'spei2y')) %>% 
  ggplot(aes(value, density_fish_m2)) + 
  geom_point(aes(label = year)) 

fish_all %>%  
  ggplot(aes(pools_mile, density_fish_m2)) + 
  geom_point() +
  geom_smooth(method = 'lm', se = F) + 
  facet_wrap(~stream_segment) + 
  resourceviz::custom_theme() 

resourceviz::cairo_view()

lmm_precip <- lmer(density_fish_m2 ~ area_acres + 
                     silt + 
                     spi5y +
                     (spei5y|stream_segment), 
                   mod_df ,
                   control = lmerControl(optimizer = 'Nelder_Mead'))
summary(lmm_precip)
jtools::summ(lmm_precip)
