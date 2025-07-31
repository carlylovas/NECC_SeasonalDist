####
## Scratch
#### 

# Libraries ----
library(here)
library(tidyverse)
library(gmRi)
library(matrixStats)
library(patchwork)
library(broom)
library(glmmTMB)
library(marginaleffects)
library(ggeffects)

# Load cleaned data ----
dat <- readRDS(here::here("Data/dat_clean.rds"))

delta_year_day <- readRDS(here::here("Data/trawl_yearday_diff.rds")) |>
  select(year, delta_year_day)

# Calculate center of biomass metrics by season ----
center_bio <- function(dat, ...){
  dat %>%
    group_by(comname, ...) %>%
    summarise(
      # Un-weighted averages
      total_biomass   = sum(total_biomass_kg),
      avg_biomass     = mean(total_biomass_kg),
      biomass_sd      = sd(total_biomass_kg),
      avg_lat         = weightedMean(lat, w = total_biomass_kg, na.rm = T),  
      avg_lon         = weightedMean(lon, w = total_biomass_kg, na.rm = T),
      .groups = "drop")
}

seas_cob <- center_bio(dat, year, season)

# Fit linear mixed model to get trend in latitude and longitude over time ----
seas_cob <- seas_cob |> 
  left_join(delta_year_day) |> 
  select(comname, year, season, avg_lat, delta_year_day)

mod <- glmmTMB(avg_lat ~ year*comname*season, data = seas_cob, REML = T)

slopes_output <- slopes(mod, variables = "year", by = c("comname", "season"), type = "response")

# define movement trends and categories
seas_cob_coeffs |>
  filter(variable == "lat") |>
  mutate(signif   = ifelse(p.value <= 0.05, "true", "false"),
         movement = ifelse(estimate > 0, "+", "-")) |>
  select(!p.value) -> seas_trends

seas_trends |> 
  filter(signif == "true" & movement == "+") |>
  mutate(trend = "northward") %>%
  full_join(seas_trends|>
              filter(signif == "true" & movement == "-") |>
              mutate(trend = "southward")) |>
  full_join(seas_trends|>
              filter(signif == "false" & movement == "-" | signif == "false" & movement == "+") |>
              mutate(trend = "stable")) %>%
  arrange(comname, season) -> seas_trends

seas_trends |> 
  select(comname, season, estimate)|>
  pivot_wider(names_from = season, values_from = estimate) |>
  group_by(comname)|> 
  summarise(difference = `Fall`-`Spring`) |> 
  left_join(seas_trends |> ungroup() |> select(comname, season, signif)) |>
  pivot_wider(names_from = season, values_from = signif) -> seas_diff

## contracting
seas_diff |>
  filter(Fall == "false" & Spring == "true" & difference < 0) |> 
  select(comname)|>
  mutate(category = "contracting") -> contracting 

## marching
### fall > 0 < spring & 0 > fall-spring < 1
seas_diff|>
  filter(Fall == "true" & Spring == "true" & abs(difference) < 1) |> # figure out this threshold
  select(comname)|>
  mutate(category = "marching")-> marching 

## expanding
### fall-spring > 0 & fall > 0 > spring
seas_trends|> # fix this
  ungroup() |>
  select(comname, season, estimate) |>
  pivot_wider(names_from = season, values_from = estimate) |>
  filter(Fall > 0 & Spring < 0) |>
  select(comname) |>
  left_join(seas_diff) |>
  filter(!Fall == "false" |!Spring == "false") |>
  select(comname) |>
  mutate(category = "expanding") -> expanding

## stable 
seas_diff |> 
  filter(Fall == "false" & Spring == "false") |>
  select(comname) |>
  mutate(category = "stable") -> stable

trends <- contracting %>%
  full_join(marching) %>%
  full_join(expanding) %>% 
  full_join(stable) 

write_rds(trends, here("Data", "movement_categories.rds"))
### haddock, northern sea robin, northern shortfin squid, and smooth skate fall out...


# Calculating distance between spring and fall centroids ----
seasonal_dist <- seas_cob %>% 
  unnest(data) %>% 
  select(comname, year, season, avg_lat, avg_lon) %>%
  group_by(comname, year ) %>% 
  nest() 

point_dist <- function(df){
  temp <- sf::st_as_sf(df,coords=c("avg_lon","avg_lat"), crs=4326, remove=FALSE)
  temp <- sf::st_transform(temp, crs = 32619)
  out  <- sf::st_distance(temp$geometry)[1,2]
  return(out)
}

seasonal_dist <- seasonal_dist |>
  mutate(dist_m = map_dbl(data, possibly(point_dist, NA)), 
         dist_km = dist_m/1000) |>
  #   filter(!is.na(dist_m)) |> 
  group_by(comname) |> 
  nest() |> 
  mutate(count = map_dbl(data, function(x){nrow(x)}), 
         season_dist = map(data, function(x){
           lm(dist_km ~ year, data = x) |> 
             tidy() |> 
             filter(term == "year")})) |> 
  unnest(season_dist) %>% 
  mutate(signif = ifelse(p.value <= 0.05, T, F)) 

# Save out the seasonal distance difference data for subsequent modeling
seasonal_dist_out <- seasonal_dist |>
  ungroup() |>
  select(comname, data) |>
  unnest(cols = c(data)) |>
  unnest(cols = c(data))
write_rds(seasonal_dist_out, here("Data","seasonal_dist.rds"))

rm(list=ls())