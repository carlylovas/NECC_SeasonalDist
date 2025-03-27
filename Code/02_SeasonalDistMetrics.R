#####
## Calculating seasonal distance metrics and generating species profiles
#####

# Libraries ----
library(here)
library(tidyverse)
library(gmRi)
library(matrixStats)
library(patchwork)
library(broom)
library(rmarkdown)
library(shiny)
library(marginaleffects)

# Load cleaned data ----
dat <- readRDS(here::here("Data/dat_clean.rds")) |>
  filter(!year %in% c(2017,2020))

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

# Find distance along coast ---- 
coastdistdat <- readRDS(here("Data","coastdistdat.rds"))

get_length <- function(df, distdf = coastdistdat){
  tmp <- distdf %>% 
    mutate(abs.diff.x2 = abs(x-df$avg_lon)^2,
           abs.diff.y2 = abs(y-df$avg_lat)^2,
           abs.diff.xy = sqrt(abs.diff.x2 + abs.diff.y2
           ))  %>% 
    filter(abs.diff.xy == min(abs.diff.xy)) %>%
    dplyr::select(lengthfromhere) %>%
    pull()
  return(tmp)
}

seas_cob |>
  group_by(comname, season, year) |>
  nest() |>
  mutate(coastdist = map_dbl(data, get_length),
         coastdist_km = coastdist/1000) |>
  select(comname, season, year, coastdist_km) -> coastdist_km

# save out for plotting
write_rds(coastdist_km, here("Data", "seas_cob.rds"))

# Linear model with interaction ----
coastdist_km |>
  group_by(comname) |> 
  nest() |>
  mutate(trends = map(data, function(x){
    mod <- lm(coastdist_km~year*season, data = x)
    slopes <- slopes(mod, variables = "year", by = "season", type = "response") |>
    select(!term) |>
    rename("term" = "season")
    return(slopes)
  })) |>
  select(!data) |>
  unnest(trends) -> seasonal_trends # name tbd

coastdist_km |>
  group_by(comname) |> 
  nest() |>
  mutate(outputs = map(data, function(x){
    mod <- lm(coastdist_km~year*season, data = x)
    slopes <- slopes(mod, variables = "year", by = "season", type = "response", hypothesis = "pairwise")
    return(slopes)
  })) |>
  select(!data) |>
  unnest(outputs) -> difference_in_season

seasonal_trends |>
  full_join(difference_in_season) |> 
  arrange(comname) -> all_outputs

# Define movement trends and categories ----
all_outputs |>
  select(comname, term, estimate, p.value) |>
  mutate(signif   = ifelse(p.value <= 0.05, "true", "false"),
         movement = ifelse(estimate > 0, "+", "-")) -> output_trends

output_trends |> 
  mutate(trend = case_when(
    term %in% c("Fall", "Spring") & signif == "true" & movement == "+" ~ "northward",
    term %in% c("Fall", "Spring") & signif == "true" & movement == "-" ~ "southward",
    term == "Fall - Spring" & signif == "true" & movement == "+" ~ "more different",
    term == "Fall - Spring" & signif == "true" & movement == "-" ~ "less different",
    signif == "false" ~ "stable")) -> all_output_trends

# Movement categories ----
all_output_trends |>
  select(comname, term, trend) |>
  pivot_wider(names_from = "term", values_from = "trend") |>
  mutate(category = case_when(
    # contracting
    `Fall` %in% c("stable", "southward") & `Spring` == "northward" ~ "contracting",
    `Fall` == "southward" & `Spring` == "stable" ~ "contracting",
    # marching
    `Fall` == "northward" & `Spring` == "northward" & `Fall - Spring` %in% c("stable", "less different") ~ "marching",
    `Fall` == "southward" & `Spring` == "southward" & `Fall - Spring` %in% c("stable", "less different") ~ "marching",
    # expanding 
    `Fall` == "northward" & `Spring` %in% c("southward", "stable") & `Fall - Spring` %in% c("stable", "more different") ~ "expanding",
    `Fall` == "stable" & `Spring` == "southward" & `Fall - Spring` %in% c("stable", "more different") ~ "expanding",
    # stable
    `Fall` == "stable" & `Spring` == "stable" ~ "stable"
  )) -> categories

write_rds(categories, here("Data", "movement_categories.rds"))

# Calculating distance between spring and fall centroids ----
seasonal_dist <- seas_cob %>% 
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
  group_by(comname) |> 
  nest() |> 
  mutate(count = map_dbl(data, function(x){nrow(x)}), 
         season_dist = map(data, function(x){
           lm(dist_km ~ year, data = x) |> 
                          tidy() |> 
                          filter(term == "year")})) |> 
  unnest(season_dist) %>% 
  mutate(signif = ifelse(p.value <= 0.05, T, F)) 

# Save out the seasonal distance difference data for subsequent modeling ----
seasonal_dist_out <- seasonal_dist |>
    ungroup() |>
    select(comname, data) |>
    unnest(cols = c(data)) |>
    unnest(cols = c(data))
write_rds(seasonal_dist_out, here("Data","seasonal_dist.rds"))

rm(list=ls())

