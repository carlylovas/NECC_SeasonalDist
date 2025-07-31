#####
## Calculating seasonal distance metrics
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

# Load cleaned data ----
dat <- readRDS(here::here("Data/dat_clean.rds"))

# Calculate center of biomass metrics by season ----
center_bio <- function(dat, ...) {
  dat %>%
    group_by(comname, ...) %>%
    summarise(
      # Un-weighted averages
      total_biomass = sum(total_biomass_kg),
      avg_biomass = mean(total_biomass_kg),
      biomass_sd = sd(total_biomass_kg),
      avg_lat = weightedMean(lat, w = total_biomass_kg, na.rm = T),
      avg_lon = weightedMean(lon, w = total_biomass_kg, na.rm = T),
      .groups = "drop"
    )
}

seas_cob <- center_bio(dat, year, season)

# Calculating distance between spring and fall centroids ----
seasonal_dist <- seas_cob %>%
  select(comname, year, season, avg_lat, avg_lon) %>%
  group_by(comname, year) %>%
  nest()

point_dist <- function(df) {
  temp <- sf::st_as_sf(df, coords = c("avg_lon", "avg_lat"), crs = 4326, remove = FALSE)
  temp <- sf::st_transform(temp, crs = 32619)
  out <- sf::st_distance(temp$geometry)[1, 2]
  return(out)
}

seasonal_dist <- seasonal_dist |>
  mutate(
    dist_m = map_dbl(data, possibly(point_dist, NA)),
    dist_km = dist_m / 1000
  ) |>
  #   filter(!is.na(dist_m)) |>
  group_by(comname) |>
  nest() |>
  mutate(
    count = map_dbl(data, function(x) {
      nrow(x)
    }),
    season_dist = map(data, function(x) {
      lm(dist_km ~ year, data = x) |>
        tidy() |>
        filter(term == "year")
    })
  ) |>
  unnest(season_dist) %>%
  mutate(signif = ifelse(p.value <= 0.05, T, F))

# Save out the seasonal distance difference data for subsequent modeling
seasonal_dist_out <- seasonal_dist |>
  ungroup() |>
  select(comname, data) |>
  unnest(cols = c(data)) |>
  unnest(cols = c(data))
write_rds(seasonal_dist_out, here("Data", "seasonal_dist.rds"))
