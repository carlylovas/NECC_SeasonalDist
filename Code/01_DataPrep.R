#####
## Loading and prepping Northeast Fisheries Science Center Bottom Trawl Survey Data
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

# Load and preliminary cleaning of raw data ----
survdat <- readRDS(here::here("Data/survdat_lw.rds"))$survdat |>
    as.data.frame()

# Some clean up
trawldat <- janitor::clean_names(survdat)

# Add in species common name
spp_classes <- readr::read_csv(here::here("Data/sppclass.csv"),
    col_types = readr::cols()
)
spp_classes <- janitor::clean_names(spp_classes)
spp_classes <- dplyr::mutate(
    .data = spp_classes, comname = stringr::str_to_lower(common_name),
    scientific_name = stringr::str_to_lower(scientific_name)
)
spp_classes <- dplyr::distinct(spp_classes, svspp, comname, scientific_name)
trawldat <- dplyr::mutate(trawldat, svspp = stringr::str_pad(svspp, 3, "left", "0"))
trawldat <- dplyr::left_join(trawldat, spp_classes, by = "svspp")

# Creating a unique tow ID column
trawldat <- dplyr::mutate(.data = trawldat, cruise6 = stringr::str_pad(
    cruise6,
    6, "left", "0"
), station = stringr::str_pad(
    station,
    3, "left", "0"
), stratum = stringr::str_pad(
    stratum,
    4, "left", "0"
), id = stringr::str_c(
    cruise6, station,
    stratum
))

# Adding a date column
trawldat <- dplyr::mutate(.data = trawldat, est_month = stringr::str_sub(
    est_towdate,
    6, 7
), est_month = as.numeric(est_month), est_day = stringr::str_sub(
    est_towdate,
    -2, -1
), est_day = as.numeric(est_day), .before = season)

# Column names/formatting
trawldat <- dplyr::mutate(.data = trawldat, comname = tolower(comname), id = format(id, scientific = FALSE), svspp = as.character(svspp), svspp = stringr::str_pad(svspp, 3, "left", "0"), season = stringr::str_to_title(season), strat_num = stringr::str_sub(stratum, 2, 3))
trawldat <- dplyr::rename(.data = trawldat, biomass_kg = biomass, length_cm = length)

# Dealing with when there is biomass/no abundance, or abundance but no biomass
trawldat <- dplyr::mutate(.data = trawldat, biomass_kg = ifelse(biomass_kg == 0 & abundance > 0, 1e-04, biomass_kg), abundance = ifelse(abundance == 0 & biomass_kg > 0, 1, abundance))
trawldat <- dplyr::filter(.data = trawldat, !is.na(biomass_kg), !is.na(abundance))

# Filtering strata not regularly sampled throughout the time series
trawldat <- dplyr::filter(.data = trawldat, stratum >= 1010, stratum <= 1760, stratum != 1310, stratum != 1320, stratum != 1330, stratum != 1350, stratum != 1410, stratum != 1420, stratum != 1490)

# Filtering species not regularly sampled (shrimps, others?)
trawldat <- dplyr::filter(.data = trawldat, !svspp %in% c(285:299, 305, 306, 307, 316, 323, 910:915, 955:961))
trawldat <- dplyr::filter(trawldat, !svspp %in% c(0, "000", 978, 979, 980, 998))

trawldat <- dplyr::filter(trawldat, year >= 1970)

# Getting distinct biomass values at the species level
dat_clean <- trawldat |>
    distinct(id, svspp, catchsex, comname, year, est_month, est_day, season, lat, lon, est_towdate, biomass_kg) |>
    group_by(id, svspp, comname, year, est_month, est_day, season, lat, lon, est_towdate) |>
    summarize("total_biomass_kg" = sum(biomass_kg)) |>
    ungroup()

# Species filtering ----
# Keep only species that were observed in at least 5 tows for each season and then in both seasons for at least 80% of survey years.
tow_spp <- dat_clean |>
    group_by(svspp, comname, year, season) |>
    summarise(tows = n_distinct(id)) |>
    filter(tows >= 5)

# 80% cut off (49 years)
cut <- (max(tow_spp$year) - min(tow_spp$year)) - floor(0.08 * (max(tow_spp$year) - min(tow_spp$year)))

tow_seas_spp <- tow_spp |>
    # 80% of years have both spring and fall
    group_by(svspp, comname, year) |>
    summarise(seasons = n_distinct(season)) |>
    filter(seasons == 2) |>
    group_by(svspp, comname) |>
    summarise(years = n_distinct(year)) |>
    filter(years >= cut)

# Summaries and saving prepped data ----
dat_out <- dat_clean |>
    filter(comname %in% tow_seas_spp$comname)

saveRDS(dat_out, here::here("Data/dat_clean.rds"))

# Median trawl day ----
year_day <- dat_out |>
    mutate(year_day = lubridate::yday(est_towdate)) |>
    select(year, season, year_day, id)

delta_year_day<- year_day |> 
  group_by(year, season) |>
  summarise(median = median(year_day)) |>
  pivot_wider(names_from = "season", values_from = "median") |> 
  mutate(delta_year_day = (Fall-Spring))

saveRDS(delta_year_day, here::here("Data/trawl_yearday_diff.rds"))

rm(list=ls())
