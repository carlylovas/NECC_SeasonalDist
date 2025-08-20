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
library(rnaturalearth)
library(rnaturalearthhires)
library(sf)

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

# LOBSTER PROBLEM
spp_classes$scientific_name[which(spp_classes$svspp == 301)] <- "Homarus americanus"

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

# Now cut 2017 and 2020
dat_out <- dat_clean |>
    filter(!year %in% c(2017, 2020))


# Summaries and saving prepped data ----
dat_out <- dat_out |>
    filter(comname %in% tow_seas_spp$comname)

# Quick time series by season/species of total catch
ts_dat <- dat_out |>
    group_by(comname, year, season) |>
    summarize("total" = sum(total_biomass_kg))

ggplot() +
    geom_point(data = ts_dat, aes(x = year, y = total, color = season)) +
    geom_line(data = ts_dat, aes(x = year, y = total, color = season)) +
    facet_wrap(~comname, scales = "free_y") +
    theme_bw()

saveRDS(dat_out, here::here("Data/dat_clean.rds"))

# Median trawl day ----
year_day <- dat_out |>
    mutate(year_day = lubridate::yday(est_towdate)) |>
    select(year, season, year_day, id)

delta_year_day <- year_day |>
    group_by(year, season) |>
    summarise(median = median(year_day)) |>
    pivot_wider(names_from = "season", values_from = "median") |>
    mutate(delta_year_day = (Fall - Spring))

saveRDS(delta_year_day, here::here("Data/trawl_yearday_diff.rds"))

rm(list = ls())

# Quick sample map of the data and timeseries of median trawl day by season
dat <- readRDS(here::here("Data/dat_clean.rds"))

delta_year_day <- readRDS(here::here("Data/trawl_yearday_diff.rds"))

# --- US states ---
us_states <- ne_states(country = "united states of america", returnclass = "sf")

northeast_states <- us_states %>%
    filter(name %in% c(
        "Maine", "New Hampshire", "Vermont", "Massachusetts",
        "Rhode Island", "Connecticut", "New York", "New Jersey",
        "Pennsylvania", "North Carolina", "Virginia", "West Virginia",
        "Maryland", "Delaware", "District of Columbia"
    ))

# --- Canada provinces ---
canada_provinces <- ne_states(country = "canada", returnclass = "sf")

atlantic_provinces <- canada_provinces %>%
    filter(name %in% c(
        "New Brunswick", "Nova Scotia", "Prince Edward Island",
        "Newfoundland and Labrador", "Quebec", "Ontario"
    ))

# --- Example sample points (replace with your own data) ---
dat_sf <- dat |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326)

# --- Plot ---
map <- ggplot() +
    geom_sf(data = dat_sf, aes(color = season), pch = 21, size = 0.5, alpha = 0.25) +
    geom_sf(data = northeast_states, fill = "grey90", color = "black", size = 0.3) +
    geom_sf(data = atlantic_provinces, fill = "grey85", color = "black", size = 0.3) +
    coord_sf(xlim = c(-76, -65), ylim = c(35, 44.575)) +
    scale_color_manual(name = "Season", values = c("#d95f02", "#1b9e77")) +
    theme_bw(base_size = 16) +
    facet_wrap(~season, ncol = 2)

ts <- ggplot(delta_year_day, aes(x = year)) +
    geom_line(aes(y = Spring, color = "Spring"), size = 1) +
    geom_point(aes(y = Spring, color = "Spring"), size = 2) +
    geom_line(aes(y = Fall, color = "Fall"), size = 1) +
    geom_point(aes(y = Fall, color = "Fall"), size = 2) +
    scale_color_manual(name = "Season", values = c("Spring" = "#1b9e77", "Fall" = "#d95f02")) +
    labs(y = "Median Trawl Day of Year") +
    theme_bw(base_size = 16)

samp_map_out <- map / ts + plot_layout(ncol = 1, heights = c(2, 1), widths = c(2, 1))
ggsave(here::here("Figures/SurveyMap_TrawlDayTimeseries.png"), samp_map_out, width = 11, height = 8, units = "in", dpi = 300)
