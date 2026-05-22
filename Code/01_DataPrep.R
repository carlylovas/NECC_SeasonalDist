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

sf_use_s2(FALSE)

# Initial data cut off for tows
cut<- 0.70 # Percent of tows to keep a species, 0.75 add cusk and lanternfish that are a problem

# Load and preliminary cleaning of raw data ----
survdat <- readRDS(here::here("Data/survdat_lw.rds"))$survdat |>
    as.data.frame()

# Some clean up
trawldat <- janitor::clean_names(survdat)

# Keep tows only within EPUs
# epus <- bind_rows(
#     st_read(here::here("Data/individual_epus/GOM.geojson"), quiet = TRUE),
#     st_read(here::here("Data/individual_epus/GB.geojson"), quiet = TRUE),
#     st_read(here::here("Data/individual_epus/MAB.geojson"), quiet = TRUE)
# )

# trawldat <- st_as_sf(trawldat, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
# trawldat <- st_join(trawldat, epus, join = st_within)
# trawldat <- st_drop_geometry(trawldat)
# trawldat <- dplyr::filter(trawldat, !is.na(EPU))    


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

# Remove problem years
trawldat<- trawldat |>
    filter(!year %in% c(2017, 2020))

# Getting distinct biomass values at the species level
dat_clean <- trawldat |>
    distinct(id, svspp, catchsex, comname, year, est_month, est_day, season, lat, lon, est_towdate, biomass_kg) |>
    group_by(id, svspp, comname, year, est_month, est_day, season, lat, lon, est_towdate) |>
    summarize("total_biomass_kg" = sum(biomass_kg)) |>
    ungroup()

# Species filtering ----
# Keep only species observed in >= 5 tows per season, in both seasons, for >= cut years
cut <- (max(dat_clean$year) - min(dat_clean$year)) - floor((1 - cut) * (max(dat_clean$year) - min(dat_clean$year)))

tow_seas_spp <- dat_clean |>
    filter(total_biomass_kg > 0) |>
    group_by(svspp, comname, year, season) |>
    summarise(tows = n_distinct(id)) |>
    filter(tows >= 5) |>
    group_by(svspp, comname, year) |>
    summarise(seasons = n_distinct(season)) |>
    filter(seasons == 2) |>
    group_by(svspp, comname) |>
    summarise(years = n_distinct(year)) |>
    filter(years >= cut)

# Summaries and saving prepped data ----
dat_out <- dat_clean |>
    filter(comname %in% tow_seas_spp$comname) # This was keeping cusk and lanternfish with problematic years
unique(dat_out$comname)

# Filter out crabs/lobster
spp_out<- c("atlantic rock crab", "jonah crab", "american lobster", "sea scallop")
dat_out<- dat_out |>
    filter(!comname %in% spp_out)
unique(dat_out$comname)
# dat_out <- dat_out |>
#     semi_join(tow_spp, by = c("svspp", "comname", "year", "season")) |>
#     filter(comname %in% tow_seas_spp$comname)

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
summary(year_day)

delta_year_day <- year_day |>
    group_by(year, season) |>
    summarise(median = median(year_day)) |>
    pivot_wider(names_from = "season", values_from = "median") |>
    mutate(delta_year_day = (Fall - Spring))

saveRDS(delta_year_day, here::here("Data/trawl_yearday_diff.rds"))

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

epu_label_df<- epus |>
    st_centroid() |>
    st_coordinates() |>
    as.data.frame() |>
    bind_cols(epus) |>
    select(EPU, X, Y)

# --- Extract shared boundary lines between EPU pairs ---
gom <- epus |> filter(EPU == "GOM")
gb  <- epus |> filter(EPU == "GB")
mab <- epus |> filter(EPU == "MAB")

# Get the full boundary of each EPU as a linestring
boundary_mab <- st_boundary(mab) |> st_cast("LINESTRING")
boundary_gb  <- st_boundary(gb)  |> st_cast("LINESTRING")

# Crop MAB boundary:
bbox_mab_crop <- st_bbox(c(
  xmin = -71.5, xmax = -69.0,
  ymin = 40.05,  ymax = 41.25
), crs = 4326) |> st_as_sfc()

line_mab_crop <- st_intersection(boundary_mab, bbox_mab_crop)

# Quick check
plot(st_geometry(mab), main = "MAB crop check")
plot(line_mab_crop, col = "red", add = TRUE, lwd = 2)

# --- Helper: connect endpoints with a straight line ---
endpoints_to_line <- function(line_sf, connect = c("north_south", "east_west")) {
  connect <- match.arg(connect)
  
  coords <- line_sf |>
    st_collection_extract("LINESTRING") |>
    st_coordinates()
  coords <- coords[, 1:2]
  
  if (connect == "north_south") {
    pt1 <- coords[which.max(coords[, 2]), ]  # furthest north
    pt2 <- coords[which.min(coords[, 2]), ]  # furthest south
  } else {
    pt1 <- coords[which.min(coords[, 1]), ]  # furthest west
    pt2 <- coords[which.max(coords[, 1]), ]  # furthest east
  }
  
  ends <- rbind(pt1, pt2)
  st_sfc(st_linestring(ends), crs = 4326)
}

line_mab_straight <- endpoints_to_line(line_mab_crop, connect = "north_south")

# Quick check
plot(st_geometry(mab), main = "MAB straight line check")
plot(line_mab_straight, col = "red", add = TRUE, lwd = 2)

# Crop GB boundary
# running roughly parallel from ~71W to ~66W
bbox_gb_crop <- st_bbox(c(
  xmin = -70, xmax = -66,
  ymin = 41.25,  ymax = 43
), crs = 4326) |> st_as_sfc()

line_gb_crop <- st_intersection(boundary_gb, bbox_gb_crop)

line_gb_straight  <- endpoints_to_line(line_gb_crop,  connect = "east_west")

line_gb_straight <- st_sfc(
  st_linestring(matrix(c(
    -70.0, 41.67,   # western anchor, keep as is
    -66.0, 42.4    # eastern end, lifted to 42.25
  ), ncol = 2, byrow = TRUE)),
  crs = 4326
)

# Quick check
plot(st_geometry(gb), main = "GB straight line check")
plot(line_gb_straight, col = "red", add = TRUE, lwd = 2)

# --- Combine for plotting ---
dividers <- bind_rows(
  st_sf(region = "MAB boundary", geometry = st_geometry(line_mab_straight)),
  st_sf(region = "GB boundary",  geometry = st_geometry(line_gb_straight))
)

# --- Sample map ---
dat_sf <- dat |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326)

# --- Plot ---
map <- ggplot() +
    geom_sf(data = dat_sf, aes(color = season), pch = 21, size = 0.5, alpha = 0.25) +
    geom_sf(data = dividers, color = "black", linewidth = 1) +
    geom_sf(data = northeast_states, fill = "grey90", color = "black", size = 0.3) +
    geom_sf(data = atlantic_provinces, fill = "grey85", color = "black", size = 0.3) +
    geom_text(data = epu_label_df, aes(x = X+0.1, y = Y, label = EPU), size = 6, fontface = "bold") +
    coord_sf(xlim = c(-76, -65), ylim = c(35, 44.575)) +
    scale_color_manual(name = "Season", values = c("#d95f02", "#1b9e77")) +
    xlab("Latitude") +  
    ylab("Longitude") +
    theme_bw(base_size = 18) +
    theme(legend.position = "none") +
    facet_wrap(~season, ncol = 2)

ts <- ggplot(delta_year_day, aes(x = year)) +
    geom_line(aes(y = Spring, color = "Spring"), size = 1) +
    geom_point(aes(y = Spring, color = "Spring"), size = 2) +
    geom_line(aes(y = Fall, color = "Fall"), size = 1) +
    geom_point(aes(y = Fall, color = "Fall"), size = 2) +
    scale_color_manual(name = "Season", values = c("Spring" = "#1b9e77", "Fall" = "#d95f02")) +
    labs(y = "Median\nTrawl Day of Year") +
    theme_bw(base_size = 18)

bottemp_ts <- trawldat |>
    distinct(id, year, season, bottemp) |>
    filter(!is.na(bottemp)) |>
    group_by(year, season) |>
    summarize(median_bottemp = median(bottemp, na.rm = TRUE), .groups = "drop")

ts_temp <- ggplot(bottemp_ts, aes(x = year, y = median_bottemp, color = season)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(name = "Season", values = c("Spring" = "#1b9e77", "Fall" = "#d95f02")) +
    labs(y = "Median\nBottom Temp. (°C)", x = "") +
    theme_bw(base_size = 18)

samp_map_out <- map / ts / ts_temp + plot_layout(ncol = 1, heights = c(2, 1, 1))
samp_map_out
ggsave(here::here("Figures/SurveyMap_TrawlDayTimeseries.png"), samp_map_out, width = 11, height = 10, units = "in", dpi = 300)


# Alternative mapping
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)

# --- Extract and crop EPU boundary segments ---

# Get the full boundary of each EPU as a linestring
boundary_mab <- st_boundary(mab) |> st_cast("LINESTRING")
boundary_gb  <- st_boundary(gb)  |> st_cast("LINESTRING")

# Crop MAB boundary: keep only the segment starting at ~71W running 
# south and east to ~69W (the northern/northeastern edge)
bbox_mab_crop <- st_bbox(c(
  xmin = -71.5, xmax = -69.0,
  ymin = 40.25,  ymax = 41.5
), crs = 4326) |> st_as_sfc()

line_mab_crop <- st_intersection(boundary_mab, bbox_mab_crop)

# Quick check
plot(st_geometry(mab), main = "MAB crop check")
plot(line_mab_crop, col = "red", add = TRUE, lwd = 2)

# Crop GB boundary: keep only the southern high-latitude edge
# running roughly parallel from ~71W to ~66W
bbox_gb_crop <- st_bbox(c(
  xmin = -71.5, xmax = -66.75,
  ymin = 41.25,  ymax = 43
), crs = 4326) |> st_as_sfc()

line_gb_crop <- st_intersection(boundary_gb, bbox_gb_crop)

plot(st_geometry(gb), main = "GB crop check")
plot(line_gb_crop, col = "red", add = TRUE, lwd = 2)

# --- Combine for plotting ---
dividers <- bind_rows(
  st_sf(region = "MAB boundary", geometry = st_geometry(line_mab_crop)),
  st_sf(region = "GB boundary",  geometry = st_geometry(line_gb_crop))
)


# --- EPU label positions (adjust X/Y after inspecting map) ---
epu_labels <- tibble(
  EPU = c("Gulf of Maine", "Georges Bank", "Mid-Atlantic Bight"),
  X   = c(-68.5, -67.0, -73.5),
  Y   = c(43.2,  40.8,  38.5)
)

# --- Map ---
map <- ggplot() +
  geom_sf(data = dat_sf, aes(color = season), pch = 21, size = 0.5, alpha = 0.25) +
  geom_sf(data = northeast_states,   fill = "grey90", color = "black", linewidth = 0.3) +
  geom_sf(data = atlantic_provinces, fill = "grey85", color = "black", linewidth = 0.3) +
  geom_sf(data = dividers, color = "black", linewidth = 1) +
  geom_text(data = epu_labels, aes(x = X, y = Y, label = EPU),
            size = 4, fontface = "bold", color = "grey30") +
  coord_sf(xlim = c(-76, -65), ylim = c(35, 44.575)) +
  scale_color_manual(name = "Season", values = c("#d95f02", "#1b9e77")) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  facet_wrap(~season, ncol = 2)

map
