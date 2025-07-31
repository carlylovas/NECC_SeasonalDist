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






# Fit linear model to get trend in latitude and longitude over time ----
seas_cob <- seas_cob |>
    group_by(comname, season) |>
    nest() |>
    mutate(
        lat = map(data, function(x) {
            lm(avg_lat ~ year, data = x) |> # 5-year rolling mean?
                tidy() |>
                filter(term == "year")
        }),
        lon = map(data, function(x) {
            lm(avg_lon ~ year, data = x) |>
                tidy() |>
                filter(term == "year")
        })
    )

# isolate coefficients
seas_cob_coeffs <- seas_cob |>
    select(!data) |>
    pivot_longer(cols = lat:lon, names_to = "variable", values_to = "data") |>
    unnest(data) |>
    select(comname, season, variable, estimate, p.value)

# define movement trends and categories
seas_cob_coeffs |>
    filter(variable == "lat") |>
    mutate(
        signif = ifelse(p.value <= 0.05, "true", "false"),
        movement = ifelse(estimate > 0, "+", "-")
    ) |>
    select(!p.value) -> seas_trends

seas_trends |>
    filter(signif == "true" & movement == "+") |>
    mutate(trend = "northward") %>%
    full_join(seas_trends |>
        filter(signif == "true" & movement == "-") |>
        mutate(trend = "southward")) |>
    full_join(seas_trends |>
        filter(signif == "false" & movement == "-" | signif == "false" & movement == "+") |>
        mutate(trend = "stable")) %>%
    arrange(comname, season) -> seas_trends

seas_trends |>
    select(comname, season, estimate) |>
    pivot_wider(names_from = season, values_from = estimate) |>
    group_by(comname) |>
    summarise(difference = `Fall` - `Spring`) |>
    left_join(seas_trends |> ungroup() |> select(comname, season, signif)) |>
    pivot_wider(names_from = season, values_from = signif) -> seas_diff

## contracting
seas_diff |>
    filter(Fall == "false" & Spring == "true" & difference < 0) |>
    select(comname) |>
    mutate(category = "contracting") -> contracting

## marching
### fall > 0 < spring & 0 > fall-spring < 1
seas_diff |>
    filter(Fall == "true" & Spring == "true" & abs(difference) < 1) |> # figure out this threshold
    select(comname) |>
    mutate(category = "marching") -> marching

## expanding
### fall-spring > 0 & fall > 0 > spring
seas_trends |> # fix this
    ungroup() |>
    select(comname, season, estimate) |>
    pivot_wider(names_from = season, values_from = estimate) |>
    filter(Fall > 0 & Spring < 0) |>
    select(comname) |>
    left_join(seas_diff) |>
    filter(!Fall == "false" | !Spring == "false") |>
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


# # Calculating distance between spring and fall centroids ----
# seasonal_dist <- seas_cob %>%
#   unnest(data) %>%
#   select(comname, year, season, avg_lat, avg_lon) %>%
#   group_by(comname, year ) %>%
#   nest()

# point_dist <- function(df){
#   temp <- sf::st_as_sf(df,coords=c("avg_lon","avg_lat"), crs=4326, remove=FALSE)
#   temp <- sf::st_transform(temp, crs = 32619)
#   out  <- sf::st_distance(temp$geometry)[1,2]
#   return(out)
# }

# seasonal_dist <- seasonal_dist |>
#   mutate(dist_m = map_dbl(data, possibly(point_dist, NA)),
#          dist_km = dist_m/1000) |>
# #   filter(!is.na(dist_m)) |>
#   group_by(comname) |>
#   nest() |>
#   mutate(count = map_dbl(data, function(x){nrow(x)}),
#          season_dist = map(data, function(x){
#            lm(dist_km ~ year, data = x) |>
#                           tidy() |>
#                           filter(term == "year")})) |>
#   unnest(season_dist) %>%
#   mutate(signif = ifelse(p.value <= 0.05, T, F))

# # Save out the seasonal distance difference data for subsequent modeling
# seasonal_dist_out <- seasonal_dist |>
#     ungroup() |>
#     select(comname, data) |>
#     unnest(cols = c(data)) |>
#     unnest(cols = c(data))
# write_rds(seasonal_dist_out, here("Data","seasonal_dist.rds"))

# rm(list=ls())
# ------------------------------------------------------------------------- #
# Species summary plots ---- AA stopped here!
# plotting rate of change in seasonal distance

seasonal_dist %>%
    mutate(z = estimate > 0) %>%
    ggplot() +
    geom_point(aes(x = comname, y = estimate, color = as.factor(z))) +
    ggrepel::geom_text_repel(aes(comname, estimate, label = comname), size = 2.8, nudge_y = 0.003) +
    theme_gmri(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none"
    ) +
    scale_color_gmri() +
    ylim(c(-4, 4)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5, color = "black") +
    ylab("Rate of change (km/year)") +
    ggtitle("Changes in seasonal distance")

## Species profile plots
# plot the distance between centroids each year for each species
seasonal_dist %>%
    mutate(plot = map2(data, comname, function(x, y) {
        ggplot(data = x) +
            geom_line(aes(x = est_year, y = dist_km)) +
            theme_gmri() +
            ggtitle(str_to_sentence(comname)) +
            xlab("Year") +
            ylab("Distance between fall and spring centroids (km)")
    })) -> seasonal_dist

seasonal_dist %>%
    select(comname, data) %>%
    mutate(dist_plot = map(data, function(x) {
        ggplot(data = x, aes(x = est_year, y = dist_km)) +
            geom_line(color = "#535353", alpha = 0.8) +
            geom_point(color = "#00608a") +
            ggtitle("Distance between spring and fall centroids") +
            ylab("Distance (km)") +
            xlab("Year") +
            theme_gmri(plot.title = element_text(size = 12))
    })) %>%
    select(!data) -> dist_plot


with_season %>%
    select(comname, data) %>%
    mutate(map = map2(data, comname, function(x, y) {
        us <- rnaturalearth::ne_states(country = "united states of america")
        x <- x %>%
            mutate(decade = 10 * est_year %/% 10)
        plot <- ggplot() +
            geom_sf(data = us) +
            coord_sf(xlim = c(-66, -76), ylim = c(36, 48)) +
            # ggtitle(str_to_sentence(comname)) +
            ggtitle("Season center of biomass by decade") +
            ylab("Latitude") +
            xlab("Longitude") +
            scale_x_continuous(breaks = c(-68, -71, -74)) +
            scale_y_continuous(breaks = c(38, 43, 48)) +
            geom_point(data = x, aes(x = avg_lon, y = avg_lat, color = season)) +
            scale_color_gmri() +
            guides(color = guide_legend(title = "Season")) +
            facet_wrap(~decade, nrow = 1) +
            theme_gmri(
                strip.background = element_rect(fill = "transparent", linetype = 1, linewidth = 1, color = "#e9e9e9"),
                strip.text = element_text(color = "black"),
                # axis.title = element_text(color = "#535353"),
                # axis.text = element_text(color = "#535353"),
                axis.line = element_line(color = "#e9e9e9"),
                axis.ticks = element_line(color = "#e9e9e9"),
                plot.title = element_text(size = 12),
                panel.border = element_rect(color = "#e9e9e9", linetype = 1, linewidth = 1),
                panel.grid.major = element_line(color = "#e9e9e9")
            )

        return(plot)
    })) -> with_season

# bsb_map <- with_season$map[[11]]
# ggsave("bsb_map.png", bsb_map)

with_season %>%
    mutate(time_series = map2(data, comname, function(x, y) {
        x <- x %>%
            mutate(decade = 10 * est_year %/% 10)
        out <- ggplot(data = x) +
            geom_line(aes(x = est_year, y = avg_lat, group = est_year), color = "#535353", alpha = 0.8) +
            geom_point(aes(x = est_year, y = avg_lat, color = season)) +
            scale_color_gmri() +
            theme_gmri(plot.title = element_text(size = 12)) +
            # ggtitle(str_to_sentence(comname)) +
            ggtitle("Seasonal center of latitude") +
            xlab("Year") +
            ylab("Avg. Latitude") +
            guides(color = guide_legend(title = "Season"))
        return(out)
    })) -> with_season

# with_season$time_series[[40]]

grouped_quantiles <- function(clean_survey, ...) {
    clean_survey %>%
        group_by(comname, ...) %>%
        summarise(
            # Un-weighted averages
            total_biomass = sum(biomass_kg),
            avg_biomass = mean(biomass_kg),
            avg_lat = mean(decdeg_beglat),
            # Weight quantiles
            `5%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.05, na.rm = T),
            `10%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.10, na.rm = T),
            `25%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.25, na.rm = T),
            `50%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.50, na.rm = T),
            `75%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.75, na.rm = T),
            `90%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.90, na.rm = T),
            `95%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.95, na.rm = T),
            .groups = "drop"
        ) %>%
        mutate(across(where(is.numeric), round, 4))
}

quantiles <- grouped_quantiles(clean_survey, est_year, season)

quantiles <- quantiles %>%
    pivot_longer(cols = 7:13, names_to = "quantile", values_to = "lat") %>%
    select(comname, est_year, season, quantile, lat) %>%
    group_by(comname, quantile) %>%
    mutate(rollmean = zoo::rollapplyr(lat, width = 5, FUN = mean, align = "center", partial = T))

quantiles %>%
    filter(quantile %in% c("5%", "95%")) %>%
    group_by(comname) %>%
    nest() %>%
    summarise(percentiles = map2(data, comname, function(x, y) {
        out <- ggplot(data = x) +
            geom_line(aes(x = est_year, y = rollmean, color = factor(quantile))) +
            geom_smooth(aes(x = est_year, y = rollmean, color = factor(quantile)), method = "lm", linetype = 2, alpha = 0.5, linewidth = 0.5, se = FALSE) +
            scale_color_manual(values = c("#ea4f12", "#00608a")) +
            ggtitle("Leading and trailing edge") +
            guides(color = guide_legend(title = "Percentile")) +
            ylab("Latitude (5-year rolling mean)") +
            xlab("Year") +
            facet_wrap(~season, ncol = 2) +
            guides(col = guide_legend(levels = c(`95%`, `5%`), title = "quantile")) +
            theme_gmri(
                strip.background = element_rect(fill = "transparent", linetype = 1, linewidth = 1, color = "transparent"),
                strip.text = element_text(color = "black", hjust = 0),
                plot.title = element_text(size = "12")
            )
        return(out)
    })) -> percentiles


library(grid)
library(gridExtra)

with_season <- with_season %>%
    left_join(percentiles) %>%
    left_join(dist_plot) %>%
    left_join(delta_sst) %>%
    left_join(delta_bt) %>%
    select(!data)

with_season %>%
    group_by(comname) %>%
    nest() %>%
    mutate(species_profile = map2(data, comname, function(x, y) {
        map <- x$map[[1]]
        times_series <- x$time_series[[1]]
        distance <- x$dist_plot[[1]]
        percentiles <- x$percentiles[[1]]
        sst <- x$sst_plot[[1]]
        bt <- x$bt_plot[[1]]

        out <- grid.arrange(map, times_series, distance, percentiles, sst, bt, nrow = 5, layout_matrix = rbind(c(1, 1), c(2, 3), c(4, 4), c(5, 6)), top = textGrob(paste(str_to_sentence(comname)), gp = gpar(col = "black", fontsize = 15, fontface = "bold")))
    })) -> species_profiles
