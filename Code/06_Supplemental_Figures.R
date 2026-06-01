library(tidyverse)
library(matrixStats)
library(gridExtra)
library(here)
library(gmRi)
library(sf)

# Case study plots ----
dat  <- read_rds(here("Data","dat_clean.rds"))
dist <- read_rds(here("Data", "seasonal_dist.rds")) 
usa  <- rnaturalearth::ne_states(country = "united states of america", returnclass = "sf")
can  <- rnaturalearth::ne_states(country = "canada", returnclass = "sf")

## Centers & Percentiles
dat |>
  # filter(comname %in% c("alewife", "black sea bass")) |>
  group_by(comname, year, season) |>
  summarise(
    `5%`  = Hmisc::wtd.quantile(lat, weights = total_biomass_kg, probs = 0.05),
    `50%` = Hmisc::wtd.quantile(lat, weights = total_biomass_kg, probs = 0.50),
    `95%` = Hmisc::wtd.quantile(lat, weights = total_biomass_kg, probs = 0.95),
    cog_lat = weightedMean(lat, w = total_biomass_kg, na.rm = T),
    cog_lon = weightedMean(lon, w = total_biomass_kg, na.rm = T)
  ) |>
  mutate(decade = 10*year%/%10) -> dat_plot

dat_plot |>
  pivot_longer(cols = c(`5%`,`50%`,`95%`), names_to = "percentile", values_to = "lat") |>
  select(comname, year, season, percentile, lat) |>
  group_by(comname, percentile, season) |>
  mutate(rmean  = zoo::rollapplyr(lat, width = 5, FUN = mean, align = "center", partial = T), # double check with Andrew that he wants rolling means
         percentile = factor(percentile, levels = c('5%','50%','95%'))) -> percentiles

## Maps
dat_plot |>
  group_by(comname) |>
  nest() |>
  mutate(map = map2(data, comname, function(x,y){
    ggplot() +
      geom_sf(data = usa) + geom_sf(data = can) + coord_sf(xlim = c(-66,-77), ylim = c(36,46)) +
      geom_point(data = x, aes(x = cog_lon, y = cog_lat, color = season)) +
      scale_x_continuous(breaks = c(-76, -72, -68)) +
      scale_color_gmri() +
      facet_wrap(~decade, nrow = 1) + 
      guides(color = guide_legend(title = "Season")) +
      labs(title = "Center of biomass") +
      theme(text = element_text("Avenir"),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.box = "vertical",
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0, face = "plain", size = 10),
            panel.grid.major = element_line(color = "#535353", linewidth = 0.1, linetype = 3),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent"),
            panel.border = element_rect(fill = "transparent", linetype = 1, linewidth = 0.5, color = "#535353"))
  })) |>
  select(!data) -> maps

maps$map[[2]]

## Seasonal centers 
dat_plot |>
  group_by(comname) |>
  nest() |>
  mutate(plot = map2(data, comname, function(x,y){
    ggplot(data = x) +
      geom_line(aes(x = year, y = cog_lat, color = season), linetype = 2) +
      geom_smooth(aes(x = year, y = cog_lat, color = season), method = "lm", se = F) + 
      labs(title = "Biomass-weighted latitude percentiles",
           subtitle = "5-year rolling average") +
      # facet_wrap(~season, nrow = 1) +
      guides(color = guide_legend("Percentiles")) +
      scale_y_continuous(labels = ~ paste0(.x, "\u00B0N")) +
      scale_color_gmri() +
      theme(text = element_text(family = "Avenir"),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.box = "vertical", 
            strip.background = element_blank(), 
            strip.text = element_text(hjust = 0, face = "plain", size = 10),
            panel.grid.major = element_line(color = "#535353", linewidth = 0.1, linetype = 3),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent"),
            panel.border = element_rect(fill = "transparent", linetype = 1, linewidth = 0.5, color = "#535353"))
  })) |>
  select(!data) -> edges

## Percentile plots
# percentiles |>
#   group_by(comname) |>
#   nest() |>
#   mutate(plot = map2(data, comname, function(x,y){
#     ggplot(data = x) +
#       geom_line(aes(x = year, y = rmean, color = percentile), linetype = 2) +
#       geom_smooth(aes(x = year, y = rmean, color = percentile), method = "lm", se = F) + 
#       labs(title = "Biomass-weighted latitude percentiles",
#            subtitle = "5-year rolling average") +
#       facet_wrap(~season, nrow = 1) +
#       guides(color = guide_legend("Percentiles")) +
#       scale_y_continuous(labels = ~ paste0(.x, "\u00B0N")) +
#       scale_color_gmri() +
#       theme(text = element_text(family = "Avenir"),
#             axis.title = element_blank(),
#             legend.position = "bottom",
#             legend.box = "vertical", 
#             strip.background = element_blank(), 
#             strip.text = element_text(hjust = 0, face = "plain", size = 10),
#             panel.grid.major = element_line(color = "#535353", linewidth = 0.1, linetype = 3),
#             panel.grid.minor = element_blank(),
#             panel.background = element_rect(fill = "transparent"),
#             panel.border = element_rect(fill = "transparent", linetype = 1, linewidth = 0.5, color = "#535353"))
#   })) |>
#   select(!data) -> edges

edges$plot[[2]]

## Distance between centroids
dat_plot |>
  group_by(comname) |>
  nest() |>
  mutate(centers = map2(data, comname, function(x,y){
    ggplot(data = x) +
      geom_line(aes(x = year, y = cog_lat, group = year), color = "#535353", alpha = 0.8) +
      geom_point(aes(x = year, y = cog_lat, color = season)) +
      guides(color = guide_legend(title = "Season")) +
      labs(title = "Center of latitude",
           subtitle = "Biomass-weighted mean latitude") + 
      scale_y_continuous(labels = ~ paste0(.x, "\u00B0N")) +
      scale_color_gmri() +
      theme(text = element_text(family = "Avenir"),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.box = "vertical", 
            strip.background = element_blank(), 
            strip.text = element_text(hjust = 0, face = "plain", size = 10),
            panel.grid.major = element_line(color = "#535353", linewidth = 0.1, linetype = 3),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent"),
            panel.border = element_rect(fill = "transparent", linetype = 1, linewidth = 0.5, color = "#535353"))
      })) |>
  select(!data) -> centers

centers$centers[[1]]

## Combine 
all_plots <- maps |>
  full_join(edges) |>
  full_join(centers) |>
  group_by(comname) |>
  nest() |>
  mutate(wrap = map2(data, comname, function(x,y){
    map        <- x$map[[1]]
    center     <- x$centers[[1]]
    percentile <- x$plot[[1]]
    out <- gridExtra::grid.arrange(map, center, percentile, nrow = 2, layout_matrix=rbind(c(1,1), c(2,3)), 
                        top = grid::textGrob(paste(str_to_sentence(comname)), 
                                             gp = grid::gpar(col = "black", fontsize = 15, fontfamily = "Avenir", fontface = "bold"),
                                             just = "center"))
    # ggsave(plot = out, filename = here("Figures", "Supplemental", paste0(comname," supplemental.png")), width = 13, height = 10)
    })) -> case_studies
plot(case_studies$wrap[[2]])
