#####
## Code for figures
#####

# Libraries ----
library(here)
library(tidyverse)
library(gmRi)
library(patchwork)

# Figure 1: conceptual figure ----
## contracting
contracting <- data.frame(year = seq(1970,2020, by = 10), spring = seq(38, 41, by = 0.6), fall = 42) %>%
  pivot_longer(cols = spring:fall, names_to = "season", values_to = "lat")

ggplot(contracting) +
  geom_line(aes(x = year, y = lat, color = season), linewidth = 1, arrow = arrow(ends = "last", angle = 45, length = unit(0.20, "cm"))) +
  geom_line(aes(x = year, y = lat, group = (year-0.5)), color = "#b94a40", alpha = 0.75, arrow = arrow(ends = "both", angle = 45, length = unit(0.20, "cm")), linetype = 2) +
  scale_color_manual(values = c("#00608a", "#ea4f12")) +
  xlim(c(1970,2020)) +
  scale_y_continuous(limits = c(36,44), breaks = c(36, 40, 44)) +
  ggtitle("Contracting", subtitle = "The spring center of biomass is outpacing the fall:\n Seasonal distance in decreasing") + 
  ylab("Center of latitude") +xlab("Year") +
  theme_gmri(plot.subtitle = element_text(size = 10)) -> plot1

## march
marching <- data.frame(year = seq(1970,2020, by = 10), spring = seq(38, 40, by = 0.4), fall = seq(42, 44, by = 0.4)) %>%
  pivot_longer(cols = spring:fall, names_to = "season", values_to = "lat")

ggplot(marching) +
  geom_line(aes(x = year, y = lat, color = season), linewidth = 1, arrow = arrow(ends = "last", angle = 45, length = unit(0.20, "cm"))) +
  geom_line(aes(x = year, y = lat, group = (year-0.5)), color = "#b94a40", alpha = 0.75, arrow = arrow(ends = "both", angle = 45, length = unit(0.20, "cm")), linetype = 2) +
  scale_color_manual(values = c("#00608a", "#ea4f12")) +
  xlim(c(1970,2020)) +
  scale_y_continuous(limits = c(36,44), breaks = c(36, 40, 44)) +
  ggtitle("Marching", subtitle = "The spring center of biomass and fall center match pace and direction:\n Seasonal distance is stable") +
  ylab("Center of latitude") +xlab("Year") +
  theme_gmri(plot.subtitle = element_text(size = 10)) -> plot2

## expanding

expanding <- data.frame(year = seq(1970,2020, by = 10), spring = seq(from = 38, to = 36, by = -0.4), fall = seq(42, 44, by = 0.4)) %>%
  pivot_longer(cols = spring:fall, names_to = "season", values_to = "lat")

ggplot(expanding) +
  geom_line(aes(x = year, y = lat, color = season), linewidth = 1, arrow = arrow(ends = "last", angle = 45, length = unit(0.20, "cm"))) +
  geom_line(aes(x = year, y = lat, group = (year-0.5)), color = "#b94a40", alpha = 0.75, arrow = arrow(ends = "both", angle = 45, length = unit(0.20, "cm")), linetype = 2) +
  scale_color_manual(values = c("#00608a", "#ea4f12")) +
  xlim(c(1970,2020)) +
  scale_y_continuous(limits = c(36,44), breaks = c(36, 40, 44)) +
  ggtitle("Expanding", subtitle = "The spring center of biomass and fall center are moving away from each other:\n Seasonal distance is increasing") +
  ylab("Center of latitude") +xlab("Year") +
  theme_gmri(plot.subtitle = element_text(size = 10)) -> plot3

## stable

stable <- data.frame(year = seq(1970,2020, by = 10), spring = 38, fall = 42) %>%
  pivot_longer(cols = spring:fall, names_to = "season", values_to = "lat")

ggplot(stable) +
  geom_line(aes(x = year, y = lat, color = season), linewidth = 1, arrow = arrow(ends = "last", angle = 45, length = unit(0.20, "cm"))) +
  geom_line(aes(x = year, y = lat, group = (year-0.5)), color = "#b94a40", alpha = 0.75, arrow = arrow(ends = "both", angle = 45, length = unit(0.20, "cm")), linetype = 2) +
  scale_color_manual(values = c("#00608a", "#ea4f12")) +
  xlim(c(1970,2020)) +
  scale_y_continuous(limits = c(36,44), breaks = c(36, 40, 44)) +
  ggtitle("Stable", subtitle = "The spring center of biomass and fall center are not changing significantly: \n Seasonal distance is stable") +
  ylab("Center of latitude") +xlab("Year") +
  theme_gmri(plot.subtitle = element_text(size = 10)) -> plot4


patchwork::wrap_plots(plot1, plot2, plot3, plot4, ncol = 2)
# ggsave(...)

# Figure 2: functional groups ----
fg_trends <- read_rds(here("Data", "functional_group_trends.rds"))

fig2 <- ggplot(fg_trends) +
  # geom_line(data = nd, aes(x = year, y = ytilda, color = comname)) +
  geom_ribbon(aes(x = year, ymin = .lower, ymax = .upper, fill = functional_group), alpha = 0.25) +
  geom_line(aes(x = year, y = ytilda, color = functional_group), linewidth = 1) + 
  scale_color_gmri() + scale_fill_gmri() +
  # facet_wrap(~functional_group) +
  ylim(c(50, 470)) +
  xlab("Year") + ylab("Distance (km)") +
  ggtitle("Predicted distance between seasonal centroids") +
          # subtitle = "dist_km ~ delta_year_day + est_year*SOE_24 +\n(1|year_fac) + (1|comname), family = Gamma(link = `log`), REML = T)") +
  theme_gmri(plot.subtitle = element_text(size = 8),
             plot.title = element_text(size = 11),
             axis.text = element_text(size= 8),
             axis.title = element_text(size = 10))

# Figure 3: individual species response ----
bayes_coeff <- read_rds(here("Data", "bayes_coeff.rds"))

# minor edits to make plotting more neat
bayes_coeff <- bayes_coeff |>
  rename("estimate" = "r_functional_group:comname") |> 
  mutate(species = str_to_sentence(str_replace_all(species, "[.]", " ")), # removes . in the middle of species names
         z = estimate > 0) # will use for label placement

fig3 <- ggplot(bayes_coeff, aes(y = fct_reorder(species,estimate), 
                        x = estimate)) +
  geom_vline(xintercept = 0.0, color = "#535353", alpha = 0.7) +
  # geom_segment(aes(y = fct_reorder(species,estimate), x = 0, xend = estimate)) +
  geom_point(size = 2, alpha = 0.8) + 
  geom_pointinterval(aes(xmin = .lower, xmax = .upper))+
  ggtitle("Individual species response") +
  scale_color_gmri() + 
  theme_gmri(plot.title = element_text(size = 10),
             plot.subtitle = element_text(size = 9),
             # axis.text.y = element_blank(),
             axis.text.x = element_text(size = 9),
             axis.title.y = element_blank(),
             axis.title.x = element_text(size = 10)) 
  # annotate("label", y = (bayes_coeff %>% filter(z == F))$species, x = 0.001,
  #          label = (bayes_coeff %>% filter(z == F))$species, hjust = "left", size = 2.5, color = "white") +
  # annotate("text", y = (bayes_coeff %>% filter(z == F))$species, x = 0.001,
  #          label = (bayes_coeff %>% filter(z == F))$species, hjust = "left", size = 2.5) +
  # annotate("label", y = (bayes_coeff %>% filter(z == T))$species, x = -0.001,
  #          label = (bayes_coeff %>% filter(z == T))$species, hjust = "right", size = 2.5, color = "white") +
  # annotate("text", y = (bayes_coeff %>% filter(z == T))$species, x = -0.001,
  #          label = (bayes_coeff %>% filter(z == T))$species, hjust = "right", size = 2.5)


# Figure 4: movement categories ----
### load data 
seas_cob <- read_rds(here("Data", "seas_cob.rds")) 
mov <- read_rds(here("Data", "movement_categories.rds"))

fig4 <- seas_cob |>
  left_join(mov) |>
  group_by(comname) |>
  nest() |>
  mutate(plot = map2(data, comname, function(x,y){
    plot <- ggplot(data = x) +
      geom_line(aes(x = year, y = coastdist_km, group = year)) +
      geom_point(aes(x = year, y = coastdist_km, color = season)) +
      scale_color_gmri() +
      ggtitle(str_to_sentence(comname)) +
      xlab("Year") + ylab("Along shelf distance (km)") +
      theme_gmri()
    return(plot)
    })) # either select individuals, group by functional group, or arrange all for supplement

# Figure that I don't think works...
bayes_coeff |>
  left_join(mov |> mutate(species = str_to_sentence(comname))) |>
  ggplot(aes(y = fct_reorder(species,estimate),
             x = estimate,
             color = category)) +
  geom_vline(xintercept = 0.0, color = "#535353", alpha = 0.7) +
  geom_point(size = 2, alpha = 0.8) + 
  geom_pointinterval(aes(xmin = .lower, xmax = .upper))+
  ggtitle("Individual species response") +
  scale_color_gmri() + 
  theme_gmri(plot.title = element_text(size = 10),
             plot.subtitle = element_text(size = 9),
             axis.text.x = element_text(size = 9),
             axis.title.y = element_blank(),
             axis.title.x = element_text(size = 10)) 
