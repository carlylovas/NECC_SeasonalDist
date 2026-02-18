#####
## Fitting Bayesian Hierarchical Models to Species Seasonal Distance Differences by Functional Group
#####

library(brms)
library(tidyverse)
library(rstanarm)
library(glmmTMB)
library(tidybayes)
library(ggrepel)
library(patchwork)
library(ggh4x)

fit <- FALSE # Set to TRUE to fit the model, starts at line 67, FALSE to read in the model
#fit <- TRUE

# Load seasonal distance data ----
dat <- readRDS(here::here("Data/seasonal_dist.rds")) # Generated in 02_SeasonalDistMetrics.R

# Duplicated seasonal distances right now, some filtering
dat <- dat |>
  distinct(comname, year, dist_km)

# Bring over the difference between seasonal median trawl survey day
delta_year_day <- readRDS(here::here("Data/trawl_yearday_diff.rds")) |>
  select(year, delta_year_day)

dat <- dat |>
  left_join(delta_year_day)

# Add in SOE species groups
# species_groupings <- ecodata::species_groupings |>
#   mutate(comname = tolower(COMNAME))
load(here::here("Data/SOE_species_list_24.RData"))
species_groupings <- species |>
  tibble() |>
  janitor::clean_names()

# Write this out for future use
write.csv(species_groupings, here::here("Data/SOE_species_list_24.csv"), row.names = FALSE)

soe_24 <- species_groupings |>
  select(comname, soe_24) |>
  mutate(across(where(is.character), tolower)) |>
  filter(comname %in% dat$comname) |>
  distinct()

# One missing?
unique(dat$comname)[!unique(dat$comname) %in% soe_24$comname]

add <- data.frame(comname = "windowpane flounder", soe_24 = "benthivore")

soe_24 <- bind_rows(soe_24, add)

dat_mod <- dat |>
  left_join(soe_24) |>
  filter(!soe_24 == "benthos") # Sea scallops -- removing these for now.

# Some final naming/tidying
dat_mod <- dat_mod |>
  rename(functional_group = soe_24) |>
  mutate(
    year_fac = factor(year),
    year.cont = (year - mean(year)) / 10, # This converts year (~1970-1923) into a scaled variable where a 1 unit change in the variable represents the shift over 1 decade.
    delta_year_day = scale(delta_year_day) # This just scales the delta_year_day variable to that everything is on similar scales.
  )
summary(dat_mod)
t <- dat_mod |> filter(is.na(dist_km))
t




mod_brms <- readRDS("Results/Fit_Mods/seasonal_dist_brms.rds")


# # Coef plot: here, I do a little wrangling with the mod_brms posteriors to get at the coefficients (e.g the slope estimate [beta = year.cont*comname]) for the change in seasonal distance OVER TIME for each species nested within functional group.
# 
# out <- mod_brms %>%
#   spread_draws(`r_functional_group:comname`[string, param_type]) %>% # Funky data wrangling. Unfortunately, this leaves a period in the common names which might curse us later, but going to avoid dealing with it for now.
#   separate(string, into = c("functional_group", "species"), sep = "_") %>%
#   filter(param_type == "year.cont")
# 
# # Build out the coeffiencient plot.
# out_table <- out %>%
#   group_by(species, functional_group) %>%
#   tidybayes::median_qi(`r_functional_group:comname`, .width = c(0.75, 0.95))
# 
# # Common name to sentence case for plotting purposes
# out_table <- out_table %>%
#   mutate(
#     comname = gsub("\\.", " ", species),
#     comname = str_to_sentence(comname)
#   )
# 
# # No colors
# out_table %>%
#   ggplot(aes(x = `r_functional_group:comname`, y = forcats::fct_reorder(species, `r_functional_group:comname`))) +
#   geom_pointinterval(aes(xmin = .lower, xmax = .upper)) +
#   geom_vline(xintercept = 0, linetype = 3, color = "gray") +
#   labs(x = "Change in seasonal distance over time", y = "Species") +
#   theme_classic()
# ggsave("Figures/Species_coefplot.png")

post_coef <- coef(mod_brms, summary = FALSE)
species_slopes <- post_coef$`functional_group:comname`[, , "year.cont"] %>%
  as.data.frame() %>%
  pivot_longer(
    cols = everything(),
    names_to = "species",
    values_to = "beta_species"
  ) %>%
  group_by(species) %>%
  median_qi(beta_species, .width = c(0.75, 0.95)) %>%
  mutate(
    pct_change = (exp(beta_species) - 1) * 100,
    pct_lower  = (exp(.lower) - 1) * 100,
    pct_upper  = (exp(.upper) - 1) * 100
  ) %>%
  separate(species, into = c("fg", "species"), sep = "_")

ggplot(species_slopes, aes(x = pct_change, y = forcats::fct_reorder(species, pct_change))) +
  geom_pointinterval(aes(xmin = pct_lower, xmax = pct_upper)) +
  geom_vline(xintercept = 0, linetype = 3, color = "gray") +
  labs(x = "Percent change in seasonal distance per decade", y = "Species") +
  theme_classic()
ggsave("Figures/Species_coefplot_correct.png")



post_coef$`functional_group`[, , "year.cont"] %>%
  as.data.frame() %>%
  pivot_longer(
    cols = everything(),
    names_to = "fg",
    values_to = "beta_fg"
  ) %>%
  group_by(fg) %>%
  median_qi(beta_fg, .width = c(0.75, 0.95))



# Adding in average seasonal distance?
# Suppose you have avg spring-fall distances per species

avg_dist <- dat %>%
  mutate(
    comname = str_to_sentence(comname)
  ) %>%
  group_by(comname) %>%
  summarize(mean_dist = mean(dist_km)) |>
  mutate(dist_group = cut(
    mean_dist,
    breaks = quantile(mean_dist, probs = c(0, 0.25, 0.75, 1), na.rm = TRUE),
    include.lowest = TRUE,
    labels = c("Low (0–25%)", "Medium (25–75%)", "High (75–100%)")
  )) # Replace spaces with periods to match the coef plot species names

# Join into slope/interval summary
out_table <- species_slopes %>%
  mutate(comname = str_to_sentence(species)) %>%
  left_join(avg_dist)

# Plot slopes + intervals with color = average distance
cols <- c("#ece2f0", "#a6bddb", "#1c9099")

p_main <- out_table %>%
  ggplot(aes(
    x = pct_change,
    y = forcats::fct_reorder(species, pct_change),
    color = dist_group
  )) +
  geom_pointinterval(aes(xmin = pct_lower, xmax = pct_upper)) +
  geom_vline(xintercept = 0, linetype = 3, color = "gray") +
  scale_color_manual(
    values = c(
      "Low (0–25%)" = cols[1],
      "Medium (25–75%)" = cols[2],
      "High (75–100%)" = cols[3]
    ),
    labels = c(
      "Low (28 – 65 km)",
      "Medium (66 – 160 km)",
      "High (161 – 358 km)"
    )
  ) +
  labs(
    x = "Percent change in seasonal distance per decade",
    y = "Species",
    color = "Avg. spring–fall\ndistance"
  ) +
  theme_classic()+
  theme(legend.position = "inside", 
        legend.position.inside = c(0.8, 0.2))
ggsave("Figures/Species_coefplot_plusdist.png", plot = p_main, height = 8, width = 11, dpi = 300)
# Save it
write.csv(out_table, "Results/seasonal_dist_coefplot.csv", row.names = FALSE)


formerge <- dat_mod |>
  distinct(comname, functional_group) # Housekeeping to generate the new data to predict on.

nd <- expand.grid(
  delta_year_day = mean(dat_mod$delta_year_day, na.rm = T),
  year = unique(dat_mod$year),
  # functional_group = unique(df$functional_group),
  comname = unique(dat_mod$comname)
) |>
  mutate(
    year.cont = (year - mean(year)) / 10,
    year_fac = as.character(year)
  ) |>
  left_join(formerge) # all this is just generating a data frame to predict on.

predictions <- rstanarm::posterior_epred(mod_brms, newdata = nd, re.form = NULL) # Here we are using the model to generate predictions. Specifically, this will generate LINEAR predictions and will have smaller variance that posterior_predict(), because the variance is only based on the uncertainty in the expected value of the posterior predictive distribution. The residual error is ignored. For true posterior PREDICTIVE intervals we could use posterior_predict(). re.form = NULL ensures that the model accounts for ALL levels of uncertainly, e.g. uncertainly due to both fixed and random effects.

# predictions[1:10, 1:10] # gut check on structure

nd$ytilda <- apply(predictions, MARGIN = 2, median) # This code just extracts the median and 95% CI's for each value of year in the new data.
nd$.lower <- apply(predictions, 2, quantile, 0.025)
nd$.upper <- apply(predictions, 2, quantile, 0.975)

species_ids <- nd %>%
  distinct(comname) %>%
  arrange(comname) %>%
  mutate(species_id = row_number())
nd <- nd %>% left_join(species_ids, by = "comname")


colors <- c("#66c2a5", "#fc8d62", "#8da0cb")

p_side <- nd %>% 
  filter(comname %in% c("alewife", "longfin squid", "jonah crab")) %>%
  ggplot(aes(x = year, y = ytilda)) +
  geom_point(
    data = dat_mod %>%
      filter(comname %in%c("alewife", "longfin squid", "jonah crab")),
    aes(x = year, y = dist_km, color = functional_group),
    size = 0.9, show.legend = FALSE, na.rm = TRUE
  ) +
  geom_line(
    aes(color = functional_group), show.legend = FALSE, na.rm = TRUE
  ) +
  geom_ribbon(
    aes(ymin = .lower, ymax = .upper), alpha = 0.1, na.rm = TRUE
  ) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = seq(from = 1970, to = 2023, by = 10)) +
  ylim(c(0, 600)) +
  facet_wrap(~comname, nrow = 3) +
  labs(y = "Distance (km)", x = "") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5))

p_final <- cowplot::plot_grid(p_main, p_side, rel_widths = c(0.7, 0.3))
ggsave("Figures/Fig3_wspecies.png", plot = p_final, height = 8, width = 12, dpi = 300, bg = "white")

