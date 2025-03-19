#####
## Modeling difference in seasonal distance
#####

# Libraries ----
library(here)
library(tidyverse)
library(broom)
library(gt)
library(gmRi)
library(matrixStats)
library(ecodata)
library(glmmTMB)
library(marginaleffects)
library(ggeffects)
library(emmeans)

# Load seasonal distance data ----
dat <- readRDS(here::here("Data/seasonal_dist.rds"))

# Duplicated seasonal distances right now, some filtering
dat <- dat |>
    distinct(comname, year, dist_km)

# Bring over the difference between seasonal median trawl survey day
delta_year_day <- readRDS(here::here("Data/trawl_yearday_diff.rds")) |>
    select(year, delta_year_day)

dat <- dat |>
    left_join(delta_year_day)

# Add in SOE species groups
species_groupings <- ecodata::species_groupings |> 
  mutate(comname = tolower(COMNAME))

soe_24 <- species_groupings |>
    select(comname, SOE.24) |>
    filter(comname %in% dat$comname) |>
    distinct()

# One missing?
unique(dat$comname)[!unique(dat$comname) %in% soe_24$comname]

add <- data.frame(comname = "windowpane flounder", SOE.24 = "Benthivore")

soe_24 <- bind_rows(soe_24, add)

dat_mod <- dat |>
    left_join(soe_24) |>
    filter(!SOE.24 == "Benthos") %>%
    mutate(year_fac = factor(year))

# add movement categories
mov <-  readRDS(here::here("Data/movement_categories.rds"))

# GLMM by group, RE for species ----
mod_group <- glmmTMB(dist_km ~ delta_year_day + year * SOE.24 + (1 | year_fac) + (1 | comname), data = dat_mod, family = Gamma(link = "log"))

m <- as.data.frame(ggpredict(mod_group, terms = c("year", "SOE.24")))

ggplot(data = m) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  geom_line(aes(x = x, y = predicted, color = group), linewidth = 1) + 
  scale_color_gmri() + scale_fill_gmri() +
  ylim(c(50, 470)) +
  xlab("Year") + ylab("Predicted (dist_km)") +
  guides(color = guide_legend(title = "Functional group"), fill = guide_legend(title = "Functional group")) +
  ggtitle("Predicted values of dist_km",
          subtitle = "dist_km ~ delta_year_day + est_year*SOE_24 +\n(1|year_fac) + (1|comname), family = Gamma(link = `log`), REML = T)") +
  theme_gmri(plot.subtitle = element_text(size = 8),
             plot.title = element_text(size = 11),
             axis.text = element_text(size= 8),
             axis.title = element_text(size = 10)) -> p1

ggsave(here("Figures", "func_groups.png"), p1, height = 5, width = 7, units = "in", bg = "white")

# GLMM by species ----
mod_species <- glmmTMB(dist_km ~ delta_year_day + year*comname + (1|year_fac), data = dat_mod, family=Gamma(link="log"), REML = T)

summary <- summary(mod_species)[["coefficients"]]
coeff <- as.data.frame(summary[["cond"]]) %>% 
  rownames_to_column() %>%
  mutate(signif = ifelse(`Pr(>|z|)` < 0.05, T, F)) %>%
  filter(rowname %in% c(str_subset(rowname, pattern = "year"))) %>%
  mutate(comname = str_remove(rowname, pattern = "year:comname"))

species_trends <- as.data.frame(emtrends(mod_species, "comname", "year")) |>
    left_join(soe_24) |>
    mutate(
        comname = str_to_sentence(comname),
        color = year.trend > 0.000
    )


# Trying out the slopes function in marginal effects function
slopes_output <- slopes(mod_species, variables = "year", by = "comname", type = "response")
# print(slopes_output)

spp_group<- dat_mod |>
  dplyr::select(comname, SOE.24) |>
  mutate(comname = str_to_lower(comname)) |>
  distinct()
slopes_output<- data.frame(slopes_output) |>
  left_join(spp_group) |>
  left_join(mov) |>
  mutate(z = estimate > 0,
         comname = str_to_sentence(comname)) 


## plot
ggplot(slopes_output) +
  geom_vline(xintercept = 0.0, color = "#535353", alpha = 0.7) +
  geom_segment(aes(y = fct_reorder(comname, estimate), x = 0, xend = estimate)) +
  geom_point(aes(y = fct_reorder(comname,estimate), x = estimate, color = category), size = 2, alpha = 0.8) + 
  xlab("slope (units?)") + 
  # xlim(c(-7,7)) +
  ggtitle("Individual species response", subtitle = "Linear mixed effects model") +
  scale_color_gmri() + 
  theme_gmri(plot.title = element_text(size = 10),
             plot.subtitle = element_text(size = 9),
             # legend.position = "none",
             axis.text.y = element_blank(),
             axis.text.x = element_text(size = 9),
             axis.title.y = element_blank(),
             axis.title.x = element_text(size = 10)) +
             # panel.grid.major = element_blank(),
             # panel.grid.minor = element_blank()) +
  guides(color = guide_legend(title = "Centroid movement")) +
  annotate("label", y = (slopes_output %>% filter(z == F))$comname, x = 0.1,
           label = (slopes_output %>% filter(z == F))$comname, hjust = "left", size = 2.5, color = "white") +
  annotate("text", y = (slopes_output %>% filter(z == F))$comname, x = 0.1,
           label = (slopes_output %>% filter(z == F))$comname, hjust = "left", size = 2.5) +
  annotate("label", y = (slopes_output %>% filter(z == T))$comname, x = -0.1,
           label = (slopes_output %>% filter(z == T))$comname, hjust = "right", size = 2.5, color = "white") +
  annotate("text", y = (slopes_output %>% filter(z == T))$comname, x = -0.1,
           label = (slopes_output %>% filter(z == T))$comname, hjust = "right", size = 2.5) -> p2
## I'm very unsure if this works in the way we want it to
## the categories seem to contradict what we are saying is happening to the seasononal distance 
## comparing apples to oranges? lm to glmm? significance testing? 

ggsave(here("Figures", "slope_plot.png"), p2, height = 6, width = 6, units = "in", bg = "white")
