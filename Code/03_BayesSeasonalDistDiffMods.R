#####
## Fitting Bayesian Hierarchical Models to Species Seasonal Distance Differences by Functional Group
#####

library(brms)
library(tidyverse)
library(rstanarm)
library(glmmTMB)
library(tidybayes)
library(ggrepel)

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

soe_24 <- species_groupings |>
    select(comname, soe_24) |>
    mutate(across(where(is.character), tolower)) |>
    filter(comname %in% dat$comname) |>
    distinct()

# One missing?
unique(dat$comname)[!unique(dat$comname) %in% soe_24$comname]

add <- data.frame(comname = "windowpane flounder", soe_24 = "Benthivore")

soe_24 <- bind_rows(soe_24, add)

dat_mod <- dat |>
    left_join(soe_24) |>
    filter(!soe_24 == "Benthos")

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

# Fit the model using brms (takes about an hour on my computer)
mod_brms <- brm(
    dist_km ~ delta_year_day + (year.cont | functional_group / comname),
    data = dat_mod,
    family = Gamma(link = "log"),
    chains = 4, # Number of MCMC chains (adjust based on your computational resources)
    cores = 4, # Number of cores (adjust as needed)
    iter = 5000, # Number of iterations per chain
    warmup = 2000, # Number of warmup iterations
    # sample_prior = "only",
    control = list(adapt_delta = 0.999, max_treedepth = 15)
)

write_rds(mod_brms, file = "Data/seasonal_dist_brms.rds", compress = "gz") # Write out the model file to disk.
# mod_brms <- readRDS("Data/mod_brms1.rds") # Read it in here if you don't want to run the model

prior_summary(mod_brms) # This is the summary of the priors for the model.

mod_brms # There are still a few divergent transitions that I wasn't able to get rid of. Its possible that adjusting the priors might help, but that seems like quite a bit of work. I would advocate for running with this model for now.

pp_check(mod_brms, ndraws = 100) # seems reasonable to me!
shinystan::launch_shinystan(mod_brms) # This will launch a GUI to interrogate the model.

#--------------------------------------------------------
## Generate predictions (e.g. figures) from the model
#---------------------------------------------------------

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

predictions[1:10, 1:10] # gut check on structure

nd$ytilda <- apply(predictions, MARGIN = 2, median) # This code just extracts the median and 95% CI's for each value of year in the new data.
nd$.lower <- apply(predictions, 2, quantile, 0.025)
nd$.upper <- apply(predictions, 2, quantile, 0.975)

species_ids <- nd %>%
    distinct(comname) %>%
    arrange(comname) %>%
    mutate(species_id = row_number())
nd <- nd %>% left_join(species_ids, by = "comname")

nd <- nd |>
    filter(!functional_group %in% c("Benthivore", "benthos"))

label_df <- nd %>%
    group_by(comname, species_id, functional_group) %>%
    filter(year == max(year)) %>%
    ungroup()

# The first predictive data frame was based at the species level. Here we build on that only focused on the higher level of the heirarchy -- the functional group level.
nd2 <- expand.grid(
    delta_year_day = mean(dat_mod$delta_year_day, na.rm = T),
    year = unique(dat_mod$year),
    functional_group = unique(dat_mod$functional_group)
) %>%
    mutate(
        year.cont = (year - mean(year)) / 10,
        year_fac = as.character(year)
    )

predictions <- rstanarm::posterior_epred(mod_brms, newdata = nd2, re.form = ~ (year.cont | functional_group)) # Here we generate the prediction BUT we specifically estimate them based on the random slope and intercept of functional group NOT the species level variances.

nd2$ytilda <- apply(predictions, MARGIN = 2, median) # Same housekeeping to get medians and CI's.
nd2$.lower <- apply(predictions, 2, quantile, 0.025)
nd2$.upper <- apply(predictions, 2, quantile, 0.975)

nd2 <- nd2 |>
    filter(!functional_group %in% c("Benthivore", "benthos"))


# This is a figure paneled by functional group w/ only the CI's for functional group.
ggplot(nd, aes(x = year, y = ytilda)) +
    geom_line(aes(color = comname), show.legend = F) +
    geom_text_repel(
        data = label_df,
        aes(label = species_id),
        direction = "y",
        # nudge_x = 1,
        size = 1.5,
        segment.color = "grey80",
        max.overlaps = 12
    ) +
    geom_line(data = nd2, aes(x = year, y = ytilda), linewidth = 1) +
    geom_ribbon(data = nd2, aes(x = year, y = ytilda, ymin = .lower, ymax = .upper), alpha = 0.1) +
    facet_wrap(~functional_group) +
    theme_bw() +
    labs(y = "Predicted distance between seasonal centroids", x = "")

ggsave("Figures/FG_temporaltrends.png")

# This one I imagine going into the supplement, but it the same figure but faceted by species and includes the CI's for each species specific trend.
ggplot(nd, aes(x = year, y = ytilda)) +
    geom_point(data = df, aes(x = year, y = dist_km, color = comname), size = 0.9, show.legend = F) +
    geom_line(aes(color = comname), show.legend = F) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.1) +
    facet_wrap(~comname) +
    theme_bw() +
    labs(y = "Predicted distance between seasonal centroids", x = "")
ggsave("Figures/Species_temporaltrends.png")


# Coef plot: here, I do a little wrangling with the mod_brms posteriors to get at the coefficients (e.g the slope estimate [beta = year.cont*comname]) for the change in seasonal distance OVER TIME for each species nested within functional group.

out <- mod_brms %>%
    spread_draws(`r_functional_group:comname`[string, param_type]) %>% # Funky data wrangling. Unfortunately, this leaves a period in the common names which might curse us later, but going to avoid dealing with it for now.
    separate(string, into = c("functional_group", "species"), sep = "_") %>%
    filter(param_type == "year.cont")

# Build out the coeffiencient plot.
out %>%
    group_by(species, functional_group) %>%
    tidybayes::median_qi(`r_functional_group:comname`, .width = c(0.75, 0.95)) %>%
    ggplot(aes(x = `r_functional_group:comname`, y = forcats::fct_reorder(species, `r_functional_group:comname`))) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper)) +
    geom_vline(xintercept = 0, linetype = 3, color = "gray") +
    labs(x = "Change in seasonal distance over time", y = "Species") +
    theme_classic()
ggsave("Figures/Species_coefplot.png")
