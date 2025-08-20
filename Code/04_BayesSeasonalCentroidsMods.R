#####
## Fitting Bayesian model to changes in latitude (longitude) overtime for each species by season
#####
library(brms)
library(tidyverse)
library(rstanarm)
library(glmmTMB)
library(tidybayes)
library(plotly)

# Load data with biomass weighted cog lat/lon ----
dat <- readRDS(here::here("Data/seasonal_dist.rds")) # Generated in 02_SeasonalDistMetrics.R

# Bring over the difference between seasonal median trawl survey day
delta_year_day <- readRDS(here::here("Data/trawl_yearday_diff.rds")) |>
    select(year, delta_year_day)

dat <- dat |>
    left_join(delta_year_day)

# Add in SOE species groups
# species_groupings <- ecodata::species_groupings |>
#   mutate(comname = tolower(COMNAME))
load(here::here("Data/SOE_species_list_24.RData"))
write_csv(species, here::here("Data/SOE_species_list_24.csv"))
species_groupings <- species |>
    tibble() |>
    janitor::clean_names()

soe_24 <- species_groupings |>
    select(comname, soe_24) |>
    mutate(across(where(is.character), tolower)) |>
    filter(comname %in% dat$comname) |>
    distinct()

table(soe_24$soe_24)
# One missing?
unique(dat$comname)[!unique(dat$comname) %in% soe_24$comname]

add <- data.frame(comname = "windowpane flounder", soe_24 = "benthivore")

soe_24 <- bind_rows(soe_24, add)

dat_mod <- dat |>
    left_join(soe_24) |>
    filter(!soe_24 == "benthos") |>
    mutate(season = factor(season, levels = c("Fall", "Spring")))
table(dat_mod$soe_24)


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

# Going to fit a model to each species
dat_mod <- dat_mod |>
    group_by(comname) |>
    nest()

fit_mod <- TRUE
if (fit_mod) {
    dat_mod <- dat_mod |>
        mutate(
            model_lat = map(data, ~ brm(
                avg_lat ~ year.cont * season,
                data = .x,
                family = gaussian(),
                chains = 4,
                cores = 4,
                iter = 5000,
                warmup = 1000,
                control = list(adapt_delta = 0.999, max_treedepth = 15)
            )),
            model_lon = map(data, ~ brm(
                avg_lon ~ year.cont * season,
                data = .x,
                family = gaussian(),
                chains = 4,
                cores = 4,
                iter = 5000,
                warmup = 1000,
                control = list(adapt_delta = 0.999, max_treedepth = 15)
            ))
        )
    saveRDS(dat_mod, here::here("Results/seasonal_centroid_models.rds"))
} else {
    dat_mod <- readRDS(here::here("Results/seasonal_centroid_models.rds"))
}

# Results summary
# Step 1: Extract draws
dat_mod <- dat_mod %>%
    mutate(
        draws_lat = map(model_lat, ~ .x %>%
            spread_draws(b_year.cont, `b_year.cont:seasonSpring`) %>%
            mutate(
                fall_trend = b_year.cont,
                spring_trend = b_year.cont + `b_year.cont:seasonSpring`,
                gap_change = fall_trend - spring_trend
            ))
    )


# Summarize credible intervals
# Step 2: Summarize draws with posterior intervals
dat_mod <- dat_mod %>%
    mutate(
        summary_lat = map(draws_lat, ~ .x %>%
            summarise(
                fall_mean = mean(fall_trend),
                fall_lower = quantile(fall_trend, 0.025),
                fall_upper = quantile(fall_trend, 0.975),
                spring_mean = mean(spring_trend),
                spring_lower = quantile(spring_trend, 0.025),
                spring_upper = quantile(spring_trend, 0.975),
                gap_mean = mean(gap_change),
                gap_lower = quantile(gap_change, 0.025),
                gap_upper = quantile(gap_change, 0.975)
            ))
    )

# Classify
dat_mod <- dat_mod %>%
    mutate(
        fall_ci_excludes_0 = map_lgl(summary_lat, ~ .x$fall_lower > 0 | .x$fall_upper < 0),
        spring_ci_excludes_0 = map_lgl(summary_lat, ~ .x$spring_lower > 0 | .x$spring_upper < 0),
        gap_ci_excludes_0 = map_lgl(summary_lat, ~ .x$gap_lower > 0 | .x$gap_upper < 0),
        shift_pattern = case_when(
            !fall_ci_excludes_0 & !spring_ci_excludes_0 ~ "Stable",
            !gap_ci_excludes_0 ~ "Same rate shift",
            gap_ci_excludes_0 & summary_lat[[1]]$gap_mean < 0 ~ "Converging",
            gap_ci_excludes_0 & summary_lat[[1]]$gap_mean > 0 ~ "Diverging",
            TRUE ~ "Unclear"
        ),
        mechanism = case_when(
            !fall_ci_excludes_0 & spring_ci_excludes_0 ~ "Spring moving only",
            fall_ci_excludes_0 & !spring_ci_excludes_0 ~ "Fall moving only",
            fall_ci_excludes_0 & spring_ci_excludes_0 ~ "Both moving",
            TRUE ~ "Unclear"
        )
    )

# Unnest and Plot
dat_plot <- dat_mod %>%
    select(comname, summary_lat, shift_pattern, mechanism) %>%
    unnest(summary_lat) |>
    left_join(soe_24)

dat_plot <- dat_plot %>%
    mutate(
        soe_24_clean = ifelse(tolower(soe_24) %in% c("benthivore", "benthos"), "benthivore", soe_24),
        shift_pattern = factor(shift_pattern, levels = c("Stable", "Same rate shift", "Converging", "Diverging"), labels = c("Stable", "Marching", "Converging", "Diverging"))
    )

# save it for later
write_csv(dat_plot, here::here("Results/seasonal_centroid_mod_summary.csv"))

ggplot(dat_plot, aes(x = fall_mean, y = spring_mean, color = mechanism, shape = soe_24_clean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
    geom_errorbar(aes(ymin = spring_lower, ymax = spring_upper), alpha = 0.4) +
    geom_errorbarh(aes(xmin = fall_lower, xmax = fall_upper), alpha = 0.4) +
    geom_point(size = 2.5) +
    facet_wrap(~shift_pattern) +
    labs(
        x = "Fall trend (° latitude/year)",
        y = "Spring trend (° latitude/year)",
        title = "Seasonal centroid shift patterns by species",
        subtitle = "With 95% credible intervals"
    ) +
    theme_minimal()
ggsave("Figures/SpeciesSeasonalTrends.jpg")

# Highlihgting species
dat_plot$alpha_highlight <- ifelse(dat_plot$comname %in% c("alewife", "longfin squid", "goosefish", "winter skate"), 1, 0.1)

t <- dat_plot |> filter(shift_pattern == "Diverging")

label_data <- subset(dat_plot, comname %in% c("alewife", "longfin squid", "goosefish", "winter skate"))
label_data$label_y <- label_data$spring_upper + 0.01 # Adjust offset as needed


ggplot(dat_plot, aes(x = fall_mean, y = spring_mean, color = mechanism, shape = soe_24_clean, alpha = alpha_highlight)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
    geom_errorbar(aes(ymin = spring_lower, ymax = spring_upper)) +
    geom_errorbarh(aes(xmin = fall_lower, xmax = fall_upper)) +
    geom_point(size = 2.5) +
    geom_text_repel(
        data = label_data,
        aes(y = label_y, label = comname), # Change to another label column if needed
        size = 4,
        max.overlaps = 20,
        segment.color = "gray40",
        box.padding = 0.4,
        show.legend = FALSE
    ) +
    facet_wrap(~shift_pattern) +
    labs(
        x = "Fall trend (° latitude/year)",
        y = "Spring trend (° latitude/year)",
        title = "Seasonal centroid shift patterns by species",
        subtitle = "With 95% credible intervals"
    ) +
    theme_minimal()
ggsave("Figures/SpeciesSeasonalTrends_CaseStudies.jpg")

# Interactive plot
# Base ggplot
p <- ggplot(dat_plot, aes(
    x = fall_mean, y = spring_mean,
    color = mechanism, shape = soe_24_clean,
    text = paste0(
        "Species: ", comname, "<br>",
        "Functional group: ", functional_group, "<br>",
        "Fall mean: ", fall_mean, "<br>",
        "Spring mean: ", spring_mean, "<br>",
        "Gap: ", gap_mean, "<br>",
        "Shift pattern: ", shift_pattern, "<br>",
        "Mechanism: ", mechanism
    )
)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
    geom_errorbar(aes(ymin = spring_lower, ymax = spring_upper), alpha = 0.4) +
    geom_errorbarh(aes(xmin = fall_lower, xmax = fall_upper), alpha = 0.4) +
    geom_point(size = 2.5) +
    facet_wrap(~shift_pattern) +
    labs(
        x = "Fall trend (° latitude/year)",
        y = "Spring trend (° latitude/year)",
        title = "Seasonal centroid shift patterns by species",
        subtitle = "With 95% credible intervals"
    ) +
    theme_minimal()

# Convert to interactive plotly
p_interactive <- ggplotly(p, tooltip = "text")

# Display interactive plot
p_interactive
