#####
## Fitting Bayesian model to changes in latitude (longitude) overtime for each species by season
#####

# Libraries ----
library(brms)
library(tidyverse)
library(rstanarm)
library(glmmTMB)
library(tidybayes)
library(plotly)
library(ggrepel)

# Set up and data prep ----
fit_mod<- TRUE # Set to false if already done to avoid refitting

# Load data with biomass weighted cog lat/lon ----
dat <- readRDS(here::here("Data/seasonal_dist_withnas.rds")) # Generated in 02_SeasonalDistMetrics.R

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

# Check for missing
unique(dat$comname)[!unique(dat$comname) %in% soe_24$comname]
add <- data.frame(comname = c("lanternfish (unclassified)", "windowpane flounder"), soe_24 = c("planktivore", "benthivore"))
#add <- data.frame(comname = c("windowpane flounder"), soe_24 = c("benthivore"))

soe_24 <- bind_rows(soe_24, add)

dat_mod <- dat |>
    left_join(soe_24) |>
    filter(!comname == "sea scallop") |>
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

# Model fitting ----
# Going to fit a model to each species
dat_mod <- dat_mod |>
    group_by(comname) |>
    nest()

priors <- c(
  set_prior("normal(0, 5)", class = "b")
)

base_fit <- brm(
  avg_lat ~ year.cont * season,
  data = dat_mod$data[[1]],  # any one dataset is fine
  family = gaussian(),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 5000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

base_fit_lon <- brm(
  avg_lon ~ year.cont * season,
  data = dat_mod$data[[1]],  # any one dataset is fine
  family = gaussian(),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 5000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)


if (fit_mod) {
    dat_mod <- dat_mod |>
        mutate(
          model_lat = map(data, ~ update(
            base_fit,
            newdata = .x,
            recompile = FALSE
          ))) |>
          mutate(
            model_lon = map(data, ~ update(
              base_fit_lon,
              newdata = .x,
              recompile = FALSE
            )))
    write_rds(dat_mod, here::here("Results/Fit_Mods/seasonal_centroid_models_nas.rds"), compress = "gz")
} else {
    dat_mod <- readRDS(here::here("Results/Fit_Mods/seasonal_centroid_models_nas.rds"))
}

# Summarize Results ----
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
            # Marching: BOTH moving AND gap includes zero
            fall_ci_excludes_0 & spring_ci_excludes_0 & !gap_ci_excludes_0 ~ "Marching",
            # Converging: gap is negative (spring faster than fall)
            gap_ci_excludes_0 & summary_lat[[1]]$gap_mean < 0 ~ "Converging",
            # Diverging: gap is positive (fall faster than spring)
            gap_ci_excludes_0 & summary_lat[[1]]$gap_mean > 0 ~ "Diverging",
            # This catches cases where only one season moves but gap not significant
            # You might want a separate category for these
            TRUE ~ "Unclear"
        ),
        Dynamics = case_when(
            !fall_ci_excludes_0 & spring_ci_excludes_0 ~ "Spring moving only",
            fall_ci_excludes_0 & !spring_ci_excludes_0 ~ "Fall moving only",
            fall_ci_excludes_0 & spring_ci_excludes_0 ~ "Both moving",
            TRUE ~ "Unclear"
        )
    )

# Unnest and Plot
dat_plot <- dat_mod %>%
    select(comname, summary_lat, shift_pattern, Dynamics) %>%
    unnest(summary_lat) |>
    left_join(soe_24)

dat_plot <- dat_plot %>%
    mutate(
        shift_pattern = factor(shift_pattern, levels = c("Stable", "Marching", "Converging", "Diverging", "Unclear"))
    )

# Work on soe_24_clean labels
dat_plot <- dat_plot %>%
    mutate(
        `Functional Group` = case_when(
            soe_24 == "planktivore" ~ "Planktivore",
            soe_24 == "piscivore" ~ "Piscivore",
            soe_24 == "benthivore" ~ "Benthivore"
        ),
        Dyanmics = Dynamics
    )


# save it for later
write_csv(dat_plot, here::here("Results/seasonal_centroid_mod_summary_nas.csv"))

# Create formatted results table
results_table <- dat_plot %>%
    mutate(
        `Fall mean (LCI - UCI)` = sprintf("%.2f (%.2f - %.2f)", fall_mean, fall_lower, fall_upper),
        `Spring mean (LCI - UCI)` = sprintf("%.2f (%.2f - %.2f)", spring_mean, spring_lower, spring_upper),
        `Fall - Spring Difference (LCI - UCI)` = sprintf("%.2f (%.2f - %.2f)", gap_mean, gap_lower, gap_upper)
    ) %>%
    select(
        `Functional Group`,
        Species = comname,
        `Fall mean (LCI - UCI)`,
        `Spring mean (LCI - UCI)`,
        `Fall - Spring Difference (LCI - UCI)`,
        `Shift Pattern` = shift_pattern,
        Dynamics
    ) %>%
    arrange(`Functional Group`, Species)

# Save formatted table
write_csv(results_table, here::here("Results/seasonal_centroid_results_formatted_nas.csv"))

ggplot(dat_plot, aes(x = fall_mean, y = spring_mean, color = Dynamics, shape = `Functional Group`)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
    geom_errorbar(aes(ymin = spring_lower, ymax = spring_upper), alpha = 0.4) +
    geom_errorbarh(aes(xmin = fall_lower, xmax = fall_upper), alpha = 0.4) +
    geom_point(size = 2.5) +
    facet_wrap(~shift_pattern, nrow = 1) +
    labs(
        x = "Fall trend (° latitude/year)",
        y = "Spring trend (° latitude/year)",
        title = "Seasonal centroid shift patterns by species",
        subtitle = "With 95% credible intervals"
    ) +
    theme_minimal()
ggsave("Figures/SpeciesSeasonalTrends.jpg", width = 11, height = 8, dpi = 300)

# Highlihgting species
dat_plot$alpha_highlight <- ifelse(dat_plot$comname %in% c("alewife", "black sea bass"), 1, 0.6)

label_data <- subset(dat_plot, comname %in% c("alewife", "black sea bass"))

# Manually set label positions - adjust these coordinates as needed
label_data$label_x <- ifelse(label_data$comname == "alewife",
                              -0.185,  # alewife x position
                              0.22)  # black sea bass x position
label_data$label_y <- ifelse(label_data$comname == "alewife",
                              0.39,   # alewife y position
                              0.48)   # black sea bass y position


out<- ggplot(dat_plot, aes(x = fall_mean, y = spring_mean, color = Dynamics, shape = `Functional Group`, alpha = alpha_highlight)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
    geom_errorbar(aes(ymin = spring_lower, ymax = spring_upper)) +
    geom_errorbarh(aes(xmin = fall_lower, xmax = fall_upper)) +
    geom_point(data = subset(dat_plot, comname %in% c("alewife", "black sea bass")),
               size = 3.5, color = "black") +
    geom_point(size = 2.5) +
    ggrepel::geom_text_repel(
        data = label_data,
        aes(label = comname, x = label_x, y = label_y),
        size = 4,
        max.overlaps = 20,
        segment.color = "gray40",
        box.padding = 0.4,
        force = 0,
        force_pull = 0,
        show.legend = FALSE
    ) +
    facet_wrap(~shift_pattern, nrow = 1) +
    labs(
        x = "Fall trend (° latitude/year)",
        y = "Spring trend (° latitude/year)"
    ) +
    scale_alpha_identity() +
    theme_minimal(base_size = 16) +
    theme(legend.position = "bottom")
out
ggsave("Figures/SpeciesSeasonalTrends_CaseStudies.jpg", out, width = 15, height = 8, dpi = 300)

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
