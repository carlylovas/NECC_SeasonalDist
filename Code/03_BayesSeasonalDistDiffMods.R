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

fit <- FALSE

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

if (fit) { # Fit the model using brms (takes about an hour on my computer)
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

    write_rds(mod_brms, file = "Results/Fit_Mods/seasonal_dist_brms.rds", compress = "gz") # Write out the model file to disk.
    # mod_brms <- readRDS("Data/mod_brms1.rds") # Read it in here if you don't want to run the model
} else {
    mod_brms <- readRDS("Results/Fit_Mods/seasonal_dist_brms.rds") # Read in the model if you don't want to run it.
}

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

# predictions[1:10, 1:10] # gut check on structure

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
labeled_plot <- ggplot(nd, aes(x = year, y = ytilda)) +
    geom_line(aes(color = comname), show.legend = F) +
    geom_text_repel(
        data = label_df,
        aes(label = species_id),
        hjust = 0, # left-aligned
        direction = "y", # repel vertically
        nudge_x = 5, # push labels outside plot area
        size = 2,
        segment.color = "grey20",
        max.overlaps = 12,
        segment.color = "black",
        segment.linetype = "dotted", # 👈 controls the line style
        segment.size = 0.3
    ) +
    geom_line(data = nd2, aes(x = year, y = ytilda), linewidth = 1) +
    geom_ribbon(data = nd2, aes(x = year, y = ytilda, ymin = .lower, ymax = .upper), alpha = 0.1) +
    facet_wrap(~functional_group) +
    theme_bw() +
    labs(y = "Predicted distance between seasonal centroids", x = "") +
    coord_cartesian(xlim = c(1970, 2025), clip = "off") + # 🔑 allows labels to go outside the panel
    theme(plot.margin = margin(5.5, 50, 5.5, 5.5)) # 👈 adds right margin space


label_df <- label_df |>
    mutate(label_plot = paste0(species_id, ": ", comname)) |>
    group_by(functional_group) |>
    arrange(ytilda, .by_group = TRUE) |> # or arrange(ytilda) for ascending
    mutate(y_pos = row_number()) |>
    ungroup()

label_table <- ggplot(label_df, aes(x = 1, y = y_pos)) +
    geom_text(aes(label = label_plot, color = comname), hjust = 1, size = 3.2) +
    facet_wrap(~functional_group, scales = "free_y") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5),
        strip.text = element_blank()
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0))) # pull text to left

out <- labeled_plot / label_table + plot_layout(heights = c(1, 0.5))
ggsave("Figures/FG_temporaltrends.png", out, width = 11, height = 11, dpi = )

# This one I imagine going into the supplement, but it the same figure but faceted by species and includes the CI's for each species specific trend.
plot_list <- vector("list", length(unique(nd$functional_group)))
names(plot_list) <- unique(nd$functional_group)
colors <- c("#1b9e77", "#d95f02", "#7570b3")
ncol <- 3
nrow <- 6
total_panels <- ncol * nrow

i <- 1
for (fg in unique(nd$functional_group)) {
    # Filter data for this group
    nd_fg <- nd %>% filter(functional_group == fg)
    dat_mod_fg <- dat_mod %>% filter(functional_group == fg)

    # Real species in this group
    real_species <- unique(nd_fg$comname)
    n_real <- length(real_species)

    # Calculate how many dummy species needed
    n_dummy <- total_panels - n_real

    if (n_dummy > 0) {
        # Create dummy species names
        dummy_species <- paste0("dummy_", seq_len(n_dummy))

        # Create dummy rows with NA for all relevant columns (adjust as needed)
        dummy_rows <- tibble(
            functional_group = fg,
            comname = dummy_species,
            year = NA_real_,
            ytilda = NA_real_,
            .lower = NA_real_,
            .upper = NA_real_
        )

        # Combine real + dummy data
        nd_fg_full <- bind_rows(nd_fg, dummy_rows) %>%
            mutate(comname = factor(comname, levels = c(real_species, dummy_species)))

        # Set factor levels on dat_mod_fg for consistent ordering (points)
        dat_mod_fg <- dat_mod_fg %>%
            mutate(comname = factor(comname, levels = c(real_species, dummy_species)))
    } else {
        # No dummy needed, just set factor for order
        nd_fg_full <- nd_fg %>%
            mutate(comname = factor(comname, levels = real_species))

        dat_mod_fg <- dat_mod_fg %>%
            mutate(comname = factor(comname, levels = real_species))
    }

    # Build the plot
    plot_list[[fg]] <- ggplot(nd_fg_full, aes(x = year, y = ytilda)) +
        geom_point(
            data = dat_mod_fg, aes(x = year, y = dist_km, color = functional_group),
            size = 0.9, show.legend = FALSE, na.rm = TRUE
        ) +
        geom_line(
            data = nd_fg_full %>% filter(!is.na(year) & !is.na(ytilda)),
            aes(color = functional_group), show.legend = FALSE, na.rm = TRUE
        ) +
        geom_ribbon(
            data = nd_fg_full %>% filter(!is.na(year) & !is.na(.lower) & !is.na(.upper)),
            aes(ymin = .lower, ymax = .upper), alpha = 0.1, na.rm = TRUE
        ) +
        scale_color_manual(values = colors[i]) +
        scale_x_continuous(breaks = seq(from = 1970, to = 2023, by = 10)) +
        ylim(c(0, 600)) +
        facet_wrap(~comname, ncol = ncol, nrow = nrow) +
        labs(y = "Predicted distance between seasonal centroids", x = "") +
        ggtitle(fg) +
        theme_minimal(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5))

    i <- i + 1
}
plot_list$benthivore + plot_list$piscivore + plot_list$planktivore + plot_layout(guides = "collect", axes = "collect")
ggsave("Figures/Species_temporaltrends.png", width = 25, height = 15, dpi = 300)

ggplot(nd, aes(x = year, y = ytilda)) +
    geom_point(data = dat_mod, aes(x = year, y = dist_km, color = functional_group), size = 0.9, show.legend = F) +
    geom_line(aes(color = functional_group), show.legend = F) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.1) +
    facet_grid(functional_group ~ comname, scales = "free_y") +
    # facet_wrap(~functional_group + comname) +
    theme_minimal(base_size = 16) +
    labs(y = "Predicted distance between seasonal centroids", x = "")
ggsave("Figures/Species_temporaltrends.png")


# Coef plot: here, I do a little wrangling with the mod_brms posteriors to get at the coefficients (e.g the slope estimate [beta = year.cont*comname]) for the change in seasonal distance OVER TIME for each species nested within functional group.

out <- mod_brms %>%
    spread_draws(`r_functional_group:comname`[string, param_type]) %>% # Funky data wrangling. Unfortunately, this leaves a period in the common names which might curse us later, but going to avoid dealing with it for now.
    separate(string, into = c("functional_group", "species"), sep = "_") %>%
    filter(param_type == "year.cont")

# Build out the coeffiencient plot.
out_table <- out %>%
    group_by(species, functional_group) %>%
    tidybayes::median_qi(`r_functional_group:comname`, .width = c(0.75, 0.95))
out_table %>%
    ggplot(aes(x = `r_functional_group:comname`, y = forcats::fct_reorder(species, `r_functional_group:comname`))) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper)) +
    geom_vline(xintercept = 0, linetype = 3, color = "gray") +
    labs(x = "Change in seasonal distance over time", y = "Species") +
    theme_classic()
ggsave("Figures/Species_coefplot.png")

# Adding in average seasonal distance?
# Suppose you have avg spring-fall distances per species
avg_dist <- dat %>%
    group_by(comname) %>%
    summarize(mean_dist = mean(dist_km)) |>
    mutate(species = gsub(" ", ".", comname)) |>
    mutate(dist_group = cut(
        mean_dist,
        breaks = quantile(mean_dist, probs = c(0, 0.25, 0.75, 1), na.rm = TRUE),
        include.lowest = TRUE,
        labels = c("Low (0–25%)", "Medium (25–75%)", "High (75–100%)")
    )) # Replace spaces with periods to match the coef plot species names

# Join into slope/interval summary
out_table <- out_table %>%
    left_join(avg_dist)

# Plot slopes + intervals with color = average distance
cols <- c("#ece2f0", "#a6bddb", "#1c9099")
out_table %>%
    ggplot(aes(
        x = `r_functional_group:comname`,
        y = forcats::fct_reorder(species, `r_functional_group:comname`),
        color = dist_group
    )) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper)) +
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
        x = "Change in seasonal distance over time",
        y = "Species",
        color = "Avg. spring–fall\ndistance"
    ) +
    theme_classic()
ggsave("Figures/Species_coefplot_plusdist.png", height = 8, width = 11, dpi = 300)
# Save it
write.csv(out_table, "Results/seasonal_dist_coefplot.csv", row.names = FALSE)
