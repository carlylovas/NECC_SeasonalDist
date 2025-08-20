#####
## Miscellaneous
#####

library(tidyverse)

# Table of relevant results
seas_diff_slopes <- read_csv("Results/seasonal_dist_coefplot.csv") |>
    filter(.width == 0.95)

seas_trends <- read_csv(here::here("Results/seasonal_centroid_mod_summary.csv"))

#
library(dplyr)
library(stringr)
library(flextable)
library(officer)

# --- Standardize species and functional group ---
seas_diff_slopes2 <- seas_diff_slopes %>%
    mutate(comname = str_replace_all(species, "\\.", " "))

# Combine with seasonal trends
combined_results <- seas_diff_slopes2 %>%
    select(comname, functional_group = functional_group, slope_lower = .lower, slope_upper = .upper) %>%
    left_join(
        seas_trends %>%
            select(comname, fall_mean, spring_mean, gap_mean, shift_pattern, mechanism),
        by = "comname"
    ) %>%
    # Standardize functional_group names
    mutate(functional_group = case_when(
        functional_group %in% c("Benthivore", "benthivore", "benthos") ~ "benthivore",
        TRUE ~ functional_group
    )) %>%
    # Remove unwanted species
    filter(!comname %in% c("jonah crab", "scallop")) %>%
    # Add significance
    mutate(significance = case_when(
        slope_lower > 0 ~ "positive",
        slope_upper < 0 ~ "negative",
        TRUE ~ "ns"
    )) %>%
    # Round numeric columns
    mutate(across(c(slope_lower, slope_upper, fall_mean, spring_mean, gap_mean), ~ round(as.numeric(.x), 3))) %>%
    # Order by functional group then species alphabetically
    arrange(functional_group, comname) %>%
    # Ensure functional_group is first column
    select(
        functional_group, comname, slope_lower, slope_upper, significance,
        fall_mean, spring_mean, gap_mean, shift_pattern, mechanism
    )

# --- Build flextable ---
ft <- flextable(combined_results) %>%
    add_header_row(
        values = c(
            "Functional Group", "Species",
            "Changes in Distance Between Seasonal Centroids",
            "Trends in Seasonal Centroids"
        ),
        colwidths = c(1, 1, 3, 5)
    ) %>%
    set_header_labels(
        slope_lower = "Lower CI",
        slope_upper = "Upper CI",
        significance = "Significance",
        fall_mean = "Fall mean",
        spring_mean = "Spring mean",
        gap_mean = "Fall-Spring gap",
        shift_pattern = "Shift pattern",
        mechanism = "Mechanism"
    ) %>%
    merge_v(j = "functional_group") %>% # group rows
    border_outer(part = "all", border = fp_border(color = "black", width = 2)) %>%
    border_inner_h(part = "body", border = fp_border(color = "black", width = 1)) %>%
    border_inner_v(part = "all", border = fp_border(color = "black", width = 1)) %>%
    align(align = "center", part = "all") %>%
    autofit() %>%
    set_table_properties(width = 1, layout = "autofit") %>%
    align(align = "center", part = "all")

# --- Export to Word ---
doc <- read_docx() %>% body_add_flextable(ft)
print(doc, target = "combined_results.docx")
