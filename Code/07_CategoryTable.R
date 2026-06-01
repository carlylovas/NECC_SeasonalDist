library(tidyverse)
library(gmRi)
library(gt)
library(gtExtras)

# Mock data
## converge
converging <- data.frame(trend = "converging", year = seq(1970,2020, by = 10), spring = seq(38, 41, by = 0.6), fall = 42) %>%
  pivot_longer(cols = spring:fall, names_to = "season", values_to = "lat")

## march
marching <- data.frame(trend = "marching", year = seq(1970,2020, by = 10), spring = seq(38, 40, by = 0.4), fall = seq(42, 44, by = 0.4)) %>%
  pivot_longer(cols = spring:fall, names_to = "season", values_to = "lat")

## diverge
diverging <- data.frame(trend = "diverging", year = seq(1970,2020, by = 10), spring = seq(from = 38, to = 36, by = -0.4), fall = seq(42, 44, by = 0.4)) %>%
  pivot_longer(cols = spring:fall, names_to = "season", values_to = "lat")

## stable
stable <- data.frame(trend = "stable", year = seq(1970,2020, by = 10), spring = 38, fall = 42) %>%
  pivot_longer(cols = spring:fall, names_to = "season", values_to = "lat")

data <- converging |>
  full_join(marching) |>
  full_join(diverging) |>
  full_join(stable) |>
  group_by(trend) |>
  rbind(data.frame(trend = "unclear")) |>
  nest() |>
  mutate(stat_crit = case_when(
    trend == "stable" ~ "95% Bayesian credible intervals (CIs) for both spring and fall rates include zero.",
    trend == "marching" ~ "CIs for both spring and fall rates exclude zero; CI for the difference in rates (fall – spring) includes zero.",
    trend == "converging" ~ "CI for the difference in rates (fall - spring) excludes zero, and the difference in rates is negative.",
    trend == "diverging" ~ "CI for the difference in rates (fall - spring) excludes zero, and the difference in rates is positive.",
    trend == "unclear" ~ "CI for a single season rate excludes zero, but the difference in seasonal rates includes zero."
  )) |>
  mutate(interp = case_when(
    trend == "stable" ~ "No directional shift; seasonal centers of latitude remain consistent over time.",
    trend == "marching" ~ "The center of latitude for both seasons shifts in parallel at similar rates.",
    trend == "converging" ~ "Spring shifts northward faster than fall, narrowing the latitudinal range between seasons.",
    trend == "diverging" ~ "Fall shifts northward faster than spring, widening the latitudinal range between seasons.",
    trend == "unclear" ~ "One season shows a directional shift, whereas the other does not; the rates are not different, so it does not meet the criteria for convergence or divergence."
  )) |>
  mutate(plot = map(data, function(x,y){
    ggplot(data = x) +
      geom_line(aes(x = year, y = lat, color = season), linewidth = 2, arrow = arrow(ends = "last", angle = 45, length = unit(0.20, "cm"))) +
      geom_line(aes(x = year, y = lat, group = (year-0.5)), color = "#b94a40", linewidth = 2, alpha = 0.75, arrow = arrow(ends = "both", angle = 45, length = unit(0.20, "cm")), linetype = 2) +
      scale_color_manual(values = c("#00608a", "#ea4f12")) +
      labs(title = str_to_title(trend)) +
      xlim(c(1970,2020)) +
      scale_y_continuous(limits = c(36,44), breaks = c(36, 40, 44)) +
      theme(legend.position = "none",
            plot.title = element_text(size = 40, family = "Avenir", face = "bold"),
            line = element_blank(),
            rect = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(), 
            plot.margin = margin_auto(0))
  }))

data |>
  select(!data) |>
  select(!plot) |>
  gt(groupname_col = NULL) |>
  text_transform(
    locations = cells_body(columns = trend),
    fn = function(x) {
      data$plot[c(1:5)] |>
        ggplot_image(height = px(200))
    }
  ) |> 
  cols_label(
    trend = "",
    stat_crit = "Statistical criteria",
    interp = "Interpretation"
  ) |>
  # opt_table_font(font = list(google_font(name = "Avenir"))) |>
  tab_style(
    style = cell_borders(sides = "right", style = "solid", color = "lightgray", weight = px(2)),
    locations = cells_body(
      columns = trend)) |>
  tab_style(
    style = cell_text(weight = "bold", size = "x-large", font = "Avenir"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_text(size = "x-large",v_align = "top", font = "Avenir"),
    locations = cells_body(column = everything())
  ) -> table

gtsave(table, here::here("Figures", "statistical_criteria.png"))  


