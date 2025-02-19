# Conceptual figure
## contract
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


patch <- patchwork::wrap_plots(plot1, plot2, plot3, plot4, ncol = 2)
ggsave(here("Figures", "concept_fig.png"), patch, width = 11, height = 8.5, units = "in")
