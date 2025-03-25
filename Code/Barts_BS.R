library(tidyverse)
# library(rstanarm)
library(glmmTMB)
library(tidybayes)
library(rstanarm)


df <- readRDS("Data/mod_Bart.rds") %>%
  rename(functional_group = SOE.24) %>%
  ungroup() %>%
  mutate(year.cont = (year - mean(year)) / 10, 
         delta_year_day = scale(delta_year_day))

summary(df)

# mod1 <- rstanarm::stan_glmer(dist_km ~ delta_year_day + year.cont + (year.cont|functional_group/comname) + (1|year_fac), data = df, family = Gamma(link = "log")) # This model had pretty considerable issues. Dropping the (1|year_fac) random effect to simplify. 

mod2 <- rstanarm::stan_glmer(dist_km ~ delta_year_day + year.cont + (year.cont|functional_group/comname), 
                             data = df,
                             family = Gamma(link = "log"), 
                             QR = T, 
                             adapt_delta = 0.99,
                             chains = 4,               # Number of chains
                             cores = 4,                # Number of cores (use a number <= the number of chains)
                             iter = 2000,              # Number of iterations
                             warmup = 1000)             # Number of warmup iterations) # This is a similar model, however I dropped the random intercept effect of year. Also I turned on the QR argument, which applies a scaled qr decomposition to the design matrix. I don't understand the math here, but the documentation recommends to have this set to True w/ multiple predictors. 

mod3 <- rstanarm::stan_glmer(dist_km ~ delta_year_day + year.cont + (year.cont|comname/functional_group), 
                             data = df,
                             family = Gamma(link = "log"), 
                             QR = T, 
                             adapt_delta = 0.99,
                             chains = 4,               # Number of chains
                             cores = 4,                # Number of cores (use a number <= the number of chains)
                             iter = 2000,              # Number of iterations
                             warmup = 1000) 

write_rds(mod2, "Data/mod2.rds", compress = "gz") # Write it out because of model fitting time.

mod2 <- readRDS("Data/mod2.rds") # Read it in from file. 

shinystan::launch_shinystan(mod2) # This opens a GUI to explore the model and check for model convergence. In particulary it can generate posterior predictive checks of the model output.

summary(mod2)
print(summary(mod2), digits = 3)
get_variables(mod2)




mod2 %>%
  spread_draws(`b[year.cont functional_group:Benthivore]`,
               `b[year.cont functional_group:Piscivore]`, 
               `b[year.cont functional_group:Planktivore]`) %>% 
  rename(beta_benthivores = `b[year.cont functional_group:Benthivore]`, 
         beta_piscivores = `b[year.cont functional_group:Piscivore]`, 
         beta_planktivores = `b[year.cont functional_group:Planktivore]`) %>%
  pivot_longer(cols = beta_benthivores:beta_planktivores) %>%
  group_by(name) %>%
  # median_qi(value) %>%
  ggplot(aes(x = value, y = name))+
  geom_halfeyeh()

mod2 %>%
  spread_draws(`b[year.cont functional_group:Benthivore]`,
               `b[year.cont functional_group:Piscivore]`, 
               `b[year.cont functional_group:Planktivore]`) %>% 
  rename(beta_benthivores = `b[year.cont functional_group:Benthivore]`, 
         beta_piscivores = `b[year.cont functional_group:Piscivore]`, 
         beta_planktivores = `b[year.cont functional_group:Planktivore]`) %>%
  pivot_longer(cols = beta_benthivores:beta_planktivores) %>%
  group_by(name) %>%
  median_qi(value)



df %>%
  modelr::data_grid(year.cont = modelr::seq_range(year.cont, n = 100),
                    comname = "alewife", 
                    functional_group = "Planktivore",
                    delta_year_day = mean(delta_year_day, na.rm = T),
                    .model = mod2) %>%
  add_epred_draws(mod2, ndraws = 100) %>%
  ggplot(aes(x = year.cont)) +
  stat_lineribbon(aes(y = .epred), .width = c( 0.95), fill = alpha("gray80", 0.5), show.legend = F)




plot(ggeffects::ggpredict(mod2, terms = ~functional_group))


formerge <- df %>% 
  distinct(comname, functional_group)

nd <- expand.grid(delta_year_day = mean(df$delta_year_day, na.rm = T), 
                  year = unique(df$year), 
                  # functional_group = unique(df$functional_group), 
                  comname = unique(df$comname)) %>%
  mutate(year.cont = (year - mean(year)) / 10, 
         year_fac = as.character(year)) %>%
  left_join(formerge)

predictions <- rstanarm::posterior_epred(mod2, newdata = nd)
predictions[1:10,1:10]

nd$ytilda <- apply(predictions, MARGIN = 2, median)
nd$.lower <- apply(predictions, 2, quantile,  0.025)
nd$.upper <- apply(predictions, 2, quantile, 0.975)

nd2 <- expand.grid(delta_year_day = mean(df$delta_year_day, na.rm = T), 
                  year = unique(df$year), 
                  functional_group = unique(df$functional_group)) %>%
  mutate(year.cont = (year - mean(year)) / 10, 
         year_fac = as.character(year))

predictions <- rstanarm::posterior_epred(mod2, newdata = nd2, re.form = ~0)

nd2$ytilda <- apply(predictions, MARGIN = 2, median)
nd2$.lower <- apply(predictions, 2, quantile,  0.025)
nd2$.upper <- apply(predictions, 2, quantile, 0.975)

ggplot(nd, aes(x = year, y = ytilda))+
  geom_line(aes(color = comname))+
  geom_line(data = nd2, aes(x = year, y = ytilda))+
  geom_ribbon(data = nd2, aes(x = year, y = ytilda, ymin = .lower, ymax = .upper), alpha = 0.1)+
  facet_wrap(~functional_group)+
  theme_bw()

rstanarm::prior_summary(mod2)

ggplot(nd2, aes(x = year, y = ytilda))+
  geom_line(aes(color = functional_group))+
  facet_wrap(~functional_group)

#------------------------------------------
## Frequentist version
#------------------------------------------

# mod1 <- glmmTMB::glmmTMB(dist_km ~ delta_year_day + year.cont + (year.cont|functional_group/comname) + (1|year_fac), data = df, family = Gamma(link = "log"))
# summary(mod1)

mod1 <- glmmTMB::glmmTMB(dist_km ~ delta_year_day + year.cont + (year.cont|comname/functional_group), data = df, family = Gamma(link = "log"))
summary(mod1)

coef(mod1)

# effects <- ranef()

nd <- expand.grid(delta_year_day = mean(df$delta_year_day, na.rm = T), 
                       year = unique(df$year), 
                       # functional_group = unique(df$functional_group), 
                       comname = unique(df$comname)[1:6]) %>%
  mutate(year.cont = (year - mean(year)) / 10, 
         year_fac = as.character(year), 
         functional_group = case_when(comname == "acadian redfish" ~ "Piscivore", 
                                      comname == "alewife" ~ "Planktivore", 
                                      comname == "american lobster" ~ "Benthivore", 
                                      comname == "american plaice" ~ "Benthivore", 
                                      comname == "american shad" ~ "Planktivore", 
                                      comname == "atlantic cod" ~ "Piscivore"))
summary(nd)


nd$predictions <- predict(mod1, newdata = nd, re.form = NULL, type = "response")


df %>%
  filter(comname %in% unique(nd$comname)) %>%
  ggplot(aes(x = year, y = dist_km))+
  geom_point(aes(color = comname))+
  geom_line(data = nd, aes(x = year, y = predictions, color = comname))+
  facet_wrap(~functional_group)


out <- ggeffects::ggpredict(mod1, term = c("year.cont", "functional_group"))
plot(out)









library(brms)

# Fit the model using brms
mod_brms <- brm(
  dist_km ~ delta_year_day + (year.cont | functional_group/comname),
  data = df,
  family = Gamma(link = "log"),
  chains = 4,   # Number of MCMC chains (adjust based on your computational resources)
  cores = 4,    # Number of cores (adjust as needed)
  iter = 5000,  # Number of iterations per chain
  warmup = 2000, # Number of warmup iterations
  # sample_prior = "only",
  control = list(adapt_delta = 0.999, max_treedepth = 15)
  )

prior_summary(mod_brms)
get_prior(mod_brms)
mod_brms

post_samples <- as_draws(mod_brms)
plot(mod_brms, variable = "r_functional_group[Planktivore,year.cont]")
plot(mod_brms, variable = "r_functional_group[Benthivore,year.cont]")
plot(mod_brms, variable = "r_functional_group[Piscivore,year.cont]")

plot(ggeffects::ggpredict(mod_brms, terms = ~year.cont*functional_group))

mod_brms %>%
  gather_draws()

predictions <- rstanarm::posterior_epred(mod_brms, newdata = nd, re.form = NULL)
predictions[1:10,1:10]

nd$ytilda <- apply(predictions, MARGIN = 2, median)
nd$.lower <- apply(predictions, 2, quantile,  0.025)
nd$.upper <- apply(predictions, 2, quantile, 0.975)

nd2 <- expand.grid(delta_year_day = mean(df$delta_year_day, na.rm = T), 
                   year = unique(df$year), 
                   functional_group = unique(df$functional_group)) %>%
  mutate(year.cont = (year - mean(year)) / 10, 
         year_fac = as.character(year))

predictions <- rstanarm::posterior_epred(mod_brms, newdata = nd2, re.form = ~(year.cont|functional_group))

nd2$ytilda <- apply(predictions, MARGIN = 2, median)
nd2$.lower <- apply(predictions, 2, quantile,  0.025)
nd2$.upper <- apply(predictions, 2, quantile, 0.975)

ggplot(nd, aes(x = year, y = ytilda))+
  geom_line(aes(color = comname), show.legend = F)+
  geom_line(data = nd2, aes(x = year, y = ytilda), linewidth = 1)+
  geom_ribbon(data = nd2, aes(x = year, y = ytilda, ymin = .lower, ymax = .upper), alpha = 0.1)+
  facet_wrap(~functional_group)+
  theme_bw()+
  labs(y = "Predicted distance between seasonal centroids", x = "")
ggsave("Figures/FG_temporaltrends.png")

ggplot(nd, aes(x = year, y = ytilda))+
  geom_point(data = df, aes(x = year, y = dist_km, color = comname), size = 0.9, show.legend = F)+
  geom_line(aes(color = comname), show.legend = F)+
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.1)+
  facet_wrap(~comname)+
  theme_bw()+
  labs(y = "Predicted distance between seasonal centroids", x = "")
ggsave("Figures/Species_temporaltrends.png")



mod_brms %>%
  spread_draws(r_functional_group[c,t]) %>%
  filter(t == "year.cont") %>%
  ggplot()+
  geom_halfeyeh(aes(x = r_functional_group, y = c))

out <- mod_brms %>%
  spread_draws(`r_functional_group:comname`[string,param_type]) %>%
  separate(string, into =c("functional_group", "species"), sep = "_") %>%
  filter(param_type == "year.cont")

out %>%
  group_by(species, functional_group) %>%
  tidybayes::median_qi(`r_functional_group:comname`, .width = c(0.75, 0.95)) %>%
  ggplot(aes(x = `r_functional_group:comname`, y = forcats::fct_reorder(species, `r_functional_group:comname`)))+
  geom_pointinterval(aes(xmin = .lower, xmax = .upper))+
  geom_vline(xintercept = 0, linetype = 3, color = "gray")+
  labs(x = "Change in seasonal distance over time", y = "Species")+
  theme_classic()
ggsave("Figures/Species_coefplot.png")




mod_brms %>%
  spread_draws(r_functional_group[,])

var_names <- get_variables(mod_brms)


mod_brms %>%
  gather_draws(!!!syms(var_names[90:95]))

, b_.*[i], regex = TRUE


