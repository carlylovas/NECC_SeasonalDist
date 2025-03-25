library(brms)
library(tidyverse)
# library(rstanarm)
library(glmmTMB)
library(tidybayes)
# library(rstanarm)


df <- readRDS("Data/mod_Bart.rds") %>%
  rename(functional_group = SOE.24) %>%
  ungroup() %>%
  mutate(year.cont = (year - mean(year)) / 10, # This converts year (~1970-1923) into a scaled variable where a 1 unit change in the variable represents the shift over 1 decade. 
         delta_year_day = scale(delta_year_day)) # This just scales the delta_year_day variable to that everything is on similar scales.

summary(df)


# Fit the model using brms (takes about an hour on my computer)
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

write_rds(mod_brms, file = "Data/mod_brms1.rds", compress = "gz") # Write out the model file to disk. 
# mod_brms <- readRDS("Data/mod_brms1.rds") # Read it in here if you don't want to run the model

prior_summary(mod_brms) # This is the summary of the priors for the model.

mod_brms # There are still a few divergent transitions that I wasn't able to get rid of. Its possible that adjusting the priors might help, but that seems like quite a bit of work. I would advocate for running with this model for now.

pp_check(mod_brms, ndraws = 100) # seems reasonable to me! 
shinystan::launch_shinystan(mod_brms)# This will launch a GUI to interrogate the model.

#--------------------------------------------------------
## Generate predictions (e.g. figures) from the model
#---------------------------------------------------------

formerge <- df %>% 
  distinct(comname, functional_group) # Housekeeping to generate the new data to predict on.

nd <- expand.grid(delta_year_day = mean(df$delta_year_day, na.rm = T), 
                  year = unique(df$year), 
                  # functional_group = unique(df$functional_group), 
                  comname = unique(df$comname)) %>%
  mutate(year.cont = (year - mean(year)) / 10, 
         year_fac = as.character(year)) %>%
  left_join(formerge) # all this is just generating a data frame to predict on. 

predictions <- rstanarm::posterior_epred(mod_brms, newdata = nd, re.form = NULL) # Here we are using the model to generate predictions. Specifically, this will generate LINEAR predictions and will have smaller variance that posterior_predict(), because the variance is only based on the uncertainty in the expected value of the posterior predictive distribution. The residual error is ignored. For true posterior PREDICTIVE intervals we could use posterior_predict(). re.form = NULL ensures that the model accounts for ALL levels of uncertainly, e.g. uncertainly due to both fixed and random effects.

predictions[1:10,1:10] # gut check on structure

nd$ytilda <- apply(predictions, MARGIN = 2, median) # This code just extracts the median and 95% CI's for each value of year in the new data.
nd$.lower <- apply(predictions, 2, quantile,  0.025)
nd$.upper <- apply(predictions, 2, quantile, 0.975)

# The first predictive data frame was based at the species level. Here we build on that only focused on the higher level of the heirarchy -- the functional group level.
nd2 <- expand.grid(delta_year_day = mean(df$delta_year_day, na.rm = T), 
                   year = unique(df$year), 
                   functional_group = unique(df$functional_group)) %>%
  mutate(year.cont = (year - mean(year)) / 10, 
         year_fac = as.character(year))

predictions <- rstanarm::posterior_epred(mod_brms, newdata = nd2, re.form = ~(year.cont|functional_group)) # Here we generate the prediction BUT we specifically estimate them based on the random slope and intercept of functional group NOT the species level variances.

nd2$ytilda <- apply(predictions, MARGIN = 2, median) # Same housekeeping to get medians and CI's.
nd2$.lower <- apply(predictions, 2, quantile,  0.025)
nd2$.upper <- apply(predictions, 2, quantile, 0.975)


# This is a figure paneled by functional group w/ only the CI's for functional group.
ggplot(nd, aes(x = year, y = ytilda))+
  geom_line(aes(color = comname), show.legend = F)+
  geom_line(data = nd2, aes(x = year, y = ytilda), linewidth = 1)+
  geom_ribbon(data = nd2, aes(x = year, y = ytilda, ymin = .lower, ymax = .upper), alpha = 0.1)+
  facet_wrap(~functional_group)+
  theme_bw()+
  labs(y = "Predicted distance between seasonal centroids", x = "")
ggsave("Figures/FG_temporaltrends.png")

# This one I imagine going into the supplement, but it the same figure but faceted by species and includes the CI's for each species specific trend.
ggplot(nd, aes(x = year, y = ytilda))+
  geom_point(data = df, aes(x = year, y = dist_km, color = comname), size = 0.9, show.legend = F)+
  geom_line(aes(color = comname), show.legend = F)+
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.1)+
  facet_wrap(~comname)+
  theme_bw()+
  labs(y = "Predicted distance between seasonal centroids", x = "")
ggsave("Figures/Species_temporaltrends.png")


# Coef plot: here, I do a little wrangling with the mod_brms posteriors to get at the coefficients (e.g the slope estimate [beta = year.cont*comname]) for the change in seasonal distance OVER TIME for each species nested within functional group.

out <- mod_brms %>%
  spread_draws(`r_functional_group:comname`[string,param_type]) %>% # Funky data wrangling. Unfortunately, this leaves a period in the common names which might curse us later, but going to avoid dealing with it for now.
  separate(string, into =c("functional_group", "species"), sep = "_") %>%
  filter(param_type == "year.cont")

# Build out the coeffiencient plot.
out %>%
  group_by(species, functional_group) %>%
  tidybayes::median_qi(`r_functional_group:comname`, .width = c(0.75, 0.95)) %>%
  ggplot(aes(x = `r_functional_group:comname`, y = forcats::fct_reorder(species, `r_functional_group:comname`)))+
  geom_pointinterval(aes(xmin = .lower, xmax = .upper))+
  geom_vline(xintercept = 0, linetype = 3, color = "gray")+
  labs(x = "Change in seasonal distance over time", y = "Species")+
  theme_classic()
ggsave("Figures/Species_coefplot.png")