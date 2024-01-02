#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------------Last updated 2022-06-25---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

library(tidyverse)
library(amt)
library(terra)

# Code for simulations (it is still in the development version of amt)
# source("https://raw.githubusercontent.com/jmsigner/amt/simulate/R/simulate.R")

set.seed(1323)

# We will use simulated data again for gps data
gps <- readRDS("Data/gps.rds")

# and the environment.
hab <- readRDS("Data/hab.rds")
hab <- c(rast(hab$forage), 
         rast(hab$temp),
         rast(hab$pred),
         rast(hab$cover),
         rast(hab$dist_to_water),
         rast(hab$dist_to_cent),
         rast(hab$rand))

# Next we need to prepare the data. 
stps <- gps %>% 
  make_track(x, y, t, crs = 32612) %>% 
  steps() %>% 
  random_steps(n_control = 20) %>% 
  extract_covariates(hab, where = "both")

# And fit a first model
m1 <- fit_issf(stps, 
               # Response
               case_ ~ 
                 # Habitat
                 forage_end + pred_end +
                 # Movement
                 sl_ + log(sl_) + cos(ta_) + 
                 # Stratum
                 strata(step_id_),
               # Need this later for model predictions
               model = TRUE)

# We can now use this model `m1` to create a redistribution kernel. This
# requires at least three arguments: 
# - x: a fitted iSSF model.
# - start: the start position (x, y) and the direction. 
# - map: the environmental covariates
k1 <- redistribution_kernel(x = m1, start = make_start(stps[2, ]), map = hab, as.rast=T)
plot(k1$redistribution.kernel)

# Finally, `k1` can be used to simulate a path. 
# This take ~ 30 seconds
p1 <- simulate_path(k1, n = 500)

# Plot results
plot(hab[["forage"]]); lines(p1)

plot(hab[["pred"]]); lines(p1)

# What did we observe
plot(hab[["forage"]]); lines(gps$x, gps$y); lines(p1, col = "red")

# Now lets adjust the model by including home-ranging behavior. 
m2 <- fit_issf(stps, 
               # Response
               case_ ~ 
                 # Habitat
                 forage_end + pred_end +
                 # Movement
                 sl_ + log(sl_) + cos(ta_) + 
                 # Home ranging
                 x2_ + y2_ + I(x2_^2 + y2_^2) +
                 # Stratum
                 strata(step_id_),
               # Need this later for model predictions
               model = TRUE)

k2 <- redistribution_kernel(x = m2, start = make_start(stps[2, ]), map = hab, as.rast=T)

# Simulate again
# Takes about 30 seconds. 
p2 <- simulate_path(k2, n = 500)

plot(hab[["forage"]]); lines(p1); lines(p2, col = "blue")
plot(hab[["pred"]]); lines(p1); lines(p2, col = "blue")

# What did we observe
plot(hab[["forage"]]); lines(gps$x, gps$y); lines(p1, col = "red"); lines(p2, col = "blue")

