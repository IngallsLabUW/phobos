library(tidyverse)
library(fuzzyjoin)

## Flexible fuzzy join test script

## Suggesetion:
# 0 - 150 ~ 5ppm
# 150 - 200 ~ 10 ppm

small.molecule.range <- c(0:100)



flexibility <- 100

Unknown <- data.frame(molecule = c("alpha", "beta", "gamma"),
                      mass = c(100, 200, 300))

Known <- data.frame(molecule = c("omega", "tau", "sigma", "rho"),
                    mass = (sample(100:300, 4)))

Current <- Unknown %>%
  fuzzyjoin::difference_left_join(Known, by = c("mass"), max_dist = flexibility)
