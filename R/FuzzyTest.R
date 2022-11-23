library(tidyverse)
library(fuzzyjoin)

## Flexible fuzzy join test script

## Functions
GetFlex <- function(mz) { # Should molecule class and flex windows be arguments? Are they gonna change?
  molecule.class <- c(100, 200, 500)
  flex <- c(0.00025, 0.0005, 0.0012)
  for (i in 1:length(molecule.class)) {
    if (mz <= molecule.class[i]) {
      return(c(flex[i]))
    }
  }
  return(c(0.0025))
}

MzCloseEnough <- function(x, y) {
  abs(x[, 1] - y) <= x[, 2]
}

## Load data
experimental <- read.csv("example_data/Example_Experimental_Data.csv") #%>%
  # select(MassFeature, mz)

theoretical <- read.csv("example_data/Example_Theoretical_Data.csv") %>% ## The theoretical data must be only two unique columns for this to work
  select(compound, mz) %>%
  unique() %>%
  drop_na()

## Sort experimental mzs into flex category
experimental <- experimental %>%
  rowwise() %>%
  mutate(flex = GetFlex(mz))

## Fuzzy join experimental and theoretical data according to flex category
flex.join <- experimental %>%
  fuzzyjoin::fuzzy_join(theoretical, multi_by = list(x=c("mz", "flex"), y="mz"),
                        multi_match_fun = MzCloseEnough, mode = "left")

