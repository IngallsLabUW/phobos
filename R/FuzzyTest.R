library(tidyverse)
library(fuzzyjoin)

## Flexible fuzzy join test script

## Ranges:
# 0 - 100 Da ~ 0.00025
# 100 - 200 Da ~ 0.0005
# 200 - 500 Da ~ 0.0012
# 500 - 1000 Da ~ 0.0025

tiny <- c(0:100)
small <- c(101:200)
medium <- c(201:500)
large <- c(501:5000)

tiny.flex <- 0.00025
small.flex <- 0.0005
medium.flex <- 0.0012
large.flex <- 0.0025

###

t <- data.frame(MassFeature = c("TinyMolecule", "LargeMolecule"),
                mz = c(89, 738580),
                molclass = c("tiny", "large"))

experimental <- read.csv("example_data/Example_Experimental_Data.csv") %>%
  select(MassFeature, mz) %>%
  mutate(molclass = case_when(
    between(mz, 0, 100) ~ "tiny",
    between(mz, 101, 200) ~ "small",
    between(mz, 201, 500) ~ "medium",
    between(mz, 501, 5000) ~ "large",
    TRUE ~ NA_character_
  )) %>%
  rbind(t)

theoretical <- read.csv("example_data/Example_Theoretical_Data.csv") %>%
  select(compound, mz) %>%
  unique()

## Currently: a join of 0.02 flexibility regardless of experimental molecule class
Current <- experimental %>%
  fuzzyjoin::difference_left_join(theoretical, by = c("mz"), max_dist = 0.02, distance_col = "distance")

## We want: a flexible, moleculeclass-based join
# Try to do fuzzy join based on molecule class, if a molecule is "small" it gets fuzzyjoined by 0.02, if it is "large" it gets fuzzyjoined by 0.05 eg

# This...
ci_str_detect <- function(x, y){str_detect(x, regex(y, ignore_case = TRUE))}
fuzzy_left_join(df1, df2, match_fun = ci_str_detect, by = c(col1 = "col4"))

# Or something like this...
fuzzy_join(experimental, theoretical,
           by = c("mz"),
           match_fun = list(function(x, y) ), ## some kind of function here, "if mz is in a certain range, join on this flex, else another flex"
           mode = "left")

## Unsuccessful loop attempt

if("small" %in% experimental$molclass) {
  testjoin1 <- experimental %>%
    fuzzyjoin::difference_left_join(theoretical, by = c("mz"), max_dist = 0.02, distance_col = "distance")
} else if("medium" %in% experimental$molclass) {
  testjoin2 <- experimental %>%
    fuzzyjoin::difference_left_join(theoretical, by = c("mz"), max_dist = 0.05, distance_col = "distance")
}

## Row by row attempt
df <- data.frame(MassFeature=character(),
                 mz=numeric(),
                 molclass=character(),
                 stringsAsFactors=FALSE)

for(i in experimental$mz) {

  if (between(i, 0, 100) == TRUE) {
    print("tiny")
    df <- experimental %>%
      fuzzyjoin::difference_left_join(theoretical, by = c("mz"), max_dist = small.flex, distance_col = "distance")

  } else if (between(i, 101, 200) == TRUE) {
    print("small")

  } else if (between(i, 201, 500) == TRUE) {
    output <- "medium"

  } else if (between(i, 501, Inf) == TRUE) {
    output <- "large"
  }
}






