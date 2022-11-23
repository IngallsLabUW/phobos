library(tidyverse)
library(fuzzyjoin)

## Flexible fuzzy join test script

experimental <- read.csv("example_data/Example_Experimental_Data.csv") %>%
  select(MassFeature, mz)

theoretical <- read.csv("example_data/Example_Theoretical_Data.csv") %>%
  select(compound, mz) %>%
  unique() %>%
  drop_na()

get_flex <- function(mz) {
  intervals <- c(100, 200, 500)
  flex <- c(0.00025, 0.0005, 0.0012)
  for (i in 1:length(intervals)) {
    if (mz <= intervals[i]) {
      return(c(flex[i]))
    }
  }
  return(c(0.0025))
}
experimental <- experimental %>% rowwise() %>% mutate(flex = get_flex(mz))

mz_close_enough <- function(x, y) {
  abs(x[,1] - y) <= x[,2]
}
Current <- experimental %>%
  fuzzyjoin::fuzzy_join(theoretical, multi_by = list(x=c("mz", "flex"), y="mz"),
                        multi_match_fun = mz_close_enough, mode = 'left')

