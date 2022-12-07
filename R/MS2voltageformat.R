library(tidyverse)
## Script to create correctly formatted MS2 data: concatenated with voltage. See example below.

## The original file is created by WKumler from the original MS2 data. He may adjust this to avoid the hypothetical second dda grab,
## therefore avoiding the superimposed lollipop issue currently handled by the "4:6", "9:11" grouping in the consensus script.

originalmsms <- read.csv("example_data/Ingalls_Lab_Standards_MSMS.csv")

## Change from df format to concatenated format
# Seems to work but has overlapping text in the "view" frame
voltagedmsms <- originalmsms %>%
  select(-filename) %>%
  group_by(voltage) %>%
  group_by(compound_name) %>%
  summarize(MS2 = paste(voltage, MS2, sep = "V ", collapse = ": ")) %>%
  as.data.frame()

## Change from concatenated to df
dfmsms <- voltagedmsms %>%
  separate_rows(MS2, sep = ": ") %>%
  separate(MS2, into = c("voltage", "MS2"), sep = "V ")


#write.csv(voltagedmsms, "~/work/phobos/example_data/Concatenated_Voltage_MSMS.csv", row.names = FALSE)
