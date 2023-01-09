library(pbapply)
library(tidyverse)

### Annotate Confidence Level 4: Everything else!


# Prepare data, experimental only! -------------------------------------------------------------------

# Experimental
ingalls.standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
                              stringsAsFactors = FALSE, header = TRUE) %>%
  select(compound_name = Compound_Name, HILIC_Mix, mz, rt = RT_minute, column = Column, z) %>%
  mutate(rt = rt * 60) %>%
  distinct()

experimental.data <- read_csv("example_data/Ingalls_Lab_Standards_MSMS.csv") %>%
  filter(str_detect(filename, "pos5|neg5")) %>%
  separate_rows(MS2, sep = "; ") %>%
  separate(MS2, into = c("mz", "int"), sep = ", ") %>%
  mutate(mz=as.numeric(mz)) %>%
  mutate(int=as.numeric(int)) %>%
  group_by(compound_name, filename, voltage) %>%
  mutate(int=int/max(int)*100) %>%
  arrange(desc(int)) %>%
  filter(int>1) %>%
  group_by(voltage, compound_name, filename) %>%
  summarise(MS2 = paste(mz, int, sep = ", ", collapse = "; "), .groups = "drop") %>%
  mutate(z = ifelse(str_detect(filename, "pos"), 1, -1)) %>%
  mutate(HILIC_Mix = ifelse(str_detect(filename, "Mix1"), "Mix1", "Mix2")) %>%
  rowwise() %>%
  mutate(MS2 = paste(voltage, MS2, sep = "V ", collapse = ": ")) %>%
  left_join(ingalls.standards, by = c("compound_name", "z", "HILIC_Mix")) %>%
  select(compound_name, voltage, mz, rt, column, z, MS2) %>%
  drop_na()


## original
  Confidence.Level.4 <- Confidence.Level.3 %>%
    dplyr::mutate(confidence_rank = ifelse(is.na(confidence_rank), 4, confidence_rank)) %>%
    dplyr::arrange(primary_key)



