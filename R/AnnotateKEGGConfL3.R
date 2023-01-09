library(pbapply)
library(tidyverse)

### Annotate Confidence Level 3: KEGG

# Prepare all data -------------------------------------------------------------------

# Theoretical
KEGG.original <- read_csv("example_data/KEGG_Data.csv") %>%
  rename(kegg_id = OtherCmpds,
         compound_name = AllNames)

# Isolate pos and neg theoretical data
KEGG.pos <- KEGG.original %>%
  select(kegg_id, compound_name, PosMZ) %>%
  mutate(z = 1) %>%
  rename(mz = PosMZ)

KEGG.neg <- KEGG.original %>%
  select(kegg_id, compound_name, NegMZ) %>%
  mutate(z = -1) %>%
  rename(mz = NegMZ)

# Combine to full theoretical dataset
KEGG.data <- KEGG.neg %>%
  rbind(KEGG.pos)

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


# Functions ---------------------------------------------------------------

CalcMzSimScore <- function(mz_exp, mz_theo, flex) {
  # Takes in experimental and theoretical m/z values in a data frame, as well as a user-defined flexibility.
  #
  # Returns a similarity score.
  similarity.score = exp(-0.5 * (((mz_exp - mz_theo) / flex) ^ 2))

  return(similarity.score)
}


CalcTotalSimScore <- function(ms1_sim) {
  # Takes in the similarity scores calculated by other functions.
  #
  # Returns: a "total" similarity score according to which scores are present.

  total.sim.score <- ms1_sim

  return(total.sim.score)
}

ConfLevel3Matches <- function(mz_i, z_i, ppm_error, theoretical_db) {
  # Pass experimental values and a theoretical data frame to this function to produce a new nested
  # column with all potential matches. Each observation in each column is an argument
  # (mz, rt, col, z, concatentenated ms2). The theoretical data frame should also be in the appropriate format.
  # Experimental values are compared to the complete set of theoretical value and matches within mz window,
  # z, and column. Calculates similarity scores for rt, mz, ms2.
  #
  # Returns: A column of nested dataframes containing all potential theoretical matches.
  potential.matches <- theoretical_db %>%
    filter(mz < mz_i + ((mz_i * ppm_error)/1e6) & mz > mz_i - ((mz_i * ppm_error)/1e6)) %>%
    filter(z == z_i) %>%
    mutate(MS1SimScore = CalcMzSimScore(mz_exp = mz_i, mz_theo = mz, flex = 5)) %>% ## TODO do we want these to not be hard coded?
    mutate(TotalSimScore = CalcTotalSimScore(MS1SimScore)) %>%
    select(compound_name, ends_with("SimScore"))

  return(potential.matches)
}


# Example ---------------------------------------------------------------
# Run the ConfLevel3Matches on a single row, which produces a data frame
# of potential matches and all similarity scores.
print(experimental.data[55, ])
single.frame <- ConfLevel3Matches(mz_i = experimental.data$mz[55], z_i = experimental.data$z[55],
                                  ppm_error = 100, theoretical_db = KEGG.data)


## Produces a dataframe with a nested column of all potential matches for each row of the experimental data.
start.time <- Sys.time()

all.matches <- experimental.data %>%
  rowwise() %>%
  mutate(matches = list(ConfLevel3Matches(mz_i = mz, z_i = z, ppm_error = 10,
                                          theoretical_db = KEGG.data))) %>%
  ungroup() %>%
  mutate(top_choice=sapply(matches, function(cmpd_matches){
    cmpd_matches %>%
      group_by(compound_name) %>%
      summarize(med_sim_overall=median(TotalSimScore, na.rm=TRUE)) %>%
      arrange(desc(med_sim_overall)) %>%
      slice(1) %>%
      pull(compound_name)
  }))

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)




