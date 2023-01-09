library(pbapply)
library(tidyverse)

### Annotate Confidence Level 3: MoNA

# Outline -------------------------------------------------------------------
# Use the downloaded MoNA relational spreadsheets as theoretical data. Downloads here: https://mona.fiehnlab.ucdavis.edu/downloads
# Use "single" MS2 data from the fifth run as experimental data.
# Find similarity and matching using only mz and z.

# Prepare all data -------------------------------------------------------------------

## Theoretical data: subtract hydrogen and mutate z columns
MoNA.Neg <- read.csv("example_data/MoNA_RelationalSpreadsheets/NEG_Spectra.csv") %>%
  mutate(mz = M_mass + 1.0072766) %>%
  mutate(z = -1)
MoNA.Pos <- read.csv("example_data/MoNA_RelationalSpreadsheets/POS_Spectra.csv") %>%
  mutate(z = 1) %>%
  mutate(mz = M_mass - 1.0072766)

MoNA <- MoNA.Neg %>%
  rbind(MoNA.Pos) %>%
  mutate(voltage = trimws(str_remove(CE, "[^ ]*$"))) %>% # icky but works
  mutate(MS2 = paste(voltage, "V ", spectrum_KRHform_filtered, sep = "")) %>%
  select(ID, compound_name = Names, voltage, mz, z, MS2)

## Experimental data: prepare the same way as in confidence level 1
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
  summarise(MS2 = paste(mz, int, sep = ", ", collapse = "; ")) %>%
  mutate(z = ifelse(str_detect(filename, "pos"), 1, -1)) %>%
  mutate(HILIC_Mix = ifelse(str_detect(filename, "Mix1"), "Mix1", "Mix2")) %>%
  rowwise() %>%
  mutate(MS2 = paste(voltage, MS2, sep = "V ", collapse = ": ")) %>%
  left_join(ingalls.standards, by = c("compound_name", "z", "HILIC_Mix")) %>%
  select(compound_name, voltage, mz, z, MS2) %>%
  drop_na() %>%
  as.data.frame()

# Functions ---------------------------------------------------------------
ConfLevel3Matches <- function(mz_i, z_i, ppm_error, theoretical_db) {
  # Pass experimental values and a theoretical data frame to this function to produce a new nested
  # column with all potential matches. Each observation in each column is an argument.
  #
  # Returns: A column of nested dataframes containing all potential theoretical matches.
  potential.matches <- theoretical_db %>%
    filter(mz < mz_i + ((mz_i * ppm_error)/1e6) & mz > mz_i - ((mz_i * ppm_error)/1e6)) %>%
    filter(z == z_i) %>%
    mutate(MS1SimScore = CalcMzSimScore(mz_exp = mz_i, mz_theo = mz, flex = 20)) %>% ## TODO: these are still hard coded, same as CL1
    mutate(TotalSimScore = CalcTotalSimScore(MS1SimScore)) %>%
    select(compound_name, voltage, ends_with("SimScore"))

  return(potential.matches)
}

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


## Example: Produces a data frame of potential matches and all similiarity scores for a single row of experimental data.
print(experimental.data[55, ])

single.frame <- ConfLevel3Matches(mz_i = experimental.data$mz[55], z_i = experimental.data$z[55],
                                  ppm_error = 1000900, theoretical_db = MoNA)


## Produces a dataframe with a nested column of all potential matches for each row of the experimental data.
start.time <- Sys.time()

all.matches <- experimental.data %>%
  rowwise() %>%
  mutate(matches = list(ConfLevel3Matches(mz_i = mz, z_i = z,
                                          ppm_error = 10000, theoretical_db = MoNA))) %>%
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

