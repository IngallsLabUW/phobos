library(pbapply)
library(tidyverse)

### Annotate Confidence Level 1: Ingalls Standards comparison

# Notes -------------------------------------------------------------------
# Prepping the data can be kind of a big step, see section for experimental.data.
# Need to make an adjustment for "n consensus files", maybe not the current 4:6/9:11 being used for the testing.
#   How many "intensity clusters" at mz points do we have? Now we're accounting
#   for two intensity clusters, but it's possible we'd have three, four...
# TODO: Randomize the 4:1 experimental:theoretical choices.
# TODO: Need consensus to assign weights in the MS2 similarity score, according to Horai et al. 2010,
#   also I believe we were going to use this as document/justification for spectra consensus.
# TODO: The flex arguments in the similarity score calculations are hard coded within the
#   big ConfLevel1 function. Should we change that?
# TODO: I have a progress bar wrapped around the MS2 Sim Score function but this could be better.
# TODO: The med_sim_overall function is choosing some incorrect top choices. Going with the max
#   sim score has fewer wrong choices, but still a decent amount. This might change with adjusted weights?

# Outline -------------------------------------------------------------------
# Use "consensed" MS2 data created from four of the five Ingalls Standards run as theoretical data.
#   See TakeMS2Consensus.R script for details.
# Use "single" MS2 data from the fifth run as experimental data.
# Create a Total Similarity Score for entries that fall within column, z, and mz filters.
# The cosine similarity will be calculated between "yellow triangle/consensus msms spectra" and experimental dot.
#   See MS2ClusterGraph.R for details.
# The below script produces a data frame with a nested column of data frames containing all potential matches for each row.
# Whole process takes about 1.5 minutes.

# Prepare all data -------------------------------------------------------------------
ingalls.standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
                                                    stringsAsFactors = FALSE, header = TRUE) %>%
  select(compound_name = Compound_Name, HILIC_Mix, mz, rt = RT_minute, column = Column, z) %>%
  mutate(rt = rt * 60) %>%
  distinct()

# Theoretical and Experimental data, in their standardized concatenated voltage input format.
# Theoretical data is the "consensed" first four runs, and experimental is the final fifth run.
theoretical.data <- read_csv("example_data/Consensed_Theoretical_FourRuns.csv") %>%
  mutate(z = ifelse(polarity == "pos", 1, -1)) %>%
  left_join(ingalls.standards, by = c("compound_name", "z")) %>%
  separate_rows(consensus_MS2, sep = ": ") %>%
  mutate(voltage = sub("\\V.*", "", consensus_MS2)) %>%
  select(compound_name, voltage, mz, rt, column, z, MS2 = consensus_MS2)

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
  select(compound_name, voltage, mz, rt, column, z, MS2) %>%
  drop_na()

# Functions ---------------------------------------------------------------
CalcMS2SimScore <- function(ms2_exp, ms2_theo, flex) {
  # Takes in concatenated experimental and theoretical MS2 values in a data frame,
  # as well as a user-defined flexibility.
  #
  # Returns a cosine similarity score.
  scan1 <- MakeScantable(ms2_exp)
  scan2 <- MakeScantable(ms2_theo)

  scans <- c(scan1, scan2)
  length.check <- lapply(scans, length) < 2

  if (TRUE %in% is.na(scans) || TRUE %in% length.check) {

    return(NA)

  } else {
    weight1 <- (scan1[, "mz"] ^ 2) * (scan1[, "intensity"] ^ 0.5)
    weight2 <- (scan2[, "mz"] ^ 2) * (scan2[, "intensity"] ^ 0.5)

    diff.matrix <- sapply(scan1[, "mz"], function(x) scan2[, "mz"] - x)
    same.index <- which(abs(diff.matrix) < flex, arr.ind = TRUE)
    cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
      (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))

    return(cosine.similarity)
  }
}

CalcMzSimScore <- function(mz_exp, mz_theo, flex) {
  # Takes in experimental and theoretical m/z values in a data frame, as well as a user-defined flexibility.
  #
  # Returns a similarity score.
  similarity.score = exp(-0.5 * (((mz_exp - mz_theo) / flex) ^ 2))

  return(similarity.score)
}

CalcRTSimScore <- function(rt_exp, rt_theo, flex) {
  # Takes in experimental and theoretical retention time values in a data frame, as well as a user-defined flexibility.
  #
  # Returns a similarity score.
  similarity.score = exp(-0.5 * (((rt_exp - rt_theo) / flex) ^ 2))

  return(similarity.score)
}

CalcTotalSimScore <- function(ms1_sim, rt_sim, ms2_sim) {
  # Takes in the similarity scores calculated by other functions.
  #
  # Returns: a "total" similarity score according to which scores are present.

  total.sim.score <- ifelse(is.na(ms2_sim), ((ms1_sim + rt_sim) / 2) * 100, ((ms1_sim + rt_sim + ms2_sim) / 3) * 100)

  return(total.sim.score)
}

ConfLevel1Matches <- function(mz_i, rt_i, col_i, z_i, MS2str_i, ppm_error, theoretical_db) {
  # Pass experimental values and a theoretical data frame to this function to produce a new nested
  # column with all potential matches. Each observation in each column is an argument
  # (mz, rt, col, z, concatentenated ms2). The theoretical data frame should also be in the appropriate format.
  # Experimental values are compared to the complete set of theoretical value and matches within mz window,
  # z, and column. Calculates similarity scores for rt, mz, ms2.
  #
  # Returns: A column of nested dataframes containing all potential theoretical matches.
  potential.matches <- theoretical_db %>%
    filter(mz < mz_i + ((mz_i * ppm_error)/1e6) & mz > mz_i - ((mz_i * ppm_error)/1e6)) %>%
    filter(column == col_i) %>%
    filter(z == z_i) %>%
    mutate(MS1SimScore = CalcMzSimScore(mz_exp = mz_i, mz_theo = mz, flex = 5)) %>% ## TODO do we want these to not be hard coded?
    mutate(RT1SimScore = CalcRTSimScore(rt_exp = rt_i, rt_theo = rt, flex = 30)) %>%
    mutate(MS2SimScore = as.numeric(ifelse(str_detect(MS2, ","),
                                           pblapply(MS2, CalcMS2SimScore, ms2_theo = MS2str_i, flex = 0.02), NA))) %>%
    mutate(TotalSimScore = CalcTotalSimScore(MS1SimScore, RT1SimScore, MS2SimScore)) %>%
    select(compound_name, voltage, ends_with("SimScore"))

  return(potential.matches)
}

MakeScantable <- function(concatenated.scan) {
  # Takes in single character string of MS2, split by voltage., and creates
  # Required format for concatenated.scan: "20V 146.11766, 100; 87.04468, 47.4; 60.08156, 26.3"
  # Format is "voltage V mass, intensity; mass, intensity"
  #
  # Returns: a mini MS2 dataframe, scaled to a max intensity of 100 and filtered to a minimum intensity of 0.5.
  requireNamespace("dplyr", quietly = TRUE)
  concatenated.scan <- trimws(gsub(".*V", "", concatenated.scan))

  if (concatenated.scan == "") {
    return(NA)

  } else {
    scantable <- read.table(text = as.character(concatenated.scan),
                            col.names = c("mz", "intensity"), fill = TRUE) %>%
      dplyr::mutate(mz = as.numeric(mz %>% stringr::str_replace(",", "")),
                    intensity = as.numeric(intensity %>% stringr::str_replace(";", "")),
                    intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
      dplyr::filter(intensity > 0.5) %>%
      dplyr::arrange(desc(intensity))

    return(scantable)
  }
}

# Example ---------------------------------------------------------------
# Run the ConfLevel1Matches on a single row, which produces a data frame
# of potential matches and all similarity scores.
print(experimental.data[71, ])
single.frame <- ConfLevel1Matches(mz_i = experimental.data$mz[71], rt_i = experimental.data$rt[71],
                                  col_i = experimental.data$column[71], z_i = experimental.data$z[71],
                                  MS2str_i = experimental.data$MS2[71], ppm_error = 100, theoretical_db = theoretical.data)


## Produces a dataframe with a nested column of all potential matches for each row of the experimental data.
start.time <- Sys.time()

all.matches <- experimental.data %>%
  rowwise() %>%
  mutate(matches = list(ConfLevel1Matches(mz_i = mz, rt_i = rt, col_i = column, z_i = z,
                                                 MS2str_i = MS2, ppm_error = 10,
                                                 theoretical_db = theoretical.data))) %>%
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


## Issues are arising from the med sim choice function: not always picking the right one.
wrong.top.match <- all.matches %>%
  filter(compound_name != top_choice)

