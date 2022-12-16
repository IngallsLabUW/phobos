library(pbapply)
library(tidyverse)

### WKumler + RML Startup script

# Notes -------------------------------------------------------------------
# Prepping the data can be kind of a big step, see preparation for single.experimental
# How do we want to handle voltage columns? Should we match on them?
# MS2 Sim Score separate_rows produced some headaches.
# Need to make an adjustment for "n consensus files", not 4:6/9:11 for the testing purposes. This could be a little softer...
# increasing suspicion as we move away from 5. How many "intensity clusters" at mz clusters do we have? Now we're accounting
# for two intensity clusters, but it's possible we'd have three, four...
# Full run is taking 51 minutes.
# Without the MS2 steps, the total process took 50 seconds.
# ...and then we dropped the "remove any single MS2 files beforehand step" which bumped it back up to 7 minutes...
# Amino Propanesulfonic acid: any relation to 3-Aminopropanesulfonate?

# Outline -------------------------------------------------------------------
# Create consensus MS2 data from four of the five MSMS runs of the Ingalls standards.
# Using the functions below, create a Total Similarity Score for entries that fall within
# column, z, and mz filters.
# The cosine similarity will be calculated between "yellow triangle/consensus msms spectra" and experimental dot

# Prepare all data -------------------------------------------------------------------
ingalls.standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
                                                    stringsAsFactors = FALSE, header = TRUE) %>%
  select(compound_name = Compound_Name, HILIC_Mix, mz, rt = RT_minute, column = Column, z) %>%
  mutate(rt = rt * 60) %>%
  distinct()

## Theoretical and Experimental data, in their standardized concatenated voltage input format
consensus.theoretical <- read_csv("example_data/Consensed_Theoretical_FourRuns.csv") %>%
  mutate(z = ifelse(polarity == "pos", 1, -1)) %>%
  rename(MS2 = consensus_MS2) %>%
  left_join(ingalls.standards, by = c("compound_name", "z")) %>%
  separate_rows(MS2, sep = ": ") %>%
  mutate(voltage = sub("\\V.*", "", MS2)) %>% ## Currently just creating this voltage column
  select(compound_name, voltage, mz, rt, column, z, MS2)

single.experimental <- read_csv("example_data/Ingalls_Lab_Standards_MSMS.csv") %>%
  filter(str_detect(filename, "pos5|neg5")) %>%
  separate_rows(MS2, sep = "; ") %>%
  separate(MS2, into = c("mz", "int"), sep = ", ") %>%
  mutate(mz=as.numeric(mz)) %>%
  mutate(int=as.numeric(int)) %>%
  group_by(compound_name, filename, voltage) %>%
  mutate(int=int/max(int)*100) %>%
  ungroup() %>%
  arrange(desc(int)) %>%
  filter(int>1) %>%
  group_by(voltage, compound_name, filename) %>% #, z) %>%
  summarise(MS2 = paste(mz, int, sep = ", ", collapse = "; ")) %>%
  mutate(z = ifelse(str_detect(filename, "pos"), 1, -1)) %>%
  mutate(HILIC_Mix = ifelse(str_detect(filename, "Mix1"), "Mix1", "Mix2")) %>%
  rowwise() %>%
  mutate(MS2 = paste(voltage, MS2, sep = "V ", collapse = ": ")) %>%
  as.data.frame() %>%
  select(-filename) %>%
  left_join(ingalls.standards, by = c("compound_name", "z", "HILIC_Mix")) %>%
  select(compound_name, voltage, mz, rt, column, z, MS2) %>%
  drop_na()

# Functions ---------------------------------------------------------------
# mz_i <- single.experimental$mz[25]
# rt_i <- single.experimental$rt[25]
# col_i <- single.experimental$column[25]
# z_i <- single.experimental$z[25]
# MS2str_i <- single.experimental$MS2[25]
# ppm_error <- 10
# theoretical_db <- consensus.theoretical

CreatePotentialMatches_1 <- function(mz_i, rt_i, col_i, z_i,
                                     MS2str_i,
                                     ppm_error, theoretical_db) {
  # Pass experimental values and a theoretical data frame to this function to produce a new nested
  # column with all potential matches. Each observation in each column is an argument
  # (mz, rt, col, z, ms2 in concatenated format). The theoretical data frame should be in the appropriate format.
  # Experimental values are compared to the complete set of theoretical value and matches within mz window,
  # z, and column. Calculates similarity scores for rt, mz, ms2.
  #
  # Returns: A column of nested dataframes containing all potential theoretical matches.
  potential.matches <- theoretical_db %>%
    filter(mz < mz_i + ((mz_i * ppm_error)/1e6) & mz > mz_i - ((mz_i * ppm_error)/1e6)) %>%
    filter(column == col_i) %>%
    filter(z == z_i) %>%
    mutate(MS1SimScore = CalculateMzSimScore_1(mz_exp = mz_i, mz_theo = mz, flex = 5)) %>%
    mutate(RT1SimScore = CalculateRTSimScore_1(rt_exp = rt_i, rt_theo = rt, flex = 30)) %>%
    mutate(MS2SimScore = as.numeric(ifelse(str_detect(MS2, ","),
                               pblapply(MS2, CalculateMS2SimScore_1, ms2_theo = MS2str_i, flex = 0.02), NA))) %>%
    mutate(TotalSimScore = CalculateTotalSimScore_1(MS1SimScore, RT1SimScore, MS2SimScore)) %>%
    select(compound_name, voltage, ends_with("SimScore")) ## Should this include voltage?

  return(potential.matches)
}

MakeScantable <- function(concatenated.scan) {
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

CalculateMzSimScore_1 <- function(mz_exp, mz_theo, flex) {
  similarity.score = exp(-0.5 * (((mz_exp - mz_theo) / flex) ^ 2))

  return(similarity.score)
} ## TODO The flex values are hardcoded in the Calctotalsimscore1 function

# ms2_exp <- single.experimental$MS2[25]
# ms2_theo <- consensus.theoretical$MS2[29]

CalculateMS2SimScore_1 <- function(ms2_exp, ms2_theo, flex) {
# Comparison will be between experimental and "consensus" spectra according to Horai et. al 2010
# as justification for using spectra
  scan1 <- MakeScantable(ms2_exp)
  scan2 <- MakeScantable(ms2_theo)

  if ((is.na(scan1) || is.na(scan2) || nrow(scan1) < 2 || nrow(scan2) < 2)) {
  #if (is.na(scan1) || is.na(scan2) || all(dim(scan1) == c(1,1)) || all(dim(scan2) == c(1,1))) {

    return(NA)

  } else {
    # Need consensus to assign weights
    weight1 <- (scan1[, "mz"] ^ 2) * (scan1[, "intensity"] ^ 0.5)
    weight2 <- (scan2[, "mz"] ^ 2) * (scan2[, "intensity"] ^ 0.5)

    diff.matrix <- sapply(scan1[, "mz"], function(x) scan2[, "mz"] - x)
    same.index <- which(abs(diff.matrix) < flex, arr.ind = TRUE)
    cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
      (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))

    return(cosine.similarity)
  }
}

CalculateRTSimScore_1 <- function(rt_exp, rt_theo, flex) {
  similarity.score = exp(-0.5 * (((rt_exp - rt_theo) / flex) ^ 2))

  return(similarity.score)
}

CalculateTotalSimScore_1 <- function(ms1_sim, rt_sim, ms2_sim) {
  total.similarity.score <- ((ms1_sim + rt_sim + ms2_sim) / 3) * 100

  return(total.similarity.score)
}


## Example: Produces dataframe of potential matches and all sim scores for a single row of experimental data.
single.frame <- CreatePotentialMatches_1(mz_i = single.experimental$mz[25], rt_i = single.experimental$rt[25],
                                col_i = single.experimental$column[25], z_i = single.experimental$z[25],
                                MS2str_i = single.experimental$MS2[25],
                                ppm_error = 100,
                                theoretical_db = consensus.theoretical)


## Produces a dataframe with a nested column
start.time <- Sys.time()
AllOutput <- single.experimental %>%
  slice(1:30) %>%
  rowwise() %>%
  mutate(matches = list(CreatePotentialMatches_1(mz_i = mz, rt_i = rt, col_i = column, z_i = z,
                                                 MS2str_i = MS2,
                                                 ppm_error = 10,
                                                 theoretical_db = consensus.theoretical))) %>%
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


### Fix the weird matching
testing <- AllOutput %>%
  filter(compound_name != top_choice)

