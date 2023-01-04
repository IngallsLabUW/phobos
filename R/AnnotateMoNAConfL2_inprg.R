library(pbapply)
library(tidyverse)

### Annotate Confidence Level 2: Comparisons with the MoNA Database

# Notes -------------------------------------------------------------------
# In previous versions we have subtracted hydrogen for reference. I've included that here.
# We are ignoring voltage since MoNA doesn't always match with the voltages we use,
#   but the columns are still included in the download and associated data.
# Should we make the filters more permissive on this level?
# TODO: The MoNA spreadsheets are still in the original KRH download, located in the example_data folder.
#   The formatting code no longer works, so that part will need to be rewritten.
# TODO: The total similarity score needs to be more flexible depending on what is present/absent. Messy right now.

# Outline -------------------------------------------------------------------
# Use the downloaded MoNA relational spreadsheets as theoretical data. Downloads here: https://mona.fiehnlab.ucdavis.edu/downloads
# Use "single" MS2 data from the fifth run as experimental data.
# Find similarity and matching using only mz, MS2, and z.


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
  select(ID, compound_name = Names, voltage, mz, z, MS2) #%>%
  #####################################################################
  #filter(str_detect(compound_name, regex("ectoine", ignore_case = TRUE)))

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
  as.data.frame() #%>%
  #########################################
  #filter(str_detect(compound_name, "Ectoine"))

# Functions ---------------------------------------------------------------
mz_i <- experimental.data[55, 3]
z_i <- experimental.data[55, 4]
MS2str_i <- experimental.data[55, 5]
theoretical_db <- MoNA
ppm_error <- 10000

test.match <- ConfLevel2Matches(mz_i = mz_i, z_i = z_i, MS2str_i = MS2str_i, ppm_error = ppm_error, theoretical_db = MoNA)

ConfLevel2Matches <- function(mz_i, z_i, MS2str_i, ppm_error, theoretical_db) {
  # Pass experimental values and a theoretical data frame to this function to produce a new nested
  # column with all potential matches. Each observation in each column is an argument.
  #
  # Returns: A column of nested dataframes containing all potential theoretical matches.
  potential.matches <- theoretical_db %>%
    filter(mz < mz_i + ((mz_i * ppm_error)/1e6) & mz > mz_i - ((mz_i * ppm_error)/1e6)) %>%
    filter(z == z_i) %>%
    mutate(MS1SimScore = CalcMzSimScore(mz_exp = mz_i, mz_theo = mz, flex = 20)) %>% ## TODO: these are still hard coded, same as CL1
    mutate(MS2SimScore = as.numeric(ifelse(str_detect(MS2, ","),
                                           pblapply(MS2, CalcMS2SimScore, ms2_theo = MS2str_i, flex = 0.02), NA))) %>%
    mutate(TotalSimScore = CalcTotalSimScore(MS1SimScore, MS2SimScore)) %>%
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

CalcMzSimScore <- function(mz_exp, mz_theo, flex) {
  # Takes in experimental and theoretical m/z values in a data frame, as well as a user-defined flexibility.
  #
  # Returns a similarity score.
  similarity.score = exp(-0.5 * (((mz_exp - mz_theo) / flex) ^ 2))

  return(similarity.score)
}

CalcMS2SimScore <- function(ms2_exp, ms2_theo, flex) {
  # Takes in experimental and theoretical MS2 values in a data frame, as well as a user-defined flexibility.
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

CalcTotalSimScore <- function(ms1_sim, ms2_sim) { ## TODO: This needs to be more flexible, depending on which similarity scores are present
  # Takes in the similarity scores calculated by other functions.
  #
  # Returns: a "total" similarity score according to which scores are present.

  total.sim.score <- ifelse(is.na(ms2_sim), ms1_sim, ((ms1_sim + ms2_sim) / 2) * 100)

  return(total.sim.score)
}


## Example: Produces a data frame of potential matches and all similiarity scores for a single row of experimental data.
print(experimental.data[55, ])

single.frame <- ConfLevel2Matches(mz_i = experimental.data$mz[55], z_i = experimental.data$z[55],
                                  MS2str_i = experimental.data$MS2[55],
                                  ppm_error = 100, theoretical_db = MoNA)


## Produces a dataframe with a nested column of all potential matches for each row of the experimental data.
start.time <- Sys.time()

all.matches <- experimental.data %>%
  rowwise() %>%
  mutate(matches = list(ConfLevel2Matches(mz_i = mz, z_i = z, MS2str_i = MS2,
                                          ppm_error = 10, theoretical_db = MoNA))) %>%
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


## None of them actually match...
wrong.matches <- all.matches %>%
  filter(compound_name != top_choice)

