library(tidyverse)

### WKumler + RML Startup script

# Notes -------------------------------------------------------------------
# MS2 Sim Score separate_rows produced some headaches.
# Need to make an adjustment for "n consensus files", not 4:6/9:11 for the testing purposes. This could be a little softer...
## increasing suspicion as we move away from 5. How many "intensity clusters" at mz clusters do we have? Now we're accounting
## for two intensity clusters, but it's possible we'd have three, four...

# Outline -------------------------------------------------------------------
# Create consensus MS2 data from four of the five MSMS runs of the Ingalls standards.
# Using the functions below, create a Total Similarity Score for entries that fall within
# column, z, and mz filters.
# The cosine similarity will be calculated between "yellow triangle/consensus msms spectra" and experimental dot

# Prepare experimental and theoretical data -------------------------------------------------------------------
ingalls.standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
                                                    stringsAsFactors = FALSE, header = TRUE) %>%
  select(compound_name = Compound_Name, mz, rt = RT_minute, column = Column, z) %>%
  mutate(rt = rt * 60) %>%
  distinct()

## Theoretical and Experimental data, in their concatenated voltage input format
four.runs.theoretical <- read_csv("example_data/Ingalls_Lab_Standards_MSMS.csv") %>%
  filter(!str_detect(filename, "pos5|neg5")) %>%
  select(-filename) %>%
  group_by(voltage, compound_name) %>%
  summarize(MS2 = paste(voltage, MS2, sep = "V ", collapse = ": ")) %>%
  as.data.frame() %>%
  left_join(ingalls.standards, by = "compound_name")

single.run.experimental <- read_csv("example_data/Ingalls_Lab_Standards_MSMS.csv") %>%
  filter(str_detect(filename, "pos5|neg5")) %>%
  select(-filename) %>%
  group_by(voltage, compound_name) %>%
  summarize(MS2 = paste(voltage, MS2, sep = "V ", collapse = ": ")) %>%
  as.data.frame() %>%
  left_join(ingalls.standards, by = "compound_name")

# Functions ---------------------------------------------------------------

CreatePotentialMatches_1 <- function(mz_i, rt_i, col_i, z_i, MS2str_i, ppm_error, theoretical_db) {
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
    ## TODO
    #mutate(MS2SimScore = ifelse(is.na(MS2), NA, as.numeric(lapply(MS2, CalculateMS2SimScore_1, ms2_theo = MS2str_i, flex = 0.02)))) %>%
    ##
    mutate(MS2SimScore = as.numeric(lapply(MS2, CalculateMS2SimScore_1, ms2_theo = MS2str_i, flex = 0.02))) %>%
    mutate(TotalSimScore = CalculateTotalSimScore_1(MS1SimScore, RT1SimScore, MS2SimScore)) %>%
    select(compound, ends_with("SimScore"))

  return(potential.matches)
}

MakeScantable <- function(concatenated.scan) {
  requireNamespace("dplyr", quietly = TRUE)

  scantable <- read.table(text = as.character(concatenated.scan),
                          col.names = c("mz", "intensity"), fill = TRUE) %>%
    dplyr::mutate(mz = as.numeric(mz %>% stringr::str_replace(",", "")),
                  intensity = as.numeric(intensity %>% stringr::str_replace(";", "")),
                  intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
    dplyr::filter(intensity > 0.5) %>%
    dplyr::arrange(desc(intensity))

  return(scantable)
}

CalculateMzSimScore_1 <- function(mz_exp, mz_theo, flex) {
  similarity.score = exp(-0.5 * (((mz_exp - mz_theo) / flex) ^ 2))

  return(similarity.score)
} ## TODO The flex values are hardcoded in the Calctotalsimscore1 function

CalculateMS2SimScore_1 <- function(ms2_exp, ms2_theo, flex) {
# Comparison will be between experimental and "consensus" spectra according to Horai et. al 2010
# as justification for using spectra
  scan1 <- MakeScantable(ms2_exp)
  scan2 <- MakeScantable(ms2_theo)

  weight1 <- (scan1[, "mz"] ^ 2) * (scan1[, "intensity"] ^ 0.5) # Need consensus to assign weights
  weight2 <- (scan2[, "mz"] ^ 2) * (scan2[, "intensity"] ^ 0.5)

  diff.matrix <- sapply(scan1[, "mz"], function(x) scan2[, "mz"] - x)
  same.index <- which(abs(diff.matrix) < flex, arr.ind = TRUE)
  cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
    (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))

  return(cosine.similarity)
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
# SingleOutput <- CreatePotentialMatches_1(mz_i = experimental.values[1, 2], rt_i = experimental.values[1, 3],
#                                 col_i = experimental.values[1, 4], z_i = experimental.values[1, 5],
#                                 MS2str_i = experimental.values[1, 6], ppm_error = 100000, theoretical_db = theoretical.values)

single.frame <- CreatePotentialMatches_1(mz_i = single.run.experimental$mz[1], rt_i = single.run.experimental$rt[1],
                                col_i = single.run.experimental$column[1], z_i = single.run.experimental$z[1],
                                MS2str_i = single.run.experimental$MS2[1], ppm_error = 100,
                                theoretical_db = four.runs.theoretical)


## Produces a dataframe with a column of dataframes all containing the information from SingleOutput
AllOutput <- experimental.values %>%
  slice(1:5) %>%
  drop_na() %>% ##THIS NEEDS TO GO
  rowwise() %>%
  mutate(newcol = list(CreatePotentialMatches_1(mz_i = mz, rt_i = rt, col_i = column, z_i = z,
                                          MS2str_i = MS2,
                                          ppm_error = 1000000,
                                          theoretical_db = theoretical.values)))


# Last step start here ---------------------------------------------------------
MyData <- AllOutput

# Inside one of the nested data frames
FilteredOutput <- MyData[[7]][[1]][which.max(MyData[[7]][[1]]$TotalSimScore),]

## Works but eeewwwww
UnnestedData <- MyData %>%
  unnest(newcol) %>%
  group_by(MassFeature) %>%
  top_n(1, TotalSimScore) %>%
  unique()

## Also works but also eeewwwww
for (i in 1:nrow(MyData)) {
  MyData$finalcol[i] = MyData[[7]][[i]][which.max(MyData[[7]][[i]]$TotalSimScore),]
  MyData$source[i] = "IngallsStandards"
}





  mutate(Cl1_choice = sapply(AnnotateCL1(TotalSimScoreDF)))


  # need to append single string column
# appends a "best match" column
# Wishlist: timestamp
