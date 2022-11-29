library(tidyverse)

### WKumler + RML Startup script

# Begin with experimental dataframe

experimental.values <- read.csv("example_data/Example_Experimental_Data.csv")
theoretical.values <- read.csv("example_data/Example_Theoretical_Data.csv")

# pass experimental df to CalcTotSim1, which takes each observation in each column as an argument, rt, mz, ms2str, col, z)
# Takes the theoretical value and ensures we are matching within mz window, z, and column.
# Calculates similarity scores for rt, mz, ms2.

TotalSimilarityScore <- function(ms1_sim, rt_sim, ms2_sim) {
  total.similarity.score <- ((ms1_sim + rt_sim + ms2_sim) / 3) * 100

  return(total.similarity.score)
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

MS1SimilarityScore <- function(mz_exp, mz_theo, flex = 5) {
  similarity.score = exp(-0.5 * (((mz_exp - mz_theo) / flex) ^ 2))

  return(similarity.score)
}

RTSimilarityScore <- function(rt_exp, rt_theo, flex = 30) {
  similarity.score = exp(-0.5 * (((rt_exp - rt_theo) / flex) ^ 2))

  return(similarity.score)
}

MS21SimilarityScore <- function(ms2_exp, ms2_theo, flex) {

  scan1 <- MakeScantable(ms2_exp)
  scan2 <- MakeScantable(ms2_theo)

  weight1 <- (scan1[, "mz"] ^ 2) * sqrt(scan1[, "intensity"])
  weight2 <- (scan2[, "mz"] ^ 2) * sqrt(scan2[, "intensity"])

  diff.matrix <- sapply(scan1[, "mz"], function(x) scan2[, "mz"] - x)
  same.index <- which(abs(diff.matrix) < flex, arr.ind = TRUE)
  cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
    (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))

  return(cosine.similarity)
}

CalculateTotalSimScore1 <- function(mz_i, rt_i, col_i, z_i, MS2str_i, ppm_error, theoretical_db) {
  output <- theoretical_db %>%
    filter(mz < mz_i + ((mz_i * ppm_error)/1e6) & mz > mz_i - ((mz_i * ppm_error)/1e6)) %>% ## 50 is ppm error
    filter(column == col_i) %>%
    filter(z == z_i) %>%
    mutate(MS1SimScore = MS1SimilarityScore(mz_exp = mz_i, mz_theo = mz, flex = 5)) %>%
    mutate(RT1SimScore = RTSimilarityScore(rt_exp = rt_i, rt_theo = rt, flex = 30)) %>%
    mutate(MS21SimScore = MS21SimilarityScore(ms2_exp = MS2str_i, ms2_theo = MS2[1], flex = 0.02)) %>% # This flex is totally made up, and is still the old function
    mutate(TotalSimScore = TotalSimilarityScore(MS1SimScore, RT1SimScore, MS21SimScore)) %>%
    select(compound, ends_with("SimScore")) # originally we had "name_db" here, did we want that to just be the compound name/whatever "id" is used for a feature?

  return(output)
}

OutputDF <- CalculateTotalSimScore1(mz_i = experimental.values[1, 2], rt_i = experimental.values[1, 3],
                                col_i = experimental.values[1, 4], z_i = experimental.values[1, 5],
                                MS2str_i = experimental.values[1, 6], ppm_error = 100000, theoretical_db = theoretical.values)

# Outputs dataframe of original unknown feature, with accompanying columns of all sim scores, and all matches.

experimental.values %>%
  mutate(TotalSimScoreDF = CalculateTotalSimScore1(mz, rt, col, z, MS2str, theoretical.values)) %>%
  mutate(Cl1_choice = sapply(AnnotateCL1(TotalSimScoreDF)))


## Annotate CL1 function needs to be written, actually makes the confidence rank decision and source (aka theoretical_db). Wishlist: timestamp
