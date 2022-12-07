library(tidyverse)


experimental.values <- read.csv("example_data/Example_Experimental_Data.csv")
theoretical.values <- read.csv("example_data/Example_Theoretical_Data.csv") %>%
  drop_na()

# All functions -----------------------------------------------------------
CalculateTotalSimScore1 <- function(mz_i, rt_i, col_i, z_i, MS2_i, ppm_error, theoretical_db) {
  output <- theoretical_db %>%
    filter(mz < mz_i + ((mz_i * ppm_error)/1e6) & mz > mz_i - ((mz_i * ppm_error)/1e6)) %>%
    filter(column == col_i) %>%
    filter(z == z_i) %>%
    mutate(MS1SimScore = MS1SimilarityScore(mz_exp = mz_i, mz_theo = mz, flex = 5)) %>%
    mutate(MS2_SimScore = lapply(MS2, MS2SimilarityScore, ms2_theo = MS2_i, flex = 0.02))
  }

MakeScantable <- function(concatenated.scan) {
  requireNamespace("dplyr", quietly = TRUE)

  print(concatenated.scan)
  scantable <- read.table(text = as.character(concatenated.scan),
                          col.names = c("mz", "intensity"), fill = TRUE) %>%
    dplyr::mutate(mz = as.numeric(mz %>% stringr::str_replace(",", "")),
                  intensity = as.numeric(intensity %>% stringr::str_replace(";", "")),
                  intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
    dplyr::filter(intensity > 0.5) %>%
    dplyr::arrange(desc(intensity))

  return(scantable)
}

MS1SimilarityScore <- function(mz_exp, mz_theo, flex) {
  similarity.score = exp(-0.5 * (((mz_exp - mz_theo) / flex) ^ 2))

  return(similarity.score)
}

MS2SimilarityScore <- function(ms2_exp, ms2_theo, flex) {

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

test <- experimental.values %>%
  mutate(newcol = lapply(MS2, MakeScantable))


# Create output dataframe -------------------------------------------------
# Produces dataframe of potential matches and all sim scores.
OutputDF <- CalculateTotalSimScore1(mz_i = experimental.values[1, 2], rt_i = experimental.values[1, 3],
                                    col_i = experimental.values[1, 4], z_i = experimental.values[1, 5],
                                    MS2_i = experimental.values[1, 6], ppm_error = 100000,
                                    theoretical_db = theoretical.values)

t <- MS2SimilarityScore(ms2_exp = experimental.values[1, 6], ms2_theo = theoretical.values[586, 7], flex = 0.02)
