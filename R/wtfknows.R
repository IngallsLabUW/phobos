library(tidyverse)


experimental.values <- read.csv("example_data/Example_Experimental_Data.csv")
theoretical.values <- read.csv("example_data/Example_Theoretical_Data.csv")


CalculateTotalSimScore1 <- function(mz_i, rt_i, col_i, z_i, MS2str_i, ppm_error, theoretical_db) {
  output <- theoretical_db %>%
    filter(mz < mz_i + ((mz_i * ppm_error)/1e6) & mz > mz_i - ((mz_i * ppm_error)/1e6)) %>%
    filter(column == col_i) %>%
    filter(z == z_i) %>%
    mutate(MS21SimScore = map_chr(MS2, function(MS2exp) { # MS2 is theoretical ms2. ms2exp is epxerimental
      #print(MS2exp)
      MS21SimilarityScore(MS2exp, MS2, 0.02)
    }))}



MS21SimilarityScore <- function(ms2_exp, ms2_theo, flex) {
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


## Produces dataframe of potential matches and all sim scores.
OutputDF <- CalculateTotalSimScore1(mz_i = experimental.values[1, 2], rt_i = experimental.values[1, 3],
                                    col_i = experimental.values[1, 4], z_i = experimental.values[1, 5],
                                    MS2str_i = experimental.values[1, 6], ppm_error = 100000, theoretical_db = theoretical.values)

