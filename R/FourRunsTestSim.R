library(tidyverse)

## Theoretical data created from four standard runs combined.
## Experimental data created from the fifth standard run.
FourRunsTheoretical <- read_csv("example_data/Ingalls_Lab_Standards_MSMS.csv") %>%
  filter(!str_detect(filename, "pos5|neg5"))
SingleRunExperimental <- read_csv("example_data/Ingalls_Lab_Standards_MSMS.csv") %>%
  filter(str_detect(filename, "pos5|neg5"))

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


testimscore <- MS21SimilarityScore(SingleRunExperimental[1,2], FourRunsTheoretical[1,2], 0.02)
