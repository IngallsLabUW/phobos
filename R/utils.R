#' Create a filtered mini dataframe from a concatenated scanlist of MS2s.
#'
#' @param scan Single observation from a dataframe containing MS2 m/z and intensity spectra, separated by semicolons.
#'
#' @return scantable: A tiny dataframe, containing columns of mz and intensity.
#' Intensity is scaled to 100 and filtered to drop all intensity values below 0.5.
#'
#' @examples
MakeScantable <- function(scan) {
  requireNamespace("dplyr", quietly = TRUE)
  scantable <- read.table(text = as.character(scan),
                          col.names = c("mz", "intensity"), fill = TRUE) %>%
    dplyr::mutate(mz = as.numeric(mz %>% stringr::str_replace(",", "")),
                  intensity = as.numeric(intensity %>% stringr::str_replace(";", "")),
                  intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
    dplyr::filter(intensity > 0.5) %>%
    dplyr::arrange(desc(intensity))

  return(scantable)
}

MS2CosineSimilarity <- function(scan1, scan2) {
  # Finds the weighted cosine similarity between two sets of MS2 spectra.
  #
  # Args
  #   scan1 & scan2: Tiny dataframes of MS2. First column is mz, second column is intensity.
  #                  These are outputs of the MakeScantable() function.
  #
  # Returns
  #   cosine.similarity: A weighted similarity score between 0 and 1, indicating the cosine
  #                      relationship of the two vectors.
  #
  mz.tolerance <- 0.02

  weight1 <- (scan1[, 1] ^ 2) * sqrt(scan1[, 2])
  weight2 <- (scan2[, 1] ^ 2) * sqrt(scan2[, 2])

  diff.matrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x)
  same.index <- which(abs(diff.matrix) < mz.tolerance, arr.ind = TRUE)
  cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
    (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))

  return(cosine.similarity)
}
