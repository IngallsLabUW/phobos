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
