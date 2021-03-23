#' Compare experimental mass features to scraped MoNA data, based on MS1 and MS2 information.
#'
#' @param MoNA.Mass Single MoNA mass, isolated from scraped MoNA df.
#' @param experimental.df experimental dataframe, islated to contain only MS1 and MS2 data.
#'
#' @return final.candidates: dataframe of experimental data, matched with information from MoNA.
#' @export
#'
#' @examples
IsolateMoNACandidates <- function(MoNA.Mass, experimental.df, potential.candidates) {
  potential.candidates <- potential.candidates %>%
    dplyr::filter(MH_mass > MoNA.Mass - 0.020,
           MH_mass < MoNA.Mass + 0.020) %>%
    fuzzyjoin::difference_inner_join(experimental.df, by = "MH_mass", max_dist = 0.02) %>%
    dplyr::rename(scan1 = spectrum_KRHform_filtered, # scan1 is MS2 from MoNA
                  scan2 = MS2_experimental,          # scan2 is MS2 from the experimental data
                  mass1 = MH_mass.x,                 # mass1 is the mass from MoNA
                  mass2 = MH_mass.y)                 # mass2 is the mass from experimental data

  if (length(potential.candidates$ID) == 0) {
    print("There are no potential candidates.")
    No.Match.Return <- Mass.Feature %>%
      dplyr::mutate(massbank_match = NA,
                    massbank_ppm = NA,
                    massbank_cosine_similarity = NA)

    return(No.Match.Return)
  }

  # Add cosine similarity scores
  print("Making potential candidates")

  potential.candidates$massbank_cosine_similarity <- apply(potential.candidates, 1, FUN = function(x) MakeMS2CosineDataframe(x))

  final.candidates <- potential.candidates %>%
    dplyr::mutate(massbank_match = paste(Names, ID, sep = " ID:"),
                  massbank_ppm = abs(mass2 - mass1) / mass1 * 10^6) %>%
    dplyr::rename(MS2_massbank = scan1,
                  mz_massbank = mass1, # should maybe name it something else because of mh?...
                  MS2_experimental = scan2,
                  mz_experimental = mass2) %>%
    unique() %>%
    dplyr::filter(massbank_ppm < 5,
                  massbank_cosine_similarity > 0.5) %>%
    dplyr::arrange(desc(massbank_cosine_similarity))

  return(final.candidates)
}

#' Create a small dataframe of
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
MakeMS2CosineDataframe <- function(df) {
  scan1 <- MakeScantable(df["scan1"])
  scan2 <- MakeScantable(df["scan2"])
  # mass1 <- df["mass1"]
  # mass2 <- df["mass2"]

  mz.tolerance <- 0.02

  similarity <- MS2CosineSimilarity(scan1 = scan1, scan2 = scan2)
  return(similarity)

  # weight1 <- (scan1[, 1] ^ 2) * sqrt(scan1[, 2])
  # weight2 <- (scan2[, 1] ^ 2) * sqrt(scan2[, 2])
  #
  # difference.matrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x)
  # same.index <- which(abs(difference.matrix) < mz.tolerance, arr.ind = TRUE)
  #
  # similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
  #   (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))
  #
  # return(similarity)
}

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

#' Finds the weighted cosine similarity between two sets of MS2 spectra.
#'
#' @param scan1 Tiny dataframe of MS2, from first set of values to be compared. Column 1 is mz, column 2 is intensity. Dataframe is the output of the MakeScantable() function.
#' @param scan2 Tiny dataframe of MS2, from second set of values to be compared. Column 1 is mz, column 2 is intensity. Dataframe is the output of the MakeScantable() function.
#'
#' @return cosine.similarity: A weighted similarity score between 0 and 1, indicating the cosinerelationship of the two vectors.
#'
#' @examples
MS2CosineSimilarity <- function(scan1, scan2) {
  mz.tolerance <- 0.02

  weight1 <- (scan1[, 1] ^ 2) * sqrt(scan1[, 2])
  weight2 <- (scan2[, 1] ^ 2) * sqrt(scan2[, 2])

  diff.matrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x)
  same.index <- which(abs(diff.matrix) < mz.tolerance, arr.ind = TRUE)
  cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
    (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))

  return(cosine.similarity)
}

#' Remove character values of "NA;" found in MARS outputs.
#'
#' @param column Character column that contains one or more "NA;" values.
#'
#' @return
#' @export
#'
#' @examples
ReplaceNA <- function(column) {
  gsub("NA; ", "", column)
}
