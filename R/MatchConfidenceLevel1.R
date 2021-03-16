#' MatchConfidenceLevel1
#'
#' Compare a dataframe of experimental values with a dataframe of all theoretical values from the Ingalls
#' Lab Standards sheet.
#'
#' @param experimental.values A dataframe of unknown experimental values, created by the user.
#' It is important to note that it must be in the following format, including capitalization and class:
#' mz (numeric), rt (numeric, in seconds), column (numeric), z (numeric), MS2 (in concatenated format, character).
#' @param theoretical.values TODO A dataframe of known, theoretical values, i.e, the Ingalls Standards list.
#' Like its partner the experimental.values dataframe, the theoretical.values dataframe columns need to be in
#' the following exact format: compound (character), mz (numeric), rt (numeric, in seconds),
#' column (all caps, character), z (numeric), OPTIONAL filename (character), and MS2 (concatenated format, character).
#' @param mz.flexibility TODO Usually defined as 0.02.
#' @param rt.flexibility Usually defined as ~ 15-30 seconds.
#'
#' @return A complete dataframe, annotated for Confidence Level 1.
#' @export
#'
#' @examples
MatchConfidenceLevel1 <- function(experimental.values, theoretical.values, mz.flexibility, rt.flexibility) {
  My.Fuzzy.Join <- theoretical.values %>%
    difference_left_join(experimental.values, by = c("mz"), max_dist = mz.flexibility) %>%
    rename(compound_theoretical = compound,
           mz_theoretical = mz.x,
           rt_sec_theoretical = rt.x,
           column_theoretical = column.x,
           z_theoretical = z.x,
           MS2_theoretical = MS2.x,
           mz_experimental = mz.y,
           rt_sec_experimental = rt.y,
           column_experimental = column.y,
           z_experimental = z.y,
           MS2_experimental = MS2.y) %>%
    select(compound_experimental, KRH_identification, compound_theoretical, mz_experimental, mz_theoretical, rt_sec_experimental, rt_sec_theoretical,
           column_experimental, column_theoretical, z_experimental, z_theoretical, MS2_experimental, MS2_theoretical)  %>%
    arrange(compound_experimental)

  # Confidence Level 1 ----------------------------------------
  Confidence.Level.1 <- My.Fuzzy.Join %>%
    filter(z_experimental == z_theoretical,
           column_experimental == column_theoretical) %>%
    mutate(mz_similarity_score = exp(-0.5 * (((mz_experimental - mz_theoretical) / mz.flexibility) ^ 2)),
           rt_similarity_score = exp(-0.5 * (((rt_sec_experimental - rt_sec_theoretical) / rt.flexibility) ^ 2))) %>%
    rowwise() %>%
    mutate(ppm_mass_error = ((abs(mz_experimental - mz_theoretical)) / mz_theoretical) * 10^6) %>%
    mutate(MS2_cosine_similarity = ifelse(is.na(MS2_experimental) | is.na(MS2_theoretical),
                                          NA, MS2CosineSimilarity(MakeScantable(MS2_experimental), MakeScantable(MS2_theoretical)))) %>%
    mutate(total_similarity_score = ifelse(is.na(MS2_cosine_similarity),
                                           ((mz_similarity_score + rt_similarity_score) / 2) * 100,
                                           ((MS2_cosine_similarity + mz_similarity_score + rt_similarity_score) / 3) * 100))


  # Sanity check -------------------------------------------------------------
  # No fuzzy match (no mz within 0.02 daltons)
  No.Fuzzy.Match <- experimental.values %>%
    select(compound_experimental) %>%
    filter(compound_experimental %in% setdiff(1:nrow(.), My.Fuzzy.Join$compound_experimental)) %>%
    pull()

  # Fuzzy match, but wrong z/column
  No.CL1.Match <- setdiff(1:nrow(experimental.values),
                          sort(c(unique(Confidence.Level.1$compound_experimental),
                                 No.Fuzzy.Match)))

  # Have any compounds been lost? Check for a TRUE output
  all.experimentals <- sort(c(No.Fuzzy.Match, No.CL1.Match, unique(Confidence.Level.1$compound_experimental)))
  length(all.experimentals) == length(experimental.values$compound_experimental)

  # Make "no match" dataframes for comparison
  No.CL1.Match.df <- experimental.values %>%
    filter(compound_experimental %in% No.CL1.Match)
  No.Fuzzy.Match.df <- experimental.values %>%
    filter(compound_experimental %in% No.Fuzzy.Match)

  # Let's land on MARS ------------------------------------------------------
  Mission.Accomplished <- Confidence.Level.1 %>%
    bind_rows(No.CL1.Match.df) %>%
    bind_rows(No.Fuzzy.Match.df) %>%
    mutate(confidence_rank = ifelse(mz_similarity_score == 1 & rt_similarity_score == 1 & ppm_mass_error < 7, 1, NA),
           confidence_source = ifelse(!is.na(confidence_rank), "Ingalls_Standards", NA)) %>%
    mutate(mz_experimental = ifelse(is.na(mz_experimental) & !is.na(mz), mz, mz_experimental),
           rt_sec_experimental = ifelse(is.na(rt_sec_experimental) & !is.na(rt), rt, rt_sec_experimental),
           column_experimental = ifelse(is.na(column_experimental) & !is.na(column), column, column_experimental),
           z_experimental = ifelse(is.na(z_experimental) & !is.na(z), z, z_experimental),
           MS2_experimental = ifelse(is.na(MS2_experimental) & !is.na(MS2), MS2, MS2_experimental)) %>%
    select(-c("mz", "rt", "column", "z", "MS2")) %>%
    arrange(compound_experimental)

  return(Mission.Accomplished)
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

#' Create a filtered mini dataframe from a concatenated scanlist of MS2s.
#'
#' @param scan Single observation from a dataframe containing MS2 m/z and intensity spectra, separated by semicolons.
#'
#' @return scantable: A tiny dataframe, containing columns of mz and intensity.
#' Intensity is scaled to 100 and filtered to drop all intensity values below 0.5.
#' @export
#'
#' @examples
MakeScantable <- function(scan) {
  requireNamespace("dplyr", quietly = TRUE)
  scantable <- read.table(text = as.character(scan),
                          col.names = c("mz", "intensity"), fill = TRUE) %>%
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", "")),
           intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
    filter(intensity > 0.5) %>%
    arrange(desc(intensity))

  return(scantable)
}

