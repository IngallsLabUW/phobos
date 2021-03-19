#' MatchConfidenceLevel1
#'
#' Compare a dataframe of experimental values with a dataframe of all theoretical values from the Ingalls
#' Lab Standards sheet.
#' @importFrom magrittr "%>%"
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
    fuzzyjoin::difference_left_join(experimental.values, by = c("mz"), max_dist = mz.flexibility) %>%
    dplyr::rename(compound_theoretical = compound,
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
    dplyr::select(compound_experimental, KRH_identification, compound_theoretical, mz_experimental, mz_theoretical, rt_sec_experimental, rt_sec_theoretical,
                  column_experimental, column_theoretical, z_experimental, z_theoretical, MS2_experimental, MS2_theoretical)  %>%
    dplyr::arrange(compound_experimental)

  # Confidence Level 1 ----------------------------------------
  Confidence.Level.1 <- My.Fuzzy.Join %>%
    dplyr::filter(z_experimental == z_theoretical,
                  column_experimental == column_theoretical) %>%
    dplyr::mutate(mz_similarity_score = exp(-0.5 * (((mz_experimental - mz_theoretical) / mz.flexibility) ^ 2)),
                  rt_similarity_score = exp(-0.5 * (((rt_sec_experimental - rt_sec_theoretical) / rt.flexibility) ^ 2))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ppm_mass_error = ((abs(mz_experimental - mz_theoretical)) / mz_theoretical) * 10^6,
                  MS2_cosine_similarity = ifelse(is.na(MS2_experimental) | is.na(MS2_theoretical),
                                                 NA, MS2CosineSimilarity(MakeScantable(MS2_experimental), MakeScantable(MS2_theoretical))),
                  total_similarity_score = ifelse(is.na(MS2_cosine_similarity),
                                                  ((mz_similarity_score + rt_similarity_score) / 2) * 100,
                                                  ((MS2_cosine_similarity + mz_similarity_score + rt_similarity_score) / 3) * 100))


  # Sanity check -------------------------------------------------------------
  # No fuzzy match (no mz within 0.02 daltons)
  No.Fuzzy.Match <- experimental.values %>%
    dplyr::select(compound_experimental) %>%
    dplyr::filter(compound_experimental %in% setdiff(1:nrow(.), My.Fuzzy.Join$compound_experimental)) %>%
    dplyr::pull()

  # Fuzzy match, but wrong z/column
  No.CL1.Match <- setdiff(1:nrow(experimental.values),
                          sort(c(unique(Confidence.Level.1$compound_experimental),
                                 No.Fuzzy.Match)))

  # Have any compounds been lost? Check for a TRUE output
  all.experimentals <- sort(c(No.Fuzzy.Match, No.CL1.Match, unique(Confidence.Level.1$compound_experimental)))
  length(all.experimentals) == length(experimental.values$compound_experimental)

  # Make "no match" dataframes for comparison
  No.CL1.Match.df <- experimental.values %>%
    dplyr::filter(compound_experimental %in% No.CL1.Match)
  No.Fuzzy.Match.df <- experimental.values %>%
    dplyr::filter(compound_experimental %in% No.Fuzzy.Match)

  # Let's land on MARS ------------------------------------------------------
  Mission.Accomplished <- Confidence.Level.1 %>%
    dplyr::bind_rows(No.CL1.Match.df) %>%
    dplyr::bind_rows(No.Fuzzy.Match.df) %>%
    dplyr::mutate(confidence_rank = ifelse(mz_similarity_score == 1 & rt_similarity_score == 1 & ppm_mass_error < 7, 1, NA),
                  confidence_source = ifelse(!is.na(confidence_rank), "Ingalls_Standards", NA)) %>%
    dplyr::mutate(mz_experimental = ifelse(is.na(mz_experimental) & !is.na(mz), mz, mz_experimental),
                  rt_sec_experimental = ifelse(is.na(rt_sec_experimental) & !is.na(rt), rt, rt_sec_experimental),
           column_experimental = ifelse(is.na(column_experimental) & !is.na(column), column, column_experimental),
           z_experimental = ifelse(is.na(z_experimental) & !is.na(z), z, z_experimental),
           MS2_experimental = ifelse(is.na(MS2_experimental) & !is.na(MS2), MS2, MS2_experimental)) %>%
    dplyr::select(-c("mz", "rt", "column", "z", "MS2")) %>%
    dplyr::arrange(compound_experimental)

  return(Mission.Accomplished)
}
