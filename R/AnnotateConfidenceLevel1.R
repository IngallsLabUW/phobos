#' AnnotateConfidenceLevel1
#'
#' Compare a dataframe of experimental values with a dataframe of all theoretical values from the Ingalls
#' Lab Standards sheet.
#' @importFrom magrittr "%>%"
#'
#' @param experimental.values A dataframe of unknown experimental values, created by the user.
#' It is important to note that it must be in the following format, including capitalization and class:
#' MassFeature (character), mz (numeric), rt (numeric, in seconds), column (numeric), z (numeric), MS2 (in concatenated format, character).
#' @param theoretical.values TODO A dataframe of known, theoretical values, i.e, the Ingalls Standards list.
#' Like its partner the experimental.values dataframe, the theoretical.values dataframe columns need to be in
#' the following exact format: compound (character), mz (numeric), rt (numeric, in seconds),
#' column (all caps, character), z (numeric), OPTIONAL filename (character), and MS2 (concatenated format, character).
#' @param mz.flexibility Flexibility for m/z matching between experimental and theoretical values. Usually defined as 0.02.
#' @param rt.flexibility Flexibility for retention time matching between experimental and theoretical values. Usually defined as ~ 15-30 seconds.
#'
#' @return A dataframe annotated for Confidence Level 1, including mz/rt/total similarity scores, as well as rank and rank source.
#' @export
#'
#' @examples
#' library(phobos)
#' # Isolate experimental data. Remember that phobos REQUIRES all experimental dataframes to be in the following format:
#' # - "MassFeature": Your unique mass feature, character.
#' # - "mz": The mz value, numeric.
#' # - "rt": The retention time, in seconds, numeric.
#' # - "column": Column the mass feature was run on, character.
#' #      There are two options for this variable: "HILIC" or "RP" (short for Reverse Phase).
#' #      In the Ingalls lab, cyano compounds are only ever run in positive mode, therefore the combination of an "RP" column observation
#' #      and a "1" z observation will results in the correct analysis.
#' # - "z": The polarity, numeric.
#' # - "MS2": MS2 data for those compounds that have it, character, in concatenated format.
#'
#' example_dir <- system.file("example_data", package = "phobos") ## This folder does not currently ship with the package, so these examples are broken
#' example_data <- list.files(example_dir, full.names = TRUE)
#' experimental.values <- read.csv(dir(example.dir, full.names = TRUE, pattern = "Example_Experimental"))
#' theoretical.values <- read.csv(dir(example.dir, full.names = TRUE, pattern = "Theoretical_Data"))
#' example.confidence.1 <- AnnotateConfidenceLevel1(experimental.values = experimental.values, theoretical.values = theoretical.values,
#' mz.flexibility = 0.02, rt.flexibility = 30)
#'
AnnotateConfidenceLevel1 <- function(experimental.values, theoretical.values, mz.flexibility, rt.flexibility) {

  # Expand and fix these tests
  if (all(colnames(experimental.values) == c("MassFeature", "mz", "rt", "column", "z", "MS2"))) {
    print("Columns are correctly named and ordered.")
  } else {
    stop("Please check your column names and order and try again!")
  }

  # if (class(experimental.values$z) != "numeric") {
  #   stop("Please check your column classes!")
  # }

  experimental.values <- experimental.values %>%
    dplyr::mutate(primary_key = 1:nrow(.))

  # Use a fuzzyjoin to combine experimental and theoretical values.
  Fuzzy.Join <- theoretical.values %>%
    fuzzyjoin::difference_left_join(experimental.values, by = c("mz"), max_dist = mz.flexibility) %>%
    dplyr::rename_with(., ~gsub("\\.y", "_experimental", .x)) %>%
    dplyr::rename_with(., ~gsub("\\.x", "_theoretical", .x)) %>%
    dplyr::select(MassFeature, primary_key, "compound_theoretical" = compound, mz_experimental, mz_theoretical,
                  "rt_sec_experimental" = rt_experimental, "rt_sec_theoretical" = rt_theoretical,
                  column_experimental, column_theoretical, z_experimental, z_theoretical, MS2_experimental, MS2_theoretical)  %>%
    dplyr::arrange(primary_key)

  # Confidence Level 1 ----------------------------------------
  Confidence.Level.1 <- Fuzzy.Join %>%
    dplyr::filter(z_experimental == z_theoretical,
                  column_experimental == column_theoretical) %>%
    dplyr::mutate(mz_similarity_score1 = CalculateSimilarityScore(mz_experimental, mz_theoretical, mz.flexibility),
                  rt_similarity_score1 = CalculateSimilarityScore(rt_sec_experimental, rt_sec_theoretical, rt.flexibility)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ppm_mass_error = ((abs(mz_experimental - mz_theoretical)) / mz_theoretical) * 10^6,
                  MS2_cosine_similarity1 = ifelse(is.na(MS2_experimental) | is.na(MS2_theoretical),
                                                  NA, MS2CosineSimilarity(MS2_experimental, MS2_theoretical)),
                  total_similarity_score = ifelse(is.na(MS2_cosine_similarity1),
                                                  mean(c(mz_similarity_score1, rt_similarity_score1)),
                                                  mean(c(MS2_cosine_similarity1, mz_similarity_score1, rt_similarity_score1) * 100)))

  # Sanity check -------------------------------------------------------------
  # No fuzzy match (no mz within 0.02 daltons)
  # No.Fuzzy.Match <- experimental.values %>%
  #   dplyr::select(primary_key) %>%
  #   dplyr::filter(primary_key %in% setdiff(1:nrow(.), Fuzzy.Join$primary_key)) %>%
  #   dplyr::pull()
  #
  # # Fuzzy match, but wrong z/column
  # No.CL1.Match <- setdiff(1:nrow(experimental.values),
  #                         sort(c(unique(Confidence.Level.1$primary_key),
  #                                No.Fuzzy.Match)))
  #
  # # Have any compounds been lost? Check for a TRUE output
  # all.experimentals <- sort(c(No.Fuzzy.Match, No.CL1.Match, unique(Confidence.Level.1$primary_key)))
  # length(all.experimentals) == length(experimental.values$primary_key)
  #
  # # Make "no match" dataframes for comparison
  # No.CL1.Match.df <- experimental.values %>%
  #   dplyr::filter(primary_key %in% No.CL1.Match)
  # No.Fuzzy.Match.df <- experimental.values %>%
  #   dplyr::filter(primary_key %in% No.Fuzzy.Match)
  #
  # # Let's land on MARS ------------------------------------------------------
  # Mission.Accomplished <- Confidence.Level.1 %>%
  #   dplyr::bind_rows(No.CL1.Match.df) %>%
  #   dplyr::bind_rows(No.Fuzzy.Match.df) %>%
  #   dplyr::mutate(confidence_rank = ifelse(mz_similarity_score1 > 0.9 & rt_similarity_score1 > 0.75 & ppm_mass_error < 7, 1, NA),
  #                 confidence_source = ifelse(!is.na(confidence_rank), "Ingalls_Standards", NA)) %>%
  #   dplyr::mutate(mz_experimental = ifelse(is.na(mz_experimental) & !is.na(mz), mz, mz_experimental),
  #                 rt_sec_experimental = ifelse(is.na(rt_sec_experimental) & !is.na(rt), rt, rt_sec_experimental),
  #          column_experimental = ifelse(is.na(column_experimental) & !is.na(column), column, column_experimental),
  #          z_experimental = ifelse(is.na(z_experimental) & !is.na(z), z, z_experimental),
  #          MS2_experimental = ifelse(is.na(MS2_experimental) & !is.na(MS2), MS2, MS2_experimental)) %>%
  #   dplyr::select(MassFeature, everything(), -c("mz", "rt", "column", "z", "MS2"), ) %>%
  #   dplyr::arrange(primary_key)
  #
  # return(Mission.Accomplished)
  return(Confidence.Level.1)
}
