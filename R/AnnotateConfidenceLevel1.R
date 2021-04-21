#' AnnotateConfidenceLevel1
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
#' @param mz.flexibility Flexibility for m/z matching between experimental and theoretical values. Usually defined as 0.02.
#' @param rt.flexibility Flexibility for retention time matching between experimental and theoretical values. Usually defined as ~ 15-30 seconds.
#'
#' @return A complete dataframe, annotated for Confidence Level 1. ## EXPAND THIS - "complete" is not descriptive
#' @export
#'
#' @examples
#' library(phobos)
#' # Isolate experimental data. Remember that phobos REQUIRES all experimental dataframes to be in the following format:
#' # - "MassFeature": Your unique mass feature, character.
#' # - "mz": The mz value, numeric.
#' # - "rt": The retention time, in seconds, numeric.
#' # - "column": Column the mass feature was run on, character. ##Possible options should be documented for this - HILIC vs RP? CYANO?
#' # - "z": The ion mode, numeric. ##This is not always known - we should use polarity ("pos" vs "neg") rather than absolute charge (z)
#' # - "MS2": MS2 data for those compounds that have it, character, in concatenated format.
#'
#' example.dir <- system.file("example_data", package = "phobos") ## This folder does not currently ship with the package, so these examples are broken
#' experimental.values <- read.csv(dir(example.dir, full.names = TRUE, pattern = "Example_Experimental"))
#' theoretical.values <- read.csv(dir(example.dir, full.names = TRUE, pattern = "Theoretical_Data"))
#' example.confidence.1 <- AnnotateConfidenceLevel1(experimental.values = experimental.values, theoretical.values = theoretical.values,
#' mz.flexibility = 0.02, rt.flexibility = 30)
AnnotateConfidenceLevel1 <- function(experimental.values, theoretical.values, mz.flexibility, rt.flexibility) {
  # Include tests here for correct shape and value of arguments
  if(class(experimental.values)!="data.frame"){
    stop("Looks like experimental.values is in the wrong format")
  }
  if(names(experimental.values)!=c("compound", "mz", "rt", "column", "z", "MS2")){
    stop("Looks like experimental.values is in the wrong format")
  }
  # Etc.

  # Variable names as a whole could use some innovation - what is the content of "My.Fuzzy.Join"?
  My.Fuzzy.Join <- theoretical.values %>%
    fuzzyjoin::difference_left_join(experimental.values, by = c("mz"), max_dist = mz.flexibility) %>%
    dplyr::rename_with(theoretical.values, ~gsub("\\.y", "_experimental", .x)) %>%
    dplyr::rename_with(theoretical.values, ~gsub("\\.x", "_theoretical", .x)) %>%
    # What *doesn't* this select? Would be clearer to drop columns explicitly since we're keeping so many
    dplyr::select(MassFeature, compound_experimental, compound_theoretical, mz_experimental, mz_theoretical, rt_sec_experimental, rt_sec_theoretical,
                  column_experimental, column_theoretical, z_experimental, z_theoretical, MS2_experimental, MS2_theoretical)  %>%
    dplyr::arrange(compound_experimental)

  # Confidence Level 1 ----------------------------------------
  Confidence.Level.1 <- My.Fuzzy.Join %>%
    dplyr::filter(z_experimental == z_theoretical,
                  column_experimental == column_theoretical) %>%
    dplyr::mutate(mz_similarity_score1 = exp(-0.5 * (((mz_experimental - mz_theoretical) / mz.flexibility) ^ 2)), # This should be a function (calcSimScore?)
                  rt_similarity_score1 = exp(-0.5 * (((rt_sec_experimental - rt_sec_theoretical) / rt.flexibility) ^ 2))) %>% # This is not how I expected the mz. and rt.flexibility to be used - needs more documentation
    dplyr::rowwise() %>%
    dplyr::mutate(ppm_mass_error = ((abs(mz_experimental - mz_theoretical)) / mz_theoretical) * 10^6,
                  MS2_cosine_similarity1 = ifelse(is.na(MS2_experimental) | is.na(MS2_theoretical),
                                                 NA, MS2CosineSimilarity(MakeScantable(MS2_experimental), MakeScantable(MS2_theoretical))),
                  # It might make sense to wrap MakeScanTable into MS2CosineSimilarity - the intuitive
                  # representation of the code here is MS2CosineSimilarity(MS2_expected, MS2_observed)
                  total_similarity_score = ifelse(is.na(MS2_cosine_similarity1),
                                                  mean(c(mz_similarity_score1, rt_similarity_score1)), # Use "mean" to be intuitive
                                                  mean(c(MS2_cosine_similarity1, mz_similarity_score1, rt_similarity_score1))))*100


  # Sanity check -------------------------------------------------------------
  # No fuzzy match (no mz within 0.02 daltons)
  No.Fuzzy.Match <- experimental.values %>%
    dplyr::select(compound_experimental) %>%
    # I don't think compound_experimental is always going to be a numeric? If it is, this should be either
    # created internally or specified in the docs
    dplyr::filter(compound_experimental %in% setdiff(1:nrow(.), My.Fuzzy.Join$compound_experimental)) %>%
    dplyr::pull()

  # Fuzzy match, but wrong z/column
  No.CL1.Match <- setdiff(1:nrow(experimental.values),
                          sort(c(unique(Confidence.Level.1$compound_experimental),
                                 No.Fuzzy.Match)))

  # Have any compounds been lost? Check for a TRUE output
  all.experimentals <- sort(c(No.Fuzzy.Match, No.CL1.Match, unique(Confidence.Level.1$compound_experimental)))
  if(!length(all.experimentals) == length(experimental.values$compound_experimental)){
    stop("Lengths don't match") # Add details to this stop message - what are the two lengths? Why should they match?
  }

  # Make "no match" dataframes for comparison
  # These data frames should be more intuitively created by anti_join
  No.CL1.Match.df <- experimental.values %>%
    dplyr::filter(compound_experimental %in% No.CL1.Match)
  No.Fuzzy.Match.df <- experimental.values %>%
    dplyr::filter(compound_experimental %in% No.Fuzzy.Match)

  # Let's land on MARS ------------------------------------------------------
  Mission.Accomplished <- Confidence.Level.1 %>%
    dplyr::bind_rows(No.CL1.Match.df) %>%
    dplyr::bind_rows(No.Fuzzy.Match.df) %>%
    dplyr::mutate(confidence_rank = ifelse(mz_similarity_score1 > 0.9 & rt_similarity_score1 > 0.75 & ppm_mass_error < 7, 1, NA),
                  confidence_source = ifelse(!is.na(confidence_rank), "Ingalls_Standards", NA)) %>%
    dplyr::mutate(mz_experimental = ifelse(is.na(mz_experimental) & !is.na(mz), mz, mz_experimental),
                  rt_sec_experimental = ifelse(is.na(rt_sec_experimental) & !is.na(rt), rt, rt_sec_experimental),
           column_experimental = ifelse(is.na(column_experimental) & !is.na(column), column, column_experimental),
           z_experimental = ifelse(is.na(z_experimental) & !is.na(z), z, z_experimental),
           MS2_experimental = ifelse(is.na(MS2_experimental) & !is.na(MS2), MS2, MS2_experimental)) %>%
    dplyr::select(MassFeature, everything(), -c("mz", "rt", "column", "z", "MS2"), ) %>%
    dplyr::arrange(compound_experimental)

  return(Mission.Accomplished)
}
