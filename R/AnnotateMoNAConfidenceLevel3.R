#' Annotate Confidence Level 3: MoNA
#'
#' @param Confidence.Level.2 A dataframe of experimental values, already annotated for Confidence Levels 1 & 2.
#' Please refer to the README or documentation for the AnnotateConfidenceLevel1 and AnnotateMoNAConfidenceLevel2 functions for details.
#' @param MassBank.Neg Spectra in negative ion mode from MassBank of North America, scraped and downloaded. A modified example is included in this package,
#' but the complete csv is available on the Ingalls Shared Drive in the MARS_Project folder, titled NEG_Spectra.csv.
#' @param MassBank.Pos Spectra in positive ion mode from MassBank of North America, scraped and downloaded. A modified example is included in this package,
#' but the complete csv is available on the Ingalls Shared Drive in the MARS_Project folder, titled POS_Spectra.csv.
#' @param mz.flexibility Flexibility for m/z matching between experimental and theoretical values. Usually defined as 0.02.
#'
#' @return A dataframe annotated for Confidence Level 3 by referencing MassBank of North America.
#' @export
#'
#' @examples
#' library(phobos)
#' # Load experimental data that has already been annotated for Confidence Levels 1 & 2
#' # using the AnnotateConfidenceLevel1 and AnnotateMoNAConfidenceLevel2 function2.
#' # Remember that dataframe format is important: please refer to README or the documentation for AnnotateConfidenceLevel1() for details.
#'
#' Confidence.Level.2 <- read.csv("example_data/Example_ConfidenceLevel2.csv")
#' MassBank.Neg <- read.csv("example_data/NEG_Spectra.csv")
#' MassBank.Pos <- read.csv("example_data/POS_Spectra.csv")
#'
#' Example_ConfidenceLevel3 <- AnnotateMoNaConfidenceLevel3(Confidence.Level.2 = Confidence.Level.2,
#' MassBank.Neg = MassBank.Neg, MassBank.Pos = MassBank.Pos, mz.flexibility = 0.02)
AnnotateMoNAConfidenceLevel3 <- function(Confidence.Level.2, MassBank.Neg, MassBank.Pos, mz.flexibility) {

  Experimental.Values <- Confidence.Level.2 %>%
    dplyr::select(MassFeature:compound_theoretical, mz_experimental, z_experimental, confidence_rank, confidence_source) %>%
    unique() %>%
    dplyr::group_by(primary_key) %>%
    dplyr::add_tally() %>%
    dplyr::mutate(temp = ifelse(n == 2 & is.na(confidence_rank), TRUE, FALSE)) %>%
    dplyr::filter(temp != TRUE) %>%
    dplyr::rename(MH_mass = mz_experimental) %>%
    dplyr::select(-n, -temp)

  MoNA.Spectra.NEG <- MassBank.Neg %>%
    dplyr::mutate(z_massbank3 = -1)
  MoNA.Spectra.Pos <- MassBank.Pos %>%
    dplyr::mutate(z_massbank3 = 1)

  MoNA.Spectra <- MoNA.Spectra.NEG %>%
    rbind(MoNA.Spectra.Pos) %>%
    dplyr::mutate(ID = paste("ID:", ID)) %>%
    tidyr::unite(massbank_match3, c(ID, Names), sep = "; ") %>%
    dplyr::mutate(MH_mass = M_mass - 1.0072766) %>%
    dplyr::select(massbank_match3, MH_mass, z_massbank3) %>%
    dplyr::mutate_all(., list(~dplyr::na_if(., ""))) %>%
    tidyr::drop_na()

  Fuzzy.Join <- MoNA.Spectra %>%
    fuzzyjoin::difference_left_join(Experimental.Values, by = c("MH_mass"), max_dist = 0.02) %>%
    dplyr::filter(z_massbank3 == z_experimental) %>%
    dplyr::rename_with(., ~gsub("\\.x", "_MoNA", .x)) %>%
    dplyr::rename_with(., ~gsub("\\.y", "_experimental", .x)) %>%
    dplyr::mutate(mz_similarity_score3 = CalculateSimilarityScore(MH_mass_experimental, MH_mass_MoNA, mz.flexibility)) %>%
    dplyr::select(primary_key, massbank_match3, MH_mass_experimental, MH_mass_MoNA,
                  z_experimental, z_massbank3, mz_similarity_score3, confidence_rank, confidence_source) %>%
    dplyr::arrange(primary_key)

  ## Sanity check
  No.Fuzzy.Match <- setdiff(1:nrow(Experimental.Values),
                            sort(unique(Fuzzy.Join$primary_key)))
  Fuzzy.Match <- sort(unique(Fuzzy.Join$primary_key))

  # Have any compounds been lost? Check for a TRUE output
  all.experimentals <- sort(c(No.Fuzzy.Match, Fuzzy.Match))
  length(all.experimentals) == length(Experimental.Values$primary_key)

  Full.Join <- Fuzzy.Join %>%
    dplyr::full_join(Confidence.Level.2,
                     by = c("primary_key", "z_experimental", "confidence_rank", "confidence_source")) %>%
    dplyr::arrange(primary_key) %>%
    unique() %>%
    dplyr::mutate(ppm_mass_error3 = ((abs(MH_mass_experimental - MH_mass_MoNA)) / mz_theoretical) * 10^6) %>%
    dplyr::mutate(confidence_rank3 = ifelse(mz_similarity_score3 > 0.9 & ppm_mass_error3 < 7, 3, NA)) %>%
    dplyr::mutate(confidence_rank = ifelse(!is.na(confidence_rank) & !is.na(confidence_rank3), paste(confidence_rank, "3", sep = "; "),
                                           ifelse(!is.na(confidence_rank3), confidence_rank3, confidence_rank)))  %>%
    dplyr::mutate(confidence_source = ifelse(stringr::str_detect(confidence_rank, "3"),
                                             apply(cbind(confidence_source, "MassBank"), 1, function(x) paste(x[!is.na(x)], collapse = "; ")),
                                             confidence_source))

  Confidence.Level.3 <- Full.Join %>%
    unique() %>%
    dplyr::select(primary_key, MassFeature, compound_theoretical, massbank_match2, massbank_match3,
                  mz_experimental:mz_massbank2, MH_mass_experimental, MH_mass_MoNA,
                  z_experimental, z_theoretical, z_massbank2, z_massbank3, rt_sec_experimental:column_theoretical,
                  MS2_experimental:MS2_massbank, ppm_mass_error1, ppm_mass_error3, massbank_ppm, mz_similarity_score1, mz_similarity_score2, mz_similarity_score3,
                  rt_similarity_score1:total_similarity_score2, everything(), -confidence_rank3) %>%
    dplyr::arrange(primary_key)

  return(Confidence.Level.3)
}
