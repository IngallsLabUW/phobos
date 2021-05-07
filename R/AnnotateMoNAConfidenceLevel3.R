#' Annotate Confidence Level 3: MoNA
#'
#' @param Confidence.Level.2 A dataframe of experimental values, already annotated for Confidence Levels 1 & 2.
#' Please refer to the README or documentation for the AnnotateConfidenceLevel1 and AnnotateMoNAConfidenceLevel2 functions for details.
#' @param MassBank.Neg Spectra in negative ion mode from MassBank of North America, scraped and downloaded. This csv is available on the Ingalls Shared Drive in the MARS_Project folder, titled NEG_Spectra.csv
#' @param MassBank.Pos Spectra in positive ion mode from MassBank of North America, scraped and downloaded. This csv is available on the Ingalls Shared Drive in the MARS_Project folder, titled POS_Spectra.csv
#' @param mz.flexibility Flexibility for m/z matching between experimental and theoretical values. Usually defined as 0.02.
#'
#' @return A complete dataframe, annotated for Confidence Level 3 using MoNA.
#' @export
#'
#' @examples
#' library(phobos)
#' # Load experimental data that has already been annotated for Confidence Levels 1 & 2
#' # using the AnnotateConfidenceLevel1 and AnnotateMoNAConfidenceLevel2 function2.
#' # Remember that dataframe format is important: please refer to README or the documentation for AnnotateConfidenceLevel1() for details.
#'
#' Confidence.Level.2 <- read.csv("example_data/Example_ConfidenceLevel2.csv")
#' MassBank.Neg <- read.csv("your/file/path/here/NEG_Spectra.csv")
#' MassBank.Pos <- read.csv("your/file/path/here/POS_Spectra.csv")
#'
#' Example_ConfidenceLevel3 <- AnnotateMoNaConfidenceLevel3(Confidence.Level.2 = Confidence.Level.2,
#' MassBank.Neg = MassBank.Neg, MassBank.Pos = MassBank.Pos, mz.flexibility = 0.02)
AnnotateMoNAConfidenceLevel3 <- function(Confidence.Level.2, MassBank.Neg, MassBank.Pos, mz.flexibility) {

  Confidence.Level.2 <- Confidence.Level.2

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
    dplyr::mutate(z_MoNA = -1)
  MoNA.Spectra.Pos <- MassBank.Pos %>%
    dplyr::mutate(z_MoNA = 1)

  MoNA.Spectra <- MoNA.Spectra.NEG %>%
    rbind(MoNA.Spectra.Pos) %>%
    dplyr::mutate(MH_mass = M_mass - 1.0072766) %>%
    dplyr::select(ID, Names, MH_mass, z_MoNA) %>%
    dplyr::mutate_all(., list(~dplyr::na_if(.,""))) %>%
    tidyr::drop_na()

  My.Fuzzy.Join <- MoNA.Spectra %>%
    fuzzyjoin::difference_left_join(Experimental.Values, by = c("MH_mass"), max_dist = 0.02) %>%
    dplyr::filter(#is.na(confidence_source),
      z_MoNA == z_experimental) %>%
    dplyr::rename(MH_mass_MoNA = MH_mass.x,
                  MH_mass_experimental = MH_mass.y) %>%
    dplyr::mutate(mz_similarity_score3 = exp(-0.5 * (((MH_mass_experimental - MH_mass_MoNA) / mz.flexibility) ^ 2))) %>%
    dplyr::select(primary_key, ID, Names, MH_mass_experimental, MH_mass_MoNA,
                  z_experimental, z_MoNA, mz_similarity_score3, confidence_rank, confidence_source) %>%
    dplyr::arrange(primary_key)

  No.Fuzzy.Match <- setdiff(1:nrow(Experimental.Values),
                            sort(unique(My.Fuzzy.Join$primary_key)))
  CL3.Match <- as.integer(unique(My.Fuzzy.Join$primary_key))

  # Have any compounds been lost? Check for a TRUE output
  all.experimentals <- sort(c(No.Fuzzy.Match, CL3.Match))
  length(all.experimentals) == length(Experimental.Values$primary_key)

  # Make "no match" dataframes for comparison
  # No.Fuzzy.Match.df <- Experimental.Values %>%
  #   filter(primary_key %in% No.Fuzzy.Match) %>%
  #   dplyr::rename(MH_mass_experimental = MH_mass)

  No.CL3.Match.df <- experimental.values %>%
    dplyr::anti_join(Confidence.Level.1) %>%
    dplyr::rename("rt_sec_experimental" = rt)


  final <- My.Fuzzy.Join %>%
    bind_rows(No.Fuzzy.Match.df) %>%
    dplyr::mutate(confidence_rank3 = ifelse(mz_similarity_score3 > 0.9, 3, NA)) %>%
    dplyr::mutate(confidence_rank = ifelse(!is.na(confidence_rank) & !is.na(confidence_rank3), paste(confidence_rank, "3", sep = "; "),
                                           ifelse(!is.na(confidence_rank3), confidence_rank3, confidence_rank)))  %>%
    dplyr::mutate(confidence_source = ifelse(str_detect(confidence_rank, "3"),
                                             apply(cbind(confidence_source, "MoNA"), 1, function(x) paste(x[!is.na(x)], collapse = "; ")),
                                             confidence_source)) %>%
    dplyr::select(primary_key, "MoNA_ID" = ID, "MoNA_Names" = Names, MH_mass_experimental, MH_mass_MoNA,
                  z_experimental, z_MoNA, mz_similarity_score3, confidence_rank, confidence_source)

  Confidence.Level.3 <- final %>%
    dplyr::left_join(Confidence.Level.2) %>%
    unique() %>%
    dplyr::select(primary_key, compound_theoretical, MoNA_ID, MoNA_Names, massbank_match, everything()) %>%
    dplyr::arrange(primary_key)

  return(Confidence.Level.3)
}
