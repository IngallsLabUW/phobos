#' Annotate Confidence Level 3: MoNA
#'
#' @param Confidence.Level.2 A dataframe of experimental values, already annotated for Confidence Levels 1 & 2.
#' Please refer to the README or documentation for the AnnotateConfidenceLevel1 and AnnotateMoNAConfidenceLevel2 functions for details.
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
#' example_dir <- system.file("example_data", package = "phobos")
#' example_data <- list.files(example_dir, full.names = TRUE)
#' Confidence.Level.2 <- read.csv(grep("Example_ConfidenceLevel2", example_data, value = TRUE))
#'
#' example_confidenceL3 <- AnnotateMoNaConfidenceLevel3(Confidence.Level.2 = Confidence.Level.2, mz.flexibility = 0.02)
AnnotateMoNAConfidenceLevel3 <- function(Confidence.Level.2, mz.flexibility) {

  Confidence.Level.2 <- Confidence.Level.2

  Experimental.Values <- Confidence.Level.2 %>%
    select(compound_experimental:KRH_identification, mz_experimental, z_experimental, confidence_rank, confidence_source) %>%
    unique() %>%
    group_by(compound_experimental) %>%
    add_tally() %>%
    mutate(temp = ifelse(n == 2 & is.na(confidence_rank), TRUE, FALSE)) %>%
    filter(temp != TRUE) %>%
    rename(MH_mass = mz_experimental) %>%
    select(-n, -temp)



  MoNA.Spectra.NEG <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Spectra.csv") %>%
    mutate(z_MoNA = -1)
  MoNA.Spectra.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Spectra.csv") %>%
    mutate(z_MoNA = 1)

  MoNA.Spectra <- MoNA.Spectra.NEG %>%
    rbind(MoNA.Spectra.Pos) %>%
    mutate(MH_mass = M_mass - 1.0072766) %>%
    select(ID, Names, MH_mass, z_MoNA) %>%
    mutate_all(., list(~na_if(.,""))) %>%
    drop_na()


  My.Fuzzy.Join <- MoNA.Spectra %>%
    difference_left_join(Experimental.Values, by = c("MH_mass"), max_dist = 0.02) %>%
    filter(is.na(confidence_source),
           z_MoNA == z_experimental) %>%
    rename(MH_mass_MoNA = MH_mass.x,
           MH_mass_experimental = MH_mass.y) %>%
    mutate(mz_similarity_score = exp(-0.5 * (((MH_mass_experimental - MH_mass_MoNA) / mz.flexibility) ^ 2))) %>%
    select(compound_experimental, KRH_identification, ID, Names, MH_mass_experimental, MH_mass_MoNA,
           z_experimental, z_MoNA, mz_similarity_score) %>% #, confidence_rank, confidence_source) %>%
    arrange(compound_experimental)


  No.Fuzzy.Match <- setdiff(1:nrow(Experimental.Values),
                            sort(unique(My.Fuzzy.Join$compound_experimental)))
  CL3.Match <- as.integer(unique(My.Fuzzy.Join$compound_experimental))

  # Have any compounds been lost? Check for a TRUE output
  all.experimentals <- sort(c(No.Fuzzy.Match, CL3.Match))
  length(all.experimentals) == length(Experimental.Values$compound_experimental)

  # Make "no match" dataframes for comparison
  No.Fuzzy.Match.df <- Experimental.Values %>%
    filter(compound_experimental %in% No.Fuzzy.Match) %>%
    rename(MH_mass_experimental = MH_mass)


  final <- My.Fuzzy.Join %>%
    bind_rows(No.Fuzzy.Match.df) %>%
    mutate(confidence_rank3 = ifelse(mz_similarity_score > 0.9, 3, NA),
           confidence_source = ifelse(is.na(confidence_rank3), confidence_source, "MoNA"),
           confidence_rank = ifelse(is.na(confidence_rank) & !is.na(confidence_rank3), confidence_rank3, confidence_rank)) %>%
    select(compound_experimental, "MoNA_ID" = ID, "MoNA_Names" = Names, KRH_identification, MH_mass_MoNA, MH_mass_experimental,
           z_MoNA, z_experimental, "mz_similarity_score_MoNa" = mz_similarity_score, confidence_rank, confidence_source)

  Confidence.Level.3 <- final %>%
    left_join(Confidence.Level.2) %>%
    unique() %>%
    select(compound_experimental, KRH_identification, compound_theoretical, ID, MoNA_ID, MoNA_Names, massbank_match, everything()) %>%
    arrange(compound_experimental)

  return(Confidence.Level.3)
}
