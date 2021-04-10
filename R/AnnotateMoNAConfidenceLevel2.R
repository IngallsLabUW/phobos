#' Annotate Confidence Level 2: MoNA
#'
#' @param Confidence.Level.1 A dataframe of experimental values, already annotated for Confidence Level 1.
#' Please refer to the AnnotateConfidenceLevel1 function for details!
#' @param mz.flexibility Flexibility for m/z matching between experimental and theoretical values. Usually defined as 0.02.
#' @param rt.flexibility Flexibility for retention time matching between experimental and theoretical values. Usually defined as ~ 15-30 seconds.
#' @param z Choose positive or negative charge for appropriate value comparison. The options are 1 for positive mode and -1 for negative.
#'
#' @return A complete dataframe, annotated for Confidence Level 1.
#' @export
#'
#' @examples
#' library(phobos)
#' # Load experimental data that has already been annotated for Confidence Level 1 using the AnnotateConfidenceLevel1 function.
#' # Remember that dataframe format is important: please refer to README or the documentation for AnnotateConfidenceLevel1() for details.
#'
#' example_dir <- system.file("example_data", package = "phobos")
#' example_data <- list.files(example_dir, full.names = TRUE)
#' Confidence.Level.1 <- read.csv(grep("Example_ConfidenceLevel1", example_data, value = TRUE))
#' z <- 1
#' example_confidenceL2 <- AnnotateMoNAConfidenceLevel2(Confidence.Level.1 = Confidence.Level.1, z = z,
#' mz.flexibility = 0.02, rt.flexibility = 30)
AnnotateMoNAConfidenceLevel2 <- function(Confidence.Level.1, mz.flexibility, rt.flexibility) {
  # Subtract hydrogen for reference database: this will be changed
  MoNA.Spectra.Neg <- read.csv("data/MoNA_RelationalSpreadsheets/NEG_Spectra.csv") %>%
    dplyr::mutate(MH_mass = M_mass - 1.0072766,
                  z_massbank = -1)
  MoNA.Spectra.Pos <- read.csv("data/MoNA_RelationalSpreadsheets/POS_Spectra.csv") %>%
    dplyr::mutate(MH_mass = M_mass + 1.0072766,
                  z_massbank = 1)

  # Tidy theoretical spectra, dropping NA MS2s
  MoNA.Spectra <- MoNA.Spectra.Neg %>%
    rbind(MoNA.Spectra.Pos) %>%
    dplyr::select("massbank_ID" = ID, Names, spectrum_KRHform_filtered, z_massbank, MH_mass) %>%
    dplyr::mutate_all(., list(~dplyr::na_if(., ""))) %>%
    tidyr::drop_na()

  # Experimental Spectra ----------------------------------------------------
  Experimental.Spectra <- Confidence.Level.1 %>%
    dplyr::filter(!is.na(MS2_experimental)) %>%
    dplyr::select(compound_experimental, z_experimental, MH_mass = "mz_experimental", MS2_experimental) %>%
    unique()

  # Reassign variables from Experimental.Spectra
  mz <- as.numeric(unlist(Experimental.Spectra["MH_mass"]))
  MS2 <- as.character(Experimental.Spectra["MS2_experimental"])
  Mass.Feature <- as.data.frame(Experimental.Spectra["compound_experimental"])

  # Assign variables for paralellized comparison
  experimental.df <- Experimental.Spectra
  potential.candidates <- MoNA.Spectra
  numCores <- parallel::detectCores()
  numCores

  # Compare
  MoNA.Matched <- parallel::mclapply(MoNA.Spectra["MH_mass"], IsolateMoNACandidates,
                                     experimental.df = experimental.df, potential.candidates = potential.candidates,
                                     mc.cores = numCores) %>%
    dplyr::bind_rows() %>%
    dplyr::full_join(Confidence.Level.1) %>%
    dplyr::select(MassFeature, compound_experimental, compound_theoretical, massbank_match, massbank_ID, mz_experimental, mz_theoretical, mz_massbank,
                  rt_sec_experimental, rt_sec_theoretical, column_experimental, column_theoretical, z_experimental, z_theoretical, z_massbank,
                  MS2_experimental, MS2_theoretical, MS2_massbank, ppm_mass_error, massbank_ppm, mz_similarity_score1, rt_similarity_score1,
                  MS2_cosine_similarity1, massbank_cosine_similarity, total_similarity_score, confidence_rank, confidence_source) %>%
    dplyr::arrange(compound_experimental)

  # Combine Confidence Level 2 with Confidence Level 1 ----------------------
  Confidence.Level.2 <- MoNA.Matched %>%
    dplyr::mutate(confidence_rank = ifelse(!is.na(massbank_match),
                                           paste(confidence_rank, "2", sep = "; "), confidence_rank),
                  confidence_source = ifelse(!is.na(massbank_match),
                                             paste(confidence_source, "MoNA", sep = "; "), confidence_source)) %>%
    dplyr::mutate(across(starts_with("confidence"), ~ReplaceNA(.x)))

  return(Confidence.Level.2)
}
