#' Annotate Confidence Level 2: MoNA
#'
#' This function takes a dataframe of experimental values, already annotated for Confidence Level 1, and compares it to theoretical values
#' from MassBank of North America (MoNA). Those compounds that have a good match (see TODO PLACE HERE for specifics on what defines a "good match")
#' to MoNA MS1 & MS2 receive a Confidence Level 2 ranking.
#' In order to run this code, external data from MoNA is required.
#' TODO The code to produce or update the scraped spectra from MassBank is contained in the TODOFUNCTION as part of the phobos package.
#' @importFrom magrittr "%>%"
#'
#' @param Confidence.Level.1 A dataframe of experimental values, already annotated for Confidence Level 1.
#' Please refer to the AnnotateConfidenceLevel1 function for details!
#' @param mz.flexibility Flexibility in Daltons for m/z matching between experimental and theoretical values. Usually defined as 0.02.
#' @param MassBank.Neg Spectra in negative ion mode from MassBank of North America, scraped and downloaded. This csv is available on the Ingalls Shared Drive in the MARS_Project folder, titled NEG_Spectra.csv
#' @param MassBank.Pos Spectra in positive ion mode from MassBank of North America, scraped and downloaded. This csv is available on the Ingalls Shared Drive in the MARS_Project folder, titled POS_Spectra.csv
#' @param rt.flexibility Flexibility for retention time matching between experimental and theoretical values. Usually defined as ~ 15-30 seconds.
#'
#' @return A dataframe annotated for Confidence Level 2, in addition to the Confidence Level 1 annotation from the previous step.
#' @export
#'
#' @examples
#' library(phobos)
#' # Load experimental data that has already been annotated for Confidence Level 1 using the AnnotateConfidenceLevel1 function.
#' # Remember that dataframe format is important: please refer to README or the documentation for AnnotateConfidenceLevel1() for details.
#'
#' Confidence.Level.1 <- read.csv("example_data/Example_ConfidenceLevel1.csv")
#' MassBank.Neg <- read.csv("your/file/path/here/NEG_Spectra.csv")
#' MassBank.Pos <- read.csv("your/file/path/here/POS_Spectra.csv")
#' Example_ConfidenceLevel2 <- AnnotateMoNAConfidenceLevel2(Confidence.Level.1 = Confidence.Level.1, MassBank.Neg = MassBank.Neg,
#' MassBank.Pos = MassBank.Pos, mz.flexibility = 0.02, rt.flexibility = 30)
AnnotateMoNAConfidenceLevel2 <- function(Confidence.Level.1, MassBank.Neg, MassBank.Pos, mz.flexibility, rt.flexibility) {

  # Subtract hydrogen for reference database
  MoNA.Spectra.Neg <- MassBank.Neg %>%
    dplyr::mutate(MH_mass = M_mass - 1.0072766,
                  z_massbank2 = -1)
  MoNA.Spectra.Pos <- MassBank.Pos %>%
    dplyr::mutate(MH_mass = M_mass + 1.0072766,
                  z_massbank2 = 1)

  # Tidy theoretical spectra, dropping NA MS2s
  MoNA.Spectra <- MoNA.Spectra.Neg %>%
    rbind(MoNA.Spectra.Pos) %>%
    dplyr::mutate(ID = paste("ID:", ID)) %>%
    dplyr::select("massbank_ID" = ID, Names, spectrum_KRHform_filtered, z_massbank2, MH_mass) %>%
    dplyr::mutate_all(., list(~dplyr::na_if(., ""))) %>%
    tidyr::drop_na()

  # Experimental Spectra ----------------------------------------------------
  Experimental.Spectra <- Confidence.Level.1 %>%
    dplyr::filter(!is.na(MS2_experimental)) %>%
    dplyr::select(primary_key, z_experimental, MH_mass = "mz_experimental", MS2_experimental) %>%
    unique()

  # Reassign variables from Experimental.Spectra
  mz <- as.numeric(unlist(Experimental.Spectra["MH_mass"]))
  MS2 <- as.character(Experimental.Spectra["MS2_experimental"])
  Mass.Feature <- as.data.frame(Experimental.Spectra["primary_key"])

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
    tidyr::unite(massbank_match2, c(massbank_ID, massbank_match), sep = "; ", na.rm = TRUE, remove = FALSE) %>%
    dplyr::select(MassFeature, primary_key, compound_theoretical, massbank_match2, mz_experimental, mz_theoretical, mz_massbank2,
                  rt_sec_experimental, rt_sec_theoretical, column_experimental, column_theoretical, z_experimental, z_theoretical, z_massbank2,
                  MS2_experimental, MS2_theoretical, MS2_massbank, ppm_mass_error1, massbank_ppm, mz_similarity_score1, rt_similarity_score1,
                  MS2_cosine_similarity1, massbank_cosine_similarity2, total_similarity_score1, confidence_rank, confidence_source) %>%
    dplyr::arrange(primary_key)

  # Combine Confidence Level 2 with Confidence Level 1 ----------------------
  Confidence.Level.2 <- MoNA.Matched %>%
    dplyr::mutate(confidence_rank = ifelse(!is.na(massbank_match2),
                                           paste(confidence_rank, "2", sep = "; "), confidence_rank),
                  confidence_source = ifelse(!is.na(massbank_match2),
                                             paste(confidence_source, "MassBank", sep = "; "), confidence_source)) %>%
    dplyr::mutate(across(starts_with(c("confidence", "massbank")), ~ReplaceNA(.x))) %>%
    unique()
  Confidence.Level.2[Confidence.Level.2 == ""] <- NA

  return(Confidence.Level.2)
}
