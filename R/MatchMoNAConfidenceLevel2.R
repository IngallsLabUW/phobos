MatchMoNAConfidenceLevel2 <- function(Confidence.Level.1, mz.flexibility, rt.flexibility, ion.mode) {
  # cosine.score.cutoff <- 0.5
  # massbank.ppm.cutoff <- 5
  # mz.flexibility <- 0.02
  # rt.flexibility <- 30 # seconds
  # ion.mode <- "negative"

  # Subtract hydrogen for reference database
  if (ion.mode == "negative") {
    MoNA.Spectra <- read.csv("test_data/MoNA_RelationalSpreadsheets/NEG_Spectra.csv") %>%
      dplyr::mutate(MH_mass = M_mass - 1.0072766)
  } else {
    MoNA.Spectra <- read.csv("test_data/MoNA_RelationalSpreadsheets/POS_Spectra.csv") %>%
      dplyr::mutate(MH_mass = M_mass + 1.0072766)
  }

  # Tidy theoretical spectra, dropping NA MS2s
  MoNA.Spectra <- MoNA.Spectra %>%
    dplyr::select(ID, Names, spectrum_KRHform_filtered, MH_mass) %>%
    dplyr::mutate_all(., list(~dplyr::na_if(.,""))) %>%
    tidyr::drop_na()

  # Experimental Spectra ----------------------------------------------------
  Confidence.Level.1 <- Confidence.Level.1
  Experimental.Spectra <- Confidence.Level.1 %>%
    dplyr::filter(!is.na(MS2_experimental),
                  z_experimental == 1) %>%
    dplyr::select(compound_experimental, MH_mass = "mz_experimental", MS2_experimental) %>%
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
                                     experimental.df = experimental.df, potential.candidates = potential.candidates, mc.cores = numCores) %>%
    dplyr::bind_rows() %>%
    dplyr::full_join(Confidence.Level.1) %>%
    dplyr::select(compound_experimental, KRH_identification, compound_theoretical, massbank_match, ID, mz_experimental, mz_theoretical, mz_massbank,
                  rt_sec_experimental, rt_sec_theoretical, column_experimental, column_theoretical, z_experimental, z_theoretical,
                  MS2_experimental, MS2_theoretical, MS2_massbank, ppm_mass_error, massbank_ppm, mz_similarity_score, rt_similarity_score,
                  MS2_cosine_similarity, total_similarity_score, massbank_cosine_similarity, confidence_rank, confidence_source) %>%
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
