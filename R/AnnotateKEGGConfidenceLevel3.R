#' Annotate Confidence Level 3: KEGG
#'
#' @param Confidence.Level.3.MoNA A dataframe, already annotated for confidence levels 1, 2 (from Massbank) and 3 (from Massbank).
#' Please refer to the README or documentation for the AnnotateConfidenceLevel1, AnnotateMoNAConfidenceLevel2, and AnnotateMoNAConfidenceLevel3 functions for details.
#' @param KEGG.Data Compound information from the KEGG database, scraped and downloaded. This csv is available on the Ingalls Shared Drive in the MARS_Project folder, titled KEGG_Data.csv
#' @param mz.flexibility Flexibility for m/z matching between experimental and theoretical values. Usually defined as 0.02.
#' @param ppm.tolerance Parts per million tolerance for matching between experimental and KEGG data. Usually defined as 15.
#'
#' @return A dataframe annotated for Confidence Levels 1, 2, and 3 using the Ingalls standards, MassBank of North America, and KEGG.
#' @export
#'
#' @examples
#' library(phobos)
#' # Load experimental data that has already been annotated for Confidence Levels 1, 2 & 3
#' # using the AnnotateConfidenceLevel1, AnnotateMoNAConfidenceLevel2 and AnnotateMoNAConfidenceLevel3 functions.
#'
#' Confidence.Level.3.MoNA <- read.csv("example_data/Example_ConfidenceLevel3MoNA.csv")
#' KEGG.Data <- read.csv("your/file/path/here/KEGG_Data.csv")
#'
#' Example_KEGG_ConfidenceLevel3 <- AnnotateKEGGConfidenceLevel3(Confidence.Level.3.MoNA = Confidence.Level.3.MoNA,
#' KEGG.Data <- KEGG_Data, mz.flexibility = 0.02, ppm.tolerance = 15)
AnnotateKEGGConfidenceLevel3 <- function(Confidence.Level.3.MoNA, KEGG.Data, mz.flexibility, ppm.tolerance) {

  KEGG.Data <- KEGG.Data %>%
    rename(Compound_KEGG = OtherCmpds)

  ToJoin <- Confidence.Level.3.MoNA %>%
    rename(mz = mz_experimental) %>%
    select(MassFeature, primary_key, mz, column_experimental, z_experimental) %>%
    unique()

  # Isolate pos and neg theoretical data
  keggPos <- KEGG.Data %>%
    select(Compound_KEGG, PosMZ) %>%
    mutate(column_KEGG = "HILIC",
           z = 1) %>%
    rename(mz = PosMZ)

  keggNeg <- KEGG.Data %>%
    select(Compound_KEGG, NegMZ) %>%
    mutate(column_KEGG = "HILIC",
           z = -1) %>%
    rename(mz = NegMZ)

  # Combine to full theoretical dataset
  keggCompounds <- keggNeg %>%
    rbind(keggPos)

  # Fuzzyjoin datasets
  My.Fuzzy.Join <- difference_inner_join(x = keggCompounds, y = ToJoin,
                                         by = "mz", max_dist = mz.flexibility, distance_col = NULL) %>%
    rename(mz_KEGG = mz.x,
           mz_experimental = mz.y,
           z_KEGG = z) %>%
    mutate(ppm = (abs(mz_KEGG - mz_experimental) / mz_KEGG *10^6),
           mz_similarity_scoreKEGG = exp(-0.5 * (((mz_experimental - mz_KEGG) / 0.02) ^ 2))) %>%
    filter(ppm < ppm.tolerance,
           z_KEGG == z_experimental) %>%
    unique() %>%
    rename(KEGGppm = ppm) %>%
    arrange(primary_key)

  # Match KEGG names to IDs
  matchedNames <- left_join(My.Fuzzy.Join, KEGG.Data) %>%
    select(Compound_KEGG, Name) %>%
    group_by(Compound_KEGG) %>%
    summarise(KEGGMatchesNames = as.character(paste(Name,  collapse=" ")))

  # Combine to full matched df
  final <- My.Fuzzy.Join %>%
    left_join(matchedNames) %>%
    select(MassFeature, primary_key, Compound_KEGG, mz_KEGG, z_KEGG, KEGGMatchesNames, KEGGppm, mz_similarity_scoreKEGG)
  Confidence.Level.3.MoNA[Confidence.Level.3.MoNA == ""] <- NA

  Confidence.Level.3 <- final %>%
    full_join(Confidence.Level.3.MoNA) %>%
    dplyr::arrange(primary_key) %>%
    unique() %>%
    dplyr::mutate(confidence_rank3 = ifelse(mz_similarity_scoreKEGG > 0.9, 3, NA)) %>%
    dplyr::mutate(confidence_rank = ifelse(!is.na(confidence_rank) & !is.na(confidence_rank3), paste(confidence_rank, "3", sep = "; "),
                                           ifelse(!is.na(confidence_rank3), confidence_rank3, confidence_rank)))  %>%
    dplyr::mutate(confidence_source = ifelse(stringr::str_detect(confidence_rank, "3"),
                                             apply(cbind(confidence_source, "KEGG"), 1, function(x) paste(x[!is.na(x)], collapse = "; ")),
                                             confidence_source)) %>%
    dplyr::select(primary_key, MassFeature, compound_theoretical, massbank_match2, massbank_match3, Compound_KEGG, KEGGMatchesNames,
                  mz_experimental:mz_massbank2, MH_mass_experimental, MH_mass_MoNA, mz_KEGG,
                  z_experimental, z_theoretical, z_massbank2, z_massbank3, z_KEGG, rt_sec_experimental:column_theoretical,
                  MS2_experimental:MS2_massbank, ppm_mass_error1, massbank_ppm, KEGGppm, mz_similarity_score1, mz_similarity_score3,
                  rt_similarity_score1:total_similarity_score1, everything(), -confidence_rank3)

  return(Confidence.Level.3)

}
