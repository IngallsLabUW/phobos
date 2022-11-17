#' Calculate value similarity score.
#'
#' The equation for calculating similarity scores is as follows:
#' \deqn{exp(0.5 x (experimental.value - theoretical.value / flexibility)^2)}
#'
#' Similarity scores between experimental and theoretical values are calculated using the above function.
#' The background hypothesis of the equations for accurate mass and retention time similarity is
#' that the differences between experimental and theoretical values follow a Gaussian distribution.
#' For more information, please refer to "MS-DIAL: Data Independing MS/MS Deconvolution for Comprehensive
#' Metabolome Analysis" Nat. Methods. 2015 June. doi: 10.1038/nmeth.3393.
#'
#'
#' @param experimental.value The experimental value to be compared to theoretical.
#' @param theoretical.value The theoretical value for comparison.
#' @param flexibility User-defined search tolerance, also used as the standard deviation value.
#'
#' @return
CalculateSimilarityScore <- function(experimental.value, theoretical.value, flexibility) {
  similarity.score = exp(-0.5 * (((experimental.value - theoretical.value) / flexibility) ^ 2))

  return(similarity.score)
}


#' Check known contaminants
#'
#' @param experimental.df Experimental dataframe, containing at least the columns primary_key (numeric), MassFeature (character),
#' mz_experimental (numeric), and z_experimental (numeric).
#' @param common.contaminants A dataframe listing known contaminants for comparison. This csv can be found on the Shared Ingalls Google Drive
#' in the MARS_project folder.
#' @param ppm.tolerance Acceptable level of parts per million (ppm). Default is 15.
#'
#' @return matchedKnownContaminants, a small dataframe containing potential contaminants. These can be matched to the
#' annotated dataframe by the primary_key column.
#' @export
#'
checkContaminants <- function(experimental.df, common.contaminants, ppm.tolerance) {
  if(missing(ppm.tolerance)) {
    ppm.tolerance <- 15
  }
  df.to.join <- experimental.df %>%
    dplyr::select(primary_key, MassFeature, mz_experimental, column_experimental, z_experimental) %>%
    dplyr::rename(mz = mz_experimental,
                  z = z_experimental) %>%
    unique()

  knownContaminants <- common.contaminants

  matchedKnownContaminants <- difference_inner_join(x = knownContaminants, y = df.to.join,
                                                         by = "mz", max_dist = 0.02, distance_col = NULL) %>%
    dplyr::rename_with(., ~gsub("\\.x", "_contam", .x)) %>%
    dplyr::rename_with(., ~gsub("\\.y", "_experimental", .x)) %>%
    dplyr::mutate(ppm = (abs(mz_contam-mz_experimental )/mz_contam *10^6)) %>%
    dplyr::filter(ppm < ppm.tolerance, z_contam == z_experimental) %>%
    dplyr::select(primary_key, MassFeature, Compound, Origin, Ion.Type, Formula, ppm)

  return(matchedKnownContaminants)
}

#' ConcatToScan
#'
#' @param concat.format MS2 data in "concatenated format": mz and intensity separated by comma,
#' and separated from the next set of fragments by a semicolon. For example, "148.0, 100; 147.0, 50;", etc.
#'
#' @return scantable: The same MS2 data separated into a mini dataframe of mz and intensity.
#' @export
ConcatToScan <- function(concat.format) {
  scantable <- cSplit(concat.format, "MS2", sep = ";") %>%
    pivot_longer(cols = starts_with("MS2")) %>%
    separate(value, sep = ",", into = c("mz", "intensity")) %>%
    rename(fragment = name)

  return(scantable)
}

#' Compare experimental mass features to scraped MoNA data, based on MS1 and MS2 information.
#'
#' @param MoNA.Mass Single MoNA mass, isolated from scraped MoNA df.
#' @param experimental.df experimental dataframe, islated to contain only MS1 and MS2 data.
#'
#' @return final.candidates: dataframe of experimental data, matched with information from MoNA.
IsolateMoNACandidates <- function(MoNA.Mass, experimental.df, potential.candidates) {
  requireNamespace("dplyr", quietly = TRUE)
  potential.candidates <- potential.candidates %>%
    # dplyr::filter(MH_mass > MoNA.Mass - 0.020, # These should not be hard-coded! Also not really sure why this is here.
    #        MH_mass < MoNA.Mass + 0.020) %>%
    fuzzyjoin::difference_inner_join(experimental.df, by = "MH_mass", max_dist = 0.02) %>%
    dplyr::filter(z_massbank2 == z_experimental) %>%
    dplyr::rename(scan1 = spectrum_KRHform_filtered, # scan1 is MS2 from MoNA
                  scan2 = MS2_experimental,          # scan2 is MS2 from the experimental data
                  mass1 = MH_mass.x,                 # mass1 is the mass from MoNA
                  mass2 = MH_mass.y)                 # mass2 is the mass from experimental data

  if (length(potential.candidates$massbank_ID) == 0) {
    print("There are no potential candidates.")
    No.Match.Return <- experimental.df %>%
      dplyr::mutate(massbank_match = NA,
                    massbank_ppm = NA,
                    MS2_cosine_similarity2 = NA)

    return(No.Match.Return)
  }

  # Add cosine similarity scores
  print("Making potential candidates")

  potential.candidates$MS2_cosine_similarity2 <- apply(potential.candidates, 1, FUN = function(x) MakeMS2CosineDataframe(x))

  final.candidates <- potential.candidates %>%
    dplyr::mutate(massbank_match = paste(Names, massbank_ID, sep = " ID:"),
                  massbank_ppm = abs(mass2 - mass1) / mass1 * 10^6) %>%
    dplyr::rename(MS2_massbank = scan1,
                  mz_massbank2 = mass1,
                  MS2_experimental = scan2,
                  mz_experimental = mass2) %>%
    unique() %>%
    dplyr::filter(massbank_ppm < 7,
                  MS2_cosine_similarity2 > 0.5) %>%
    dplyr::arrange(desc(MS2_cosine_similarity2))

  return(final.candidates)
}

#' Create a filtered mini dataframe (a scantable) from two columns in a user-defined dataframe.
#'
#' @param df A dataframe containing "scan1" and "scan2" columns.
#'
#' @return
MakeMS2CosineDataframe <- function(df) {
  scan1 <- MakeScantable(df["scan1"])
  scan2 <- MakeScantable(df["scan2"])

  weight1 <- (scan1[, 1] ^ 2) * sqrt(scan1[, 2])
  weight2 <- (scan2[, 1] ^ 2) * sqrt(scan2[, 2])

  diff.matrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x)
  diff.matrix <- as.matrix(diff.matrix)
  same.index <- which(abs(diff.matrix) < 0.02, arr.ind = TRUE)
  cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
    (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))

  return(cosine.similarity)
}

#' Create a filtered mini dataframe from a concatenated scanlist of MS2s.
#'
#' @param concatenated.scan Single observation from a dataframe containing MS2 m/z and intensity spectra, separated by semicolons.
#'
#' @return scantable: A tiny dataframe, containing columns of mz and intensity.
#' Intensity is scaled to 100 and filtered to drop all intensity values below 0.5.
#'
MakeScantable <- function(concatenated.scan) {
  requireNamespace("dplyr", quietly = TRUE)

  scantable <- read.table(text = as.character(concatenated.scan),
                          col.names = c("mz", "intensity"), fill = TRUE) %>%
    dplyr::mutate(mz = as.numeric(mz %>% stringr::str_replace(",", "")),
                  intensity = as.numeric(intensity %>% stringr::str_replace(";", "")),
                  intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
    dplyr::filter(intensity > 0.5) %>%
    dplyr::arrange(desc(intensity))

  return(scantable)
}

#' Finds the weighted cosine similarity between two sets of MS2 spectra.
#'
#' @param scan1 Tiny dataframe of MS2, from first set of values to be compared. Column 1 is mz, column 2 is intensity. Dataframe is the output of the MakeScantable() function.
#' @param scan2 Tiny dataframe of MS2, from second set of values to be compared. Column 1 is mz, column 2 is intensity. Dataframe is the output of the MakeScantable() function.
#' @param mz.flexibility User-defined search tolerance for mz matching. This value should be consistent throughout the analysis.
#'
#' @return cosine.similarity: A weighted similarity score between 0 and 1, indicating the cosine relationship of the two vectors.
#'
MS2CosineSimilarity <- function(scan1, scan2, mz.flexibility) {

  mz.flexibility <- mz.flexibility
  scan1 <- MakeScantable(scan1)
  scan2 <- MakeScantable(scan2)

  weight1 <- (scan1[, "mz"] ^ 2) * sqrt(scan1[, "intensity"])
  weight2 <- (scan2[, "mz"] ^ 2) * sqrt(scan2[, "intensity"])

  diff.matrix <- sapply(scan1[, "mz"], function(x) scan2[, "mz"] - x)
  same.index <- which(abs(diff.matrix) < 0.02, arr.ind = TRUE)
  cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
    (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))

  return(cosine.similarity)
}

#' Paste items together, ignoring NAs
#'
#' @param ... List of items to be pasted
#' @param sep Defined separator of ";" to remain consistent with MARS concatenation protocols
#'
#' @return Pasted item, ignoring the character "NA"s produced by regular paste.
#' @export
#'
pasteWithoutNA <- function(..., sep = "; ") {
  my.list <- list(...)
  my.list <- lapply(my.list, function(x) {x[is.na(x)] <- ""; x})
  product <- gsub(paste0("(^", sep, "|", sep, "$)"), "",
                  gsub(paste0(sep, sep), sep,
                       do.call(paste, c(my.list, list(sep = sep)))))
  is.na(product) <- product == ""
  product
}


#' Remove character values of "NA;" found in MARS outputs.
#'
#' @param column Character column that contains one or more "NA; " values.
#'
#' @return
ReplaceNA <- function(column) {
  gsub("NA; ", "", column)
}

#' Scale Concatenated MS2 data
#'
#' @param concatenated.scan MS2 data in concatenated format: paired mz and intensity values
#' separated by a semicolon.
#'
#' @return scantable: MS2 data scaled so the highest intensity corresponds to 100, and all
#' intensity values below 0.5 are thrown out.
#' @export
#'
ScaleMS2 <- function(concatenated.scan) {
  scantable <- read.table(text = as.character(concatenated.scan),
                          col.names = c("mz", "intensity"), fill = TRUE) %>%
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", "")),
           intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
    filter(intensity > 0.5) %>%
    arrange(desc(intensity)) %>%
    unite("MS2", mz:intensity, sep = ", ") %>%
    mutate(MS2 = paste(MS2, collapse = "; ")) %>%
    unique()

  return(scantable)
}

#' ScanToConcat
#'
#' @param scantable.format MS2 data in scantable format: a mini dataframe of mz and intensity.
#'
#' @return concatenated: the same MS2 data converted to a concatenated format of paured mz and intensity
#' values separated by semicolons.
#' @export
ScanToConcat <- function(scantable.format) {
  concatenated <- scantable.format %>%
    unique() %>%
    group_by(compound_name, filename) %>%
    unite("MS2", mz:intensity, sep = ",") %>%
    ungroup() %>%
    group_by(compound_name, voltage) %>%
    mutate(MS2 = paste(MS2, collapse = "; ")) %>%
    select(-fragment) %>%
    unique() %>%
    as.data.frame()

  return(concatenated)
}
