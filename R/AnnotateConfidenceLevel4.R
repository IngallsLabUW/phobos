#' Annotate Confidence Level 4: Everything else!
#'
#' @param Confidence.Level.3 A dataframe annotated for Confidence Levels 1, 2, & 3. Level 3 must have
#' received annotation from both MassBank and KEGG! Please refer to the README for details.
#'
#' @return A fully annotated dataframe, annotated for Confidence Levels 1, 2, 3, & 4.
#' @export
#'
#' @examples
#' #' library(phobos)
#' # Load experimental data that has already been annotated for Confidence Levels 1, 2, & 3.
#
#' Confidence.Level.3_KEGG <- read.csv("example_data/Example_ConfidenceLevel3_KEGG.csv")
#'
#' example_confidenceL4 <- AnnotateConfidenceLevel4(Confidence.Level.3 = Confidence.Level.3_KEGG)
AnnotateConfidenceLevel4 <- function(Confidence.Level.3) {

  Confidence.Level.4 <- Confidence.Level.3 %>%
    dplyr::mutate(confidence_rank = ifelse(is.na(confidence_rank), 4, confidence_rank)) %>%
    dplyr::arrange(primary_key)
  # We then immediately sort by primary key in the demo - this arrange() can
  # probably be dropped.

  return(Confidence.Level.4)
}

