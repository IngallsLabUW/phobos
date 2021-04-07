#' Annotate Confidence Level 4: Everything else!
#'
#' @param Confidence.Level.3 A dataframe annotated for Confidence Levels 1, 2, & 3.
#' Please refer to the README for details.
#'
#' @return A fully annotated dataframe, annotated for Confidence Levels 1, 2, 3, & 4.
#' @export
#'
#' @examples
#' #' library(phobos)
#' # Load experimental data that has already been annotated for Confidence Levels 1, 2, & 3.
#'
#' example_dir <- system.file("example_data", package = "phobos")
#' example_data <- list.files(example_dir, full.names = TRUE)
#' Confidence.Level.3 <- read.csv(grep("Example_ConfidenceLevel3", example_data, value = TRUE))
#'
#' example_confidenceL4 <- AnnotateConfidenceLevel4(Confidence.Level.3 = Confidence.Level.3)
AnnotateConfidenceLevel4 <- function(Confidence.Level.3) {

  Confidence.Level.4 <- Confidence.Level.3 %>%
    mutate(confidence_rank = ifelse(is.na(confidence_rank), 4, confidence_rank))

  return(Confidence.Level.4)
}

