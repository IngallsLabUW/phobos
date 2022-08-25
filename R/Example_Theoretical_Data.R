#' Liquid Chromatography with Tandem Mass Spectrometry data
#'
#' A dataset containing a subset of LCMS information, including mass to charge ratio, retention time, column, charge, and MS2 data for a range of organic compounds.
#'
#' @format A data frame with 1632 rows and 7 variables:
#' \describe{
#'   \item{compound}{compound name}
#'   \item{mz}{mass to charge ratio}
#'   \item{rt}{retention time, seconds}
#'   \item{column}{mass spectrometry column time}
#'   \item{z}{mode, either positive (1) or negative (-1)}
#'   \item{file_name}{unique identifier of raw data}
#'   \item{MS2}{MS2 data, concatenated to mz, intensity;}
#' }
#' @source \url{https://github.com/IngallsLabUW/Ingalls_Standards}

