% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ScaleMS2}
\alias{ScaleMS2}
\title{Scale Concatenated MS2 data}
\usage{
ScaleMS2(concatenated.scan)
}
\arguments{
\item{concatenated.scan}{MS2 data in concatenated format: paired mz and intensity values
separated by a semicolon.}
}
\value{
scantable: MS2 data scaled so the highest intensity corresponds to 100, and all
intensity values below 0.5 are thrown out.
}
\description{
Scale Concatenated MS2 data
}
