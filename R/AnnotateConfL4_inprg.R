library(pbapply)
library(tidyverse)

### Annotate Confidence Level 4: Everything else!

# Notes -------------------------------------------------------------------
# Should we have some kind of filter for level 4? As in,
#   this is definitely a consistent peak, shows up? Or does that come ahead of time?

# Outline -------------------------------------------------------------------

# Prepare all data -------------------------------------------------------------------

AnnotateConfidenceLevel4 <- function(Confidence.Level.3) {

  Confidence.Level.4 <- Confidence.Level.3 %>%
    dplyr::mutate(confidence_rank = ifelse(is.na(confidence_rank), 4, confidence_rank)) %>%
    dplyr::arrange(primary_key)

  return(Confidence.Level.4)
}

