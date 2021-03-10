MakeScantable <- function(scan) {
  # Create a filtered MS2 scantable from a concatenated scanlist of MS2s.
  #
  # Args
  #   scan: Single observation from a dataframe containing MS2 mz and intensity
  #         spectra, separated by semicolons (;).
  #
  # Returns
  #   scantable: Tiny dataframe, containing columns of mz and intensity.
  #              Intensity is scaled to 100 and filtered to drop all intensity
  #              values below 0.5.
  #
  scantable <- read.table(text = as.character(scan),
                          col.names = c("mz", "intensity"), fill = TRUE) %>%
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", "")),
           intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
    filter(intensity > 0.5) %>%
    arrange(desc(intensity))

  return(scantable)
}
