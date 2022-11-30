library(tidyverse)

### WKumler + RML Startup script

## MAIN QUESTIONS FOR WILL
# Clustering: Will do cosine similarity between "yellow triangle/consensus msms spectra" and experimental dot
# MS2 Sim Score is not working like the others in terms of producing a column upon comparison to the experimental values.
# MS2 Sim Score separate_rows also produced some headaches.
# Output looks ok but haven't completed the final step of modifying the original dataframe (?) input or actually selecting the best value.
## Need to make an adjustment for "n consensus files", not 4:6/9:11 for the testing purposes. This could be a little softer...
## increasing suspicion as we move away from 5. How many "intensity clusters" at mz clusters do we have? Now we're accounting
## for two intensity clusters, but it's possible we'd have three, four...

experimental.values <- read.csv("example_data/Example_Experimental_Data.csv")
theoretical.values <- read.csv("example_data/Example_Theoretical_Data.csv")

# Pass experimental values to the CalculateTotalSimScore1 function, which takes each observation in each column as an argument (mz, rt, col, z, ms2str)
# Compares to the complete set of theoretical value and ensures we are matching within mz window, z, and column.
# Calculates similarity scores for rt, mz, ms2.

CalculateTotalSimScore1 <- function(mz_i, rt_i, col_i, z_i, MS2str_i, ppm_error, theoretical_db) {
  output <- theoretical_db %>%
    filter(mz < mz_i + ((mz_i * ppm_error)/1e6) & mz > mz_i - ((mz_i * ppm_error)/1e6)) %>% ## 50 is ppm error
    filter(column == col_i) %>%
    filter(z == z_i) %>%
    mutate(MS1SimScore = MS1SimilarityScore(mz_exp = mz_i, mz_theo = mz, flex = 5)) %>%
    mutate(RT1SimScore = RTSimilarityScore(rt_exp = rt_i, rt_theo = rt, flex = 30)) %>%
    # TODO MS2 flex is totally made up, and is still the old function, and currently isn't working like the Ms1 and Rt sim scores (won't apply to full df)
    mutate(MS21SimScore = MS21SimilarityScore(ms2_exp = MS2str_i, ms2_theo = MS2[1], flex = 0.02)) %>%
    mutate(TotalSimScore = TotalSimilarityScore(MS1SimScore, RT1SimScore, MS21SimScore)) %>%
    select(compound, ends_with("SimScore"))

  return(output)
}

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
} ## TODO this is still the old function- the separate_rows produces an error

MS1SimilarityScore <- function(mz_exp, mz_theo, flex) {
  print(mz_theo)
  similarity.score = exp(-0.5 * (((mz_exp - mz_theo) / flex) ^ 2))

  return(similarity.score)
}

MS21SimilarityScore <- function(ms2_exp, ms2_theo, flex) {
# Comparison will be between experimental and "consensus" spectra according to Horai et. al 2010 as justification for using spectra
  scan1 <- MakeScantable(ms2_exp)
  scan2 <- MakeScantable(ms2_theo)

  weight1 <- (scan1[, "mz"] ^ 2) * (scan1[, "intensity"] ^ 0.5) # Need consensus to assign weights
  weight2 <- (scan2[, "mz"] ^ 2) * (scan2[, "intensity"] ^ 0.5)

  diff.matrix <- sapply(scan1[, "mz"], function(x) scan2[, "mz"] - x)
  same.index <- which(abs(diff.matrix) < flex, arr.ind = TRUE)
  cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
    (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))

  return(cosine.similarity)
}

RTSimilarityScore <- function(rt_exp, rt_theo, flex) {
  print(rt_theo)
  similarity.score = exp(-0.5 * (((rt_exp - rt_theo) / flex) ^ 2))

  return(similarity.score)
}

TotalSimilarityScore <- function(ms1_sim, rt_sim, ms2_sim) {
  total.similarity.score <- ((ms1_sim + rt_sim + ms2_sim) / 3) * 100

  return(total.similarity.score)
}


## Produces dataframe of potential matches and all sim scores.
OutputDF <- CalculateTotalSimScore1(mz_i = experimental.values[1, 2], rt_i = experimental.values[1, 3],
                                col_i = experimental.values[1, 4], z_i = experimental.values[1, 5],
                                MS2str_i = experimental.values[1, 6], ppm_error = 100000, theoretical_db = theoretical.values)

# Outputs dataframe of

experimental.values %>%
  mutate(TotalSimScoreDF = lapply(CalculateTotalSimScore1(mz, rt, col, z, MS2str, theoretical.values))) %>%
  ## Annotate CL1 function needs to be written, actually makes the confidence rank decision and source (aka theoretical_db).
  mutate(Cl1_choice = sapply(AnnotateCL1(TotalSimScoreDF))) # will append single string oclumn
#appends a "best match" column


# Wishlist: timestamp
