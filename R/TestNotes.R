### WKumler + RML Startup script

## Walkthrough Layout

# Begin with experimental dataframe

pmppm <- function(mass, ppm=4)c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))

experimental.values <- read.csv("example_data/Example_Experimental_Data.csv")

# pass experimental df to CalcTotSim1, which takes each observation in each column as an argument, rt, mz, ms2str, col, z)
# Takes the theoretical value and ensures we are matching within mz window, z, and column.
# Calculates similarity scores for rt, mz, ms2.
CalculateTotalSimScore1 <- function(mz_i, rt_i, col_i, z_i, MS2str_i, theoretical_db) {
  theoretical_db %>%
    filter(db_mz < mz_i + ((mz_i * 50)/1e6) & db_mz > mz_i - ((mz_i * 50)/1e6)) %>% ## 50 is ppm error, maybe should be an argument?
    filter(col_db == col_i) %>%
    filter(z_db == z_i) %>%
    mutate(MS1SimScore = CalculateMS1SimScore(mz_i, mz_db, flex = 5)) %>% # 5ppm? argument?
    mutate(RT1SimScore = CalculateRT1SimScore(rt_i, rt_db, flex = 30)) %>% # Ditto
    mutate(MS21SimScore = CalcualteMS21SimScore(MS2str_i, MS2str_db, flex = 0.02)) %>% # This flex is totally made up
    mutate(TotalSimScore = CalculateTotalSimScore(MS1SimScore, RT1SimScore, MS2SimScore)) %>%
    select(name_db, ends_with("SimScore"))

  return(theoretical_db)
}

# Outputs dataframe of original unknown feature, with accompanying columns of all sim scores, and all matches.

experimental.df %>%
  mutate(TotalSimScoreDF = CalculateTotalSimScore1(mz, rt, col, z, MS2str, theoretical)) %>%
  mutate(Cl1_choice = sapply(AnnotateCL1(TotalSimScoreDF)))


## Annotate CL1 function needs to be written, actually makes the confidence rank decision and source (aka theoretical_db). Wishlist: timestamp
#dplyr::mutate(confidence_rank = ifelse(mz_similarity_score1 > 0.9 & rt_similarity_score1 > 0.75 & ppm_mass_error1 < 7, 1, NA),
#               confidence_source = ifelse(!is.na(confidence_rank), "Ingalls_Standards", NA)) %>%
