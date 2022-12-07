library(tidyverse)

## Theoretical data created from four standard runs combined.
## "Experimental data" created from the fifth standard run.
FourRunsTheoretical <- read_csv("example_data/Ingalls_Lab_Standards_MSMS.csv") %>%
  filter(!str_detect(filename, "pos5|neg5"))
SingleRunExperimental <- read_csv("example_data/Ingalls_Lab_Standards_MSMS.csv") %>%
  filter(str_detect(filename, "pos5|neg5"))



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
MS2SimilarityScore <- function(ms2_exp, ms2_theo, flex) {
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
mz_group <- function(mz_vals, ppm) {
  ## Will's mz grouping function: exclusively groups on mz
  group_vec <- numeric(length(mz_vals))
  group_num <- 1L
  init_vec <- mz_vals
  names(init_vec) <- seq_along(init_vec)
  while(length(init_vec)>0){
    mz_i <- init_vec[1]
    err <- mz_i*ppm/1000000
    mz_idxs <- init_vec>mz_i-err & init_vec<mz_i+err # locate indexes falling within defined error range
    group_vec[as.numeric(names(mz_idxs)[mz_idxs])] <- group_num
    init_vec <- init_vec[!mz_idxs]
    group_num <- group_num+1L
  }
  group_vec
}

testimscore <- MS21SimilarityScore(SingleRunExperimental[1,2], FourRunsTheoretical[1,2], 0.02)




# Test the grouping function on a single compound
dat <- init_dat %>%
  filter(compound_name == "L-Alanine") %>%
  separate_rows(MS2, sep = "; ") %>%
  separate(MS2, into = c("mz", "int"), sep = ", ") %>%
  mutate(mz=as.numeric(mz)) %>%
  mutate(int=as.numeric(int)) %>%
  group_by(filename, voltage) %>%
  mutate(int=int/max(int)*100) %>%
  ungroup() %>%
  arrange(desc(int)) %>%
  mutate(mz_group=mz_group(mz, ppm = 10)) %>%
  filter(mz_group%in%sort(unique(mz_group))[1:20]) %>% # Take only the first 20 mz groups
  ggplot() +
  geom_point(aes(x=mz, y=int, color=factor(mz_group), group=filename)) +
  facet_wrap(~voltage, ncol=1)
plotly::ggplotly()


#' General plan:
#' 1. Find and parse fragments for a given compound
#' 2. Group those fragments into clusters (currently done via custom mz_group)
#' 3. Identify "good" fragments - those that are robust across multiple files
#'   - Done by hierarchical clustering (hclust)
#'   - Cluster fragments in m/z space
#'     - Could be intensity too but not currently implemented
#'   - Cut the tree into clusters based on a cut height that produces the most
#'     groups with 5 +/- 1 fragment
#'  4. Calculate mean m/z and intensity for each cluster of "good" fragments
#'  5. Repeat for each voltage separately
#'  6. Repeat for each compound separately
# Takes about a minute to run
consensus_MS2 <- mock.experimental %>%
  mutate(polarity=str_extract(filename, "pos|neg")) %>%
  distinct(compound_name, polarity) %>%
  # Loop over each compound separately, expect a single character vector back
  mutate(consensus_MS2 = map_chr(compound_name, function(cmpd_name) {
    print(cmpd_name)
    # Grab the fragments associated with a given file
    # Parse into tidy format (cols for name, mz, int)
    # Normalize to % of largest fragment in a given file/voltage
    # Group in m/z space
    # Remove fragments less than 1% the intensity of the largest one
    meth_frags <- init_dat %>%
      filter(compound_name==cmpd_name) %>%
      separate_rows(MS2, sep = "; ") %>% ## standardize and scale.
      separate(MS2, into = c("mz", "int"), sep = ", ") %>%
      mutate(mz=as.numeric(mz)) %>%
      mutate(int=as.numeric(int)) %>%
      group_by(filename, voltage) %>%
      mutate(int=int/max(int)*100) %>%
      ungroup() %>%
      arrange(desc(int)) %>% ## finish standardizing and scaling
      mutate(mz_group=mz_group(mz, ppm = 10)) %>% ## group using the mz grouping function
      filter(int>1)

    ## Within in standardized, grouped compound,
    ## loop over each voltage separately, expect single character vector back
    all_volts <- sort(unique(meth_frags$voltage))
    map_chr(all_volts, function(volts_i){
      good_frags <- meth_frags %>%
        filter(voltage==volts_i)

      # Create hierarchical clustering tree
      # Using Euclidean distance in m/z space
      hc <- good_frags %>%
        select(mz) %>%
        dist() %>%
        hclust()

      # Test multiple possible cut heights
      # Typically around 0.01 seems to work well but does vary (why?)
      cutree_output <- cutree(hc, h = 10^c(-5:5)) %>%
        as.data.frame() %>%
        cbind(good_frags, .) %>%
        pivot_longer(cols = -everything(names(good_frags)), names_to = "cut_height",
                     values_to = "cluster") %>%
        mutate(cut_height=as.numeric(cut_height))
      # Choose the best possible cut height based on the one that produces the most
      # fragment clusters of 5 (because there's 5 files and we want "good"
      # fragments to show up in all 5 files)
      best_cut_height <- cutree_output %>%
        group_by(cut_height, cluster) %>%
        count() %>%
        group_by(cut_height) %>%
        filter(n%in%c(4:6, 9:11)) %>%
        count() %>%
        ungroup() %>%
        arrange(desc(n)) %>%
        filter(n==n[1]) %>%
        summarize(cut_height=10^floor(mean(log10(cut_height)))) %>%
        pull(cut_height)
      # Use ideal cut height to separate tree
      # Only keep clusters that are in 5 +/- 1 file (one dups/missing ok),
      # or 10 +/- 1 for "bimodal superipmosed lollipop effect"
      consistent_frags <- cutree_output %>%
        filter(cut_height==best_cut_height) %>%
        group_by(mz_group) %>%
        filter(n()%in%c(4:6, 9:11))
      if(nrow(consistent_frags)==0){
        return("") # Sometimes you get zero frags that are "good", oh well
      }

      # Calculate average statistics for each fragment
      # Mean for m/z because it feels right?
      # Median for intensity because we want to avoid skew due to outliers
      # Paste into consistent format string ("mz1, int1; mz2, int2; " etc.)
      consistent_frags %>%
        summarise(mz=round(mean(mz), 5),
                  int=round(median(int), 1),
                  voltage=unique(voltage)) %>%
        group_by(voltage) %>%
        summarise(MS2=paste(mz, int, sep = ", ", collapse = "; ")) %>%
        pull(MS2)
    }) %>%
      # Paste into consistent format string
      # "voltage_1 MS2s_from_voltage_1: voltage_2 MS2s_from_voltage_2: " etc.)
      paste(all_volts, ., sep="V ", collapse = ": ")
  }))
