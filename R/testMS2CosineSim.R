library(plotly)
library(tidyverse)
source("R/utils.R")

AssignCluster <- function(mz, int, theoretical) {
  ## Function that takes in single experimental fragment
  ## and assigns to cluster of theoretical data
  cor_mz <- mz/7.331762e-05
  closestcluster <- theoretical %>%
    group_by(cluster) %>%
    mutate(mz = mz/7.331762e-05) %>%
    summarise(mean_mz = mean(mz), mean_intensity = mean(intensity)) %>%
    mutate(mz_dist=cor_mz-mean_mz, int_dist=int-mean_intensity) %>%
    mutate(total_dist=sqrt(mz_dist^2+int_dist^2)) %>%
    filter(total_dist < 10) %>% # Will need to be fiddled with
    arrange(total_dist) %>%
    slice(1) %>%
    pull(cluster)
  if(length(closestcluster)==0){
    return(0)
  }
  closestcluster
}


## Load and tidy spectra. Combine theoretical to one dataframe.
experimental <- c("53.1943054199219, 7666.30419921875; 56.5771636962891, 7860.208984375; 60.0564270019531, 212678.671875; 61.9458084106445, 6289.337890625; 70.0658721923828, 159965.25; 72.0815277099609, 9107.90625; 77.6530532836914, 6911.32958984375; 83.9745559692383, 7009.45849609375; 94.0988616943359, 8369.44921875; 100.966766357422, 8566.7353515625; 106.529121398926, 8185.25439453125; 110.386734008789, 6960.23876953125; 112.087326049805, 22621.90234375; 114.102821350098, 7501.333984375; 116.070930480957, 247170.75; 119.085800170898, 10558.1015625; 130.097702026367, 93940; 157.108703613281, 23185.689453125; 158.092559814453, 108007.5078125; 160.639587402344, 7734.41455078125; 175.086959838867, 21641.658203125; 175.112426757812, 37277.390625; 175.119094848633, 815061.9375; 175.148025512695, 15538.94140625") %>%
  MakeScantable()

theoretical1 <- c("58.48081, 4161; 59.18518, 3224; 60.05643, 239657; 60.44617, 3626; 63.6882, 3972; 67.56001, 3296; 70.06588, 153105; 72.08139, 14450; 77.71708, 3770; 80.89391, 3722; 91.95894, 3489; 111.85516, 4268; 112.08733, 26194; 114.1029, 6346; 115.08686, 4902; 116.07099, 261445; 119.88438, 3608; 120.86121, 3727; 130.09776, 102881; 132.01106, 4013; 137.13593, 3684; 141.06642, 3932; 157.10863, 29600; 158.09254, 143890; 175.11913, 1141570; 191.26741, 3905") %>%
  MakeScantable()
theoretical2 <- c("60.05638, 284683; 69.57538, 3850; 70.06581, 176540; 70.52679, 3517; 72.08147, 15337; 76.89897, 3626; 90.39793, 4155; 107.40523, 3599; 112.08733, 27236; 114.10276, 12173; 115.08661, 7533; 116.07089, 318679; 125.31926, 3882; 125.42632, 3896; 130.09766, 116939; 133.09732, 4621; 133.31738, 3728; 140.08226, 6550; 157.10857, 40142; 158.09248, 156906; 175.11899, 1233636; 183.13173, 3683; 185.21806, 4119") %>%
  MakeScantable()
theoretical3 <- c("60.05642, 309212; 68.92125, 3973; 69.25908, 4062; 70.06586, 189675; 72.08142, 18393; 86.44745, 3826; 112.08734, 35409; 113.0714, 4917; 114.10283, 11231; 115.08681, 8407; 116.07095, 346820; 124.10611, 3765; 130.09773, 133121; 144.36554, 3986; 146.70622, 3969; 146.77718, 3276; 146.91052, 3815; 157.1086, 38035; 158.09256, 166139; 172.14598, 4010; 175.11908, 1313840; 187.29317, 3736") %>%
  MakeScantable()
theoretical4 <- c("60.05646, 302463; 64.17995, 4225; 68.11596, 3495; 70.06591, 182791; 72.08153, 15981; 97.02829, 4277; 106.29856, 4475; 112.08753, 28404; 114.10301, 10143; 115.08712, 7610; 116.07103, 362345; 124.09547, 4112; 130.09779, 140969; 149.39337, 4322; 152.52785, 3761; 157.10854, 38673; 158.09264, 162281; 164.09802, 3710; 164.94202, 3501; 169.94133, 3633; 175.11919, 1338598; 186.36243, 3805")%>%
  MakeScantable()
theoretical5 <- c("50.68667, 3597; 60.05649, 325893; 63.70285, 3851; 70.06595, 214896; 72.08157, 23573; 83.3022, 3745; 94.37806, 3664; 94.64034, 3611; 97.88653, 3667; 98.02602, 3807; 104.53662, 4044; 107.66678, 3768; 112.08753, 29813; 113.0714, 4161; 114.10293, 14074; 115.08696, 10128; 116.07111, 377898; 126.07323, 3911; 130.09787, 130216; 131.95659, 4670; 135.25867, 3535; 135.69394, 3750; 135.97992, 3850; 157.10873, 46256; 158.09279, 177553; 160.53001, 4332; 175.11931, 1471160") %>%
  MakeScantable()

complete.theoretical <- bind_rows(list(a = theoretical1, b = theoretical2, c = theoretical3, d = theoretical4, e = theoretical5), .id = "id")

## Create distance matrix and perform hierarchical cluster analysis
hc <- complete.theoretical[ ,c("mz", "intensity")] %>%
  mutate(intensity = intensity * 7.331762e-05) %>% # TODO This value will need to get fiddled with
  dist() %>%
  hclust()
hc %>%
  as.dendrogram() %>%
  plot(ylim =c(0, 0.01))
complete.theoretical$cluster <- hc %>%
  cutree(h = 0.001) # TODO This value will need to get fiddled with

## Plot clustered theoretical data including the mean theoretical values,
## and add experimental points.
complete.theoretical %>%
  group_by(cluster) %>%
  mutate(mz = mz/7.331762e-05) %>% ## TODO Ditto value check
  mutate(mean_mz = mean(mz),
         mean_intensity = mean(intensity)) %>%
  ggplot(aes(x = mz, y = intensity, color = factor(cluster))) +
  geom_jitter(width = 0, height = 0) +
  geom_point(data = experimental %>% mutate(mz = mz/7.331762e-05), color = "black") +
  geom_point(aes(x = mean_mz, y = mean_intensity), color = "yellow", shape = 2, size = 3)
ggplotly()


complete.theoretical %>% # TODO What is the point of this step?
  filter(mz > 125, mz < 150) %>%
  summarise(sd_mz = sd(mz),
            sd_int = sd(intensity))

## Assign cluster to experimental data and plot
exp_assigned <- experimental %>%
  rowwise() %>%
  mutate(cluster=AssignCluster(mz, intensity, complete.theoretical))
complete.theoretical %>%
  group_by(cluster) %>%
  mutate(mz = mz/7.331762e-05) %>%
  mutate(mean_mz = mean(mz),
         mean_intensity = mean(intensity)) %>%
  ggplot(aes(x = mz, y = intensity, color = factor(cluster))) +
  geom_jitter(width = 0, height = 0) +
  geom_point(data = exp_assigned %>% mutate(mz = mz/7.331762e-05), size=3) +
  geom_point(aes(x = mean_mz, y = mean_intensity), color = "yellow", shape = 2, size = 3)


## original
weight1 <- (scan1[, "mz"] ^ 2) * sqrt(scan1[, "intensity"]) # high mz, high intensity = larger weight according to current function
weight2 <- (scan2[, "mz"] ^ 2) * sqrt(scan2[, "intensity"])

diff.matrix <- sapply(scan1[, "mz"], function(x) scan2[, "mz"] - x) # calculate difference between all mzs
## use dist here instead??

same.index <- which(abs(diff.matrix) < 0.02, arr.ind = TRUE) # filter for mzs within given range (currently 0.02, which is including a lot)

cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
  (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))
