init_msms <- data.frame(MS2="20V 118.0865, 100; 112.1234, 50; 97.4962, 10: 35V 118.0865, 100; 186.8121, 35; 133.7070, 5")
msmsdf <- init_msms %>%
  separate_rows(MS2, sep = ": ") %>%
  separate(MS2, into = c("voltage", "MS2"), sep = "V ") %>%
  separate_rows(MS2, sep = "; ") %>%
  separate(MS2, into = c("mz", "int"), sep = ", ") %>%
  mutate(across(.fns = as.numeric))

msmsdf %>%
  group_by(voltage) %>%
  summarise(MS2=paste(mz, int, sep = ", ", collapse = "; ")) %>%
  summarise(MS2=paste(voltage, MS2, sep="V ", collapse = ": ")) %>%
  as.data.frame()
