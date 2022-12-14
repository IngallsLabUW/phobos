
# Profvis checks which line of code is running every 20 microseconds
# and returns it in a nice, interactive profile that can be opened
# in RStudio. The bottlenecks in your code and/or memory are highlighted
# as black bars on the right side - the higher the bar, the larger the
# relative time/usage required.

# When the profile tab opens automatically, make sure to check out both the
# flame graph tabs AND the data tabs.

library(profvis)

profvis({
  newfile <- tempfile(fileext = ".csv")
  write.csv(x = `names<-`(data.frame(rep(LETTERS, each=1000), 1:26000, rnorm(26000)),
                          c("Letters", "Numbers", "Random")), file = newfile)
  
  data <- read.csv(newfile)
  
  summary(data)
  
  str(data)
  
  data$Random2 <- rnorm(1:26000)
  
  cor(x = data$Random, y = data$Random2)
  
  lm(data$Random2~data$Random)
  
  plot(data$Random, data$Random2)
  
  file.remove(newfile)
})



# Wraps nicely into functions
# But the profile will need to be expanded by clicking, otherwise it's
# simply the time taken by the function as a whole
timeTaker <- function(iterations){
  newfile <- tempfile(fileext = ".csv")
  write.csv(x = `names<-`(data.frame(rep(LETTERS, each=iterations), 
                                     1:(26*iterations), rnorm(26*iterations)),
                          c("Letters", "Numbers", "Random")), file = newfile)
  
  data <- read.csv(newfile)
  
  summary(data)
  
  str(data)
  
  data$Random2 <- rnorm(1:(26*iterations))
  
  cor(x = data$Random, y = data$Random2)
  
  lm(data$Random2~data$Random)
  
  plot(data$Random, data$Random2)
  
  file.remove(newfile)
}

profvis({
  for(i in 1:10){
    timeTaker(100)
  }
})

profvis({
  sapply(1:10, timeTaker)
})



# Check overhead due to sampling
init_time <- Sys.time()
timeTaker(1000)
Sys.time() - init_time
#About 1.5 seconds on my laptop

init_time <- Sys.time()
profvis({timeTaker(1000)})
Sys.time() - init_time
#About 5 seconds on my laptop

# Drop the interval to increase resolution
# Throws a warning though
init_time <- Sys.time()
profvis({timeTaker(10000)}, interval = 0.001)
Sys.time() - init_time
#About 6 seconds


# Increase sample size to check linearity
init_time <- Sys.time()
timeTaker(10000)
Sys.time() - init_time
#About 10 seconds on my laptop

init_time <- Sys.time()
profvis({timeTaker(10000)})
Sys.time() - init_time
#About 50 seconds on my laptop

# Drop the interval to increase resolution
# Throws a warning though
init_time <- Sys.time()
profvis({timeTaker(10000)}, interval = 0.001)
Sys.time() - init_time
#About 50 seconds





# Alternatively, you can use RStudio's built-in version
# Use the "Profile -> Start profiling" tab at the top of RStudio
# Then copy and paste the below code into the terminal and run
# BUT you won't get the interactive flame graph, though the data is preserved

newfile <- tempfile(fileext = ".csv")
write.csv(x = `names<-`(data.frame(rep(LETTERS, each=1000), 1:26000, rnorm(26000)),
                        c("Letters", "Numbers", "Random")), file = newfile)

data <- read.csv(newfile)

summary(data)

str(data)

data$Random2 <- rnorm(1:26000)

cor(x = data$Random, y = data$Random2)

lm(data$Random2~data$Random)

plot(data$Random, data$Random2)

file.remove(newfile)

# Then click the "Stop profiling" button at the top of the terminal or in the tab