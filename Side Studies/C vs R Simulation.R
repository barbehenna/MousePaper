# Copyright (c) 2017 Alton Barbehenn

# This script is to compare the two versions of simulating the trapping data. It looks
# like the C function is significantly faster and that the resulting distributions are 
# identical, so I didn't mess up the transcription from R to C.

library(raster)

ts = 1.5
delta = 0.5
fs = 20
np = 1000
nv = 4

# Compare the performance of the two methods
system.time(old <- studySim(ts, fs, np, delta, nv)) #~8min for np=100000
system.time({new <- trapSim1(ts, fs, np, delta, nv)
             new <- data.frame(trap=new[,1], day=new[,2])
             new <- na.omit(new)})


# Check more percisely that they give the same result by setting a common
# seed value and looking at the results.
set.seed(42)
old <- studySim(ts, fs, np, delta, nv)

set.seed(42)
new <- trapSim1(ts, fs, np, delta, nv)
new <- data.frame(trap=new[,1], day=new[,2])
new <- na.omit(new)
new$trap <- new$trap + 1 # Adjust for difference between c++ and R indexing
new$day <- new$day + 1

# They look the same!
View(old)
View(new)

nrow(old) == nrow(new) # Check for the same number of rows
diff <- NULL
for (i in 1:nrow(old)) {
  diff <- rbind(diff, c(abs(new$trap[i]-old$trap[i]), abs(new$day[i] - old$day[i])))
}
sum(diff[,1]) == 0 # Check for different values between the data frames
sum(diff[,2]) == 0


# Other, less rigorous tests to make sure the methods are the same
# Check that the days get roughtly the correct distributes
table(old$day)
table(new$day)
# Look good to me

# Look at data by trap visually
old_traps <- as.data.frame(table(old$trap))
old_mat <- matrix(old_traps$Freq, nrow = 8, ncol = 8)
old_rast <- raster(xmn = 0, xmx = 8, ymn = 0, ymx = 8, nrows = 8, ncols = 8)
old_rast[] <- old_mat
plot(old_rast, main = "R method")

new_traps <- as.data.frame(table(new$trap))
new_mat <- matrix(new_traps$Freq, nrow = 8, ncol = 8)
new_rast <- raster(xmn = 0, xmx = 8, ymn = 0, ymx = 8, nrows = 8, ncols = 8)
new_rast[] <- new_mat
plot(new_rast, main = "C method")


# Check that the distiubutions are the same by ring
rings <- c(4,4,4,4,4,4,4,4, 
           4,3,3,3,3,3,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,3,3,3,3,3,4, 
           4,4,4,4,4,4,4,4)
old_traps$ring <- rings
new_traps$ring <- rings

t.test(old_traps$Freq[old_traps$ring == 4], new_traps$Freq[new_traps$ring == 4])
t.test(old_traps$Freq[old_traps$ring == 3], new_traps$Freq[new_traps$ring == 3])
t.test(old_traps$Freq[old_traps$ring < 3], new_traps$Freq[new_traps$ring < 3])


