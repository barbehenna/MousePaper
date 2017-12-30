# Copyright (c) 2017 Alton Barbehenn

# This script is to compare the two versions of simulating the trapping data. It looks
# like the C function is significantly faster and that the resulting distributions are 
# identical, so I didn't mess up the transcription from R to C.

library(raster)

ts = 1.5
delta = 0.5
fs = 20
np = 100000
nv = 4

rings <- c(4,4,4,4,4,4,4,4, 
           4,3,3,3,3,3,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,3,3,3,3,3,4, 
           4,4,4,4,4,4,4,4)

system.time(old <- studySim(ts, fs, np, delta, nv)) #~8min for np=100000
system.time(new <- trapSim1(ts, fs, np, delta, nv))
new <- data.frame(trap=new[,1], day=new[,2])
new <- na.omit(new)

# check that the days get roughtly the correct distributes
table(old$day)
table(new$day)
# look good to me

# look at data by trap visually
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


# check that the distiubutions are the same by ring
old_traps$ring <- rings
new_traps$ring <- rings

t.test(old_traps$Freq[old_traps$ring == 4], new_traps$Freq[new_traps$ring == 4])
t.test(old_traps$Freq[old_traps$ring == 3], new_traps$Freq[new_traps$ring == 3])
t.test(old_traps$Freq[old_traps$ring < 3], new_traps$Freq[new_traps$ring < 3])



