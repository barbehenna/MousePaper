# Test script for the full simulation backend

# Test mouse foraging
nsamples = 1000
plot(GenMouse(nforages = nsamples, fieldSize = 13.5)) # forages as a bi-variate normal
nortest::lillie.test(rnorm(nsamples)) # calabrate
nortest::lillie.test(runif(nsamples))
nortest::lillie.test(GenMouse(nforages = nsamples, fieldSize = 13.5)[,1]) # test normality of x coordinates
nortest::lillie.test(GenMouse(nforages = nsamples, fieldSize = 13.5)[,2]) # test normality of y coordinates




# Test mice home-ranges
nsamples = 1000
plot(t(sapply(seq(nsamples), function(x) GenMouse(nforages = 1, fieldSize = 13.5))))
abline(v = c(-6.75, 6.75), col = 'red') # add field boarder
abline(h = c(-6.75, 6.75), col = 'red')
# the reason some points are outside is because this is the first forage location, not the home,
# so they can be +- a few sigma from their home, most are inside square. If you modify the function
# to return just NumericVector of homes, it works perfectly.




# Test is mouse caught
mouse <- GenMouse(4, fieldSize = 13.5)
Traps <- GenTraps(nrings = 8, trapspacing = 0.5)

# day 1 (increase mouse index 1 for other days)
dxdy <- t(apply(Traps, MARGIN = 1, FUN = function(x) {return(abs(x - mouse[1,]))}))
which(apply(dxdy, MARGIN = 1, FUN = function(x) {max(x) < 0.3})) - 1 # which trap(s) catch a mouse
Traps[apply(dxdy, MARGIN = 1, FUN = function(x) {max(x) < 0.3}), ]
isCaught(forages = mouse, Traps = Traps, catchRadius = 0.3)
# iterated enough times to check outputs that I'm happy the answers match for isCaught()