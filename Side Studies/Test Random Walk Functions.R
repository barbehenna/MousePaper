library(Rcpp)
Rcpp::sourceCpp('Documents/Projects/MousePaper/RandomWalk.cpp')

st <- c(0,0)
now <- c(-1,0)

tmp <- t(replicate(10000, step(now, st, 0.25)))
tmp <- as.data.frame(tmp)
names(tmp) <- c("x","y")

ggplot(tmp, aes(x,y)) + geom_density2d()


################################
pose <- c(0.0,0.0)
coords <- matrix(pose, ncol = 2)
for (i in seq(10000)) {
  pose <- step(pose, c(0.0,0.0), 0.25)
  coords <- rbind(coords, pose)
}

coords <- as.data.frame(coords)
names(coords) <- c("x", "y")
coords$num <- seq(nrow(coords))

ggplot(coords, aes(x, y, colour = num)) + geom_path()


################################
ts <- 2
cr <- 0.5
b <- (ts/2)*(2*c(0:7)-7)
gx <- rep(b, times=length(b))
gy <- rep(-b, each=length(b))
trap_coords <- matrix(c(gx,gy), ncol = 2)


nearTrap_Square(c(-6.75, 6.75), trap_coords, cr)
nearTrap_Square(c(2.55, 6.93), trap_coords, cr)
nearTrap_Square(c(9,9), trap_coords, cr)


################################
plot(density(replicate(10000, von_mises(2,2))))


################################
sign(-2)
sign(2)
sign(0)


