# Copyright (c) 2017 Alton Barbehenn and Michael Barbehenn

#!/usr/bin/r

## in the final program Density * FieldSize determines number mice
## and 4 days of trapping determines number of visits
## TODO: look at number caught per grid ring (4 on this 8x8 grid)
## TODO: plot of how well the denisty estimation scales with the number of traps


studySim <- function(ts, fs, np, delta, nv) {
  ## ------ Fixed Stuff ------
  gx <- NULL
  gy <- NULL
  b <- (ts/2)*(2*c(0:7)-7)
  gx <- rep(b, times=length(b))
  gy <- rep(-b, each=length(b))
  
  ## ------ Variable Stuff ------
  mx <- NULL
  my <- NULL
  vx <- NULL
  vy <- NULL
  sx <- matrix(nrow=nv+1,ncol=np)
  sy <- matrix(nrow=nv+1,ncol=np)
  for (ii in c(1:np)) { #we can do this without the for-loop by generating batches of rv's (np, nv, and/or np*nv)
    x <- fs*(0.5-runif(1))
    y <- fs*(0.5-runif(1))
    fx <- x+rnorm(nv) #forage coordinates
    fy <- y+rnorm(nv)
    vx <- append(vx, fx) #all forage locations
    vy <- append(vy, fy)
  }
  
  mice <- data.frame(trap=NULL, day=NULL)
  for (mouse in 1:np) {
    v_range <- ((mouse-1)*nv + 1):(mouse*nv)
    possible_traps <- traps(trap_x = gx, trap_y = gy, forage_x = vx[v_range], forage_y = vy[v_range], delta = delta, nv = nv)
    trap <- possible_traps$trap[!is.na(possible_traps$trap)][1] #first trap the mouse gets caught in (trap = NA indicates that the mouse didn't get caught)
    day <- which(possible_traps$trap==trap)[1] #caught on first day mouse got to that trap
    mice <- rbind(mice, c(trap, day))
  }
  names(mice) <- c("trap", "day")
  mice <- na.omit(mice)
  
  return(mice)
}


oneRun <- function(print=FALSE)
{
  ## ------ Fixed Stuff ------
  gx <- NULL
  gy <- NULL
  b <- (ts/2)*(2*c(0:7)-7)
  # for (ii in b) {
  #   for (jj in b) {
  #     gx <- append(gx,ii)
  #     gy <- append(gy,jj)
  #   } 
  # }
  gx <- rep(b, times=length(b))
  gy <- rep(-b, each=length(b))
  
  # could also use rect or box instead of lines
  ax <- c(-fs/2, fs/2,fs/2,-fs/2,-fs/2)
  ay <- c(-fs/2,-fs/2,fs/2, fs/2,-fs/2)

  ## ------ Variable Stuff ------
  mx <- NULL
  my <- NULL
  vx <- NULL
  vy <- NULL
  sx <- matrix(nrow=nv+1,ncol=np)
  sy <- matrix(nrow=nv+1,ncol=np)
  for (ii in c(1:np)) { #we can do this without the for-loop by generating batches of rv's (np, nv, and/or np*nv)
    x <- fs*(0.5-runif(1))
    y <- fs*(0.5-runif(1))
    mx <- append(mx, x) #mouse location
    my <- append(my, y)
    fx <- x+rnorm(nv) #forage coordinates
    fy <- y+rnorm(nv)
    vx <- append(vx, fx) #all forage locations
    vy <- append(vy, fy)
    sx[,ii] <- c(x,fx) #all forage locations for each mouse
    sy[,ii] <- c(y,fy) #first value is the mouse's location, the rest are it's forage location
  }
  
  x_min = min(ax,gx,mx,vx)
  x_max = max(ax,gx,mx,vx)
  y_min = min(ay,gy,my,vy)
  y_max = max(ay,gy,my,vy)
  
  
  ## ------ Determine Capture ---
  ## This approach is statistically the same as checking whether each forage is caught at the time of the forage
  catch_count <- rep(0, times = length(gx))
  for (mouse in 1:np) {
    v_range <- ((mouse-1)*nv + 1):(mouse*nv)
    possible_traps <- traps(trap_x = gx, trap_y = gy, forage_x = vx[v_range], forage_y = vy[v_range], delta = delta, nv = nv)
    #trap <- possible_traps[!is.na(possible_traps)][1] #first trap the mouse gets caught in (trap = NA indicates that the mouse didn't get caught)
    trap <- possible_traps$trap[!is.na(possible_traps$trap)][1]
    
    if (!is.na(trap)) {
      catch_count[trap] = catch_count[trap] + 1
    }
  }
  catch_count <- matrix(catch_count, nrow = 8, ncol = 8, byrow = TRUE)


  ## ------ Draw Results ------
  
  # apparently c(vx,fx) does the same as append
  if (print) {
    pdf()
  }
  
  ## start new plot base layer (build up from visits)
  plot(c(x_min,x_max),c(y_min,y_max),type="n",xlab="",ylab="")
  title(sprintf("%s=%d %s=%d %s=%3.2f", "np",np,"nv",nv,"ts",ts))
  rect(-fs/2,-fs/2,fs/2,fs/2,border="green")
  
  points(vx,vy,col="red",pch=",")
  
  # create independent clusters
  for (ii in c(1:np)) {
    tx <- NULL
    ty <- NULL
    for (jj in c(2:(nv+1))) {
      tx <- append(tx, c(sx[1,ii], sx[jj,ii]))
      ty <- append(ty, c(sy[1,ii], sy[jj,ii]))
    }
    lines(tx,ty,col="red")
  }
  
  # overlay HRC, Grid
  points(mx,my,col="blue",pch=20)
  points(gx,gy,col="black",lwd=2,pch=0)

  if (print) {
    dev.off()
  }
  
  return(catch_count)
}


#This function is exactly the same as oneRun (above), but it does not attempt to print anything
#or even generate the data to print. This is to speed it up.
matrixSim <- function() {
  ## ------ Fixed Stuff ------
  gx <- NULL
  gy <- NULL
  b <- (ts/2)*(2*c(0:7)-7)
  gx <- rep(b, times=length(b))
  gy <- rep(-b, each=length(b))
  
  ## ------ Variable Stuff ------
  vx <- NULL
  vy <- NULL
  for (ii in c(1:np)) { #we can do this without the for-loop by generating batches of rv's (np, nv, and/or np*nv)
    x <- fs*(0.5-runif(1))
    y <- fs*(0.5-runif(1))
    fx <- x+rnorm(nv) #forage coordinates
    fy <- y+rnorm(nv)
    vx <- append(vx, fx) #all forage locations
    vy <- append(vy, fy)
  }

  ## ------ Determine Capture ---
  ## This approach is statistically the same as checking whether each forage is caught at the time of the forage
  catch_count <- rep(0, times = length(gx))
  for (mouse in 1:np) {
    v_range <- ((mouse-1)*nv + 1):(mouse*nv)
    possible_traps <- traps(trap_x = gx, trap_y = gy, forage_x = vx[v_range], forage_y = vy[v_range], delta = delta, nv = nv)
    trap <- possible_traps$trap[!is.na(possible_traps$trap)][1] #first trap the mouse gets caught in (trap = NA indicates that the mouse didn't get caught)
    
    if (!is.na(trap)) {
      catch_count[trap] = catch_count[trap] + 1
    }
  }
  catch_count <- matrix(catch_count, nrow = 8, ncol = 8, byrow = TRUE)
  
  return(catch_count)
}


#Returns the traps, if any, that the mouse gets caught in for forage
#Pretty sure this approach maintains the foraging order
traps <- function(trap_x, trap_y, forage_x, forage_y, delta, nv){
  dx <- lapply(forage_x, function(x) which(abs(x - trap_x) <= delta))
  dx <- lapply(dx, function(x) unique(trap_x[x]))
  dx <- lapply(dx, function(x) ifelse(length(x)>0, x, NA))
  dx <- unlist(dx)

  dy <- lapply(forage_y, function(x) which(abs(x - trap_y) <= delta))
  dy <- lapply(dy, function(x) unique(trap_y[x]))
  dy <- lapply(dy, function(x) ifelse(length(x)>0, x, NA))
  dy <- unlist(dy)
  
  possible_traps <- lapply(1:length(dx), function(x) if(is.na(dx[x]) || is.na(dy[x])) {NA} else {intersect(which(trap_x==dx[x]), which(trap_y==dy[x]))})
  possible_traps <- unlist(possible_traps)
  possible_traps <- data.frame(trap=possible_traps, day=1:nv)
  
  return(possible_traps)
}

