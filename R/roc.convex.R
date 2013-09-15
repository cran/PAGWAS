roc.convex <- function(sens,spec) {
  n <- length(spec)
  ind <- i <- 1
  area <- 0
  while (i < n) {
    j <- (i+1):n
    if (max(sens[j] - sens[i]) > 1e-12) {
      slope <- (spec[j] - spec[i])/(sens[j] - sens[i])
      new.i <- i + which.max(slope)
      area <- area + 0.5*(spec[new.i]+spec[i])*(sens[new.i]-sens[i])
    } else {
      new.i <- n
    }
    i <- new.i
    ind <- c(ind,i)
  }
area	
	}


