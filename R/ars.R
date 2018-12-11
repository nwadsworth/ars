# get h(x)=log(g(x)) ----------------------------------------------------------
h <- function(f, x = 0) {
  t <- log(f(x))
  return(t)
}

# get h'(x) -------------------------------------------------------------------
#' @importFrom pracma fderiv
d <- function(f, x = 0) {
  return(fderiv(f, x, n = 1, method = "central")/f(x))
}

# initialize the data matrix from user input ----------------------------------
#' @importFrom utils head
#' @importFrom utils tail
initial <- function(formulas, min, max, xinit) {
  xinit=sort(unique(xinit))
  invalue = 1:(4 * length(xinit))
  data <- matrix(invalue, nrow = length(xinit), ncol = 4)
  data[, 1] <- xinit[order(xinit)]
  data[, 2] <- h(formulas, data[, 1])
  data[, 3] <- d(formulas, data[, 1])
  zt <- (tail(data[, 2], -1) - head(data[, 2], -1) - tail(data[, 1], -1) * tail(data[, 3], -1) + head(data[, 1], -1) * head(data[, 3], -1)) /
    (head(data[, 3], -1) - tail(data[, 3], -1))
  data[, 4] <- c(zt, max)
  return(data)
}

# get u_k(x) ------------------------------------------------------------------
u <- function(data, x) {
  i <- findInterval(x, data[, 4]) + 1
  return(data[i, 2] + (x - data[i, 1]) * data[i, 3])
}

# get l_k(x) ------------------------------------------------------------------
l <- function(data, x) {
  i <- findInterval(x, data[, 1])
  if (i == nrow(data) || i == 0)
    return(-Inf)
  return(((data[i + 1, 1] - x) * data[i, 2] + (x - data[i, 1]) * data[i + 1, 2])/(data[i + 1, 1] - data[i, 1]))
}

# get x* ----------------------------------------------------------------------
#' @importFrom stats runif
#' @importFrom utils head
#' @importFrom utils tail
get_sample <- function(f, y, hp, zt, min) {
  z = c(min, zt)
  # make sure hp has no 0, otherwise will be mistakes
  intgration <- f(y)/hp * exp(-hp * y) * (exp(hp * tail(z, -1)) - exp(hp * head(z, -1)))
  intgration <- c(0, intgration)
  tot = sum(intgration)
  cumint <- cumsum(intgration)/tot
  t <- runif(1)
  index <- findInterval(t, cumint)

  s <- log((t - cumint[index]) * tot * exp(hp[index] * y[index]) * hp[index]/f(y[index]) + exp(hp[index] * z[index]))/hp[index]
  return(s)
}

# get u_k(x) ------------------------------------------------------------------
#' @importFrom stats runif
test <- function(data, f, x) {
  # the first return value is accept, the second is update or not.
  w = runif(1)
  if (w <= exp(l(data, x) - u(data, x))) { # squeezing test
    return(c(TRUE, FALSE))
  } else if (w <= exp(h(f, x) - u(data, x))) { # rejection test
    return(c(TRUE, TRUE))
  } else {
    return(c(FALSE, TRUE))
  }
}

# add x*, h(x*), h'(x*) to data matrix ----------------------------------------
update_step <- function(data, f, x_star, min, max) {
  newrow <- c(x_star, h(f, x_star), d(f, x_star), 0)  # the new row will be added
  pos <- findInterval(x_star, data[, 1])
  if (pos == length(data[, 1])) {
    # x_star is the largest number
    newrow[4] <- max
    data[pos, 4] <- (h(f, x_star) - h(f, data[pos, 1]) - x_star * d(f, x_star) + data[pos, 1] * d(f, data[pos, 1])) /
      (d(f, data[pos, 1]) - d(f, x_star))
  } else if (pos == 0) {
    # x_star is the smallest number
    newrow[4] <- (h(f, data[pos + 1, 1]) - h(f, x_star) - data[pos + 1, 1] * d(f, data[pos + 1, 1]) + x_star * d(f, x_star)) /
      (d(f, x_star) - d(f, data[pos + 1, 1]))
  } else {
    # x_star is not the smallest or the largest number
    newrow[4] <- (h(f, data[pos + 1, 1]) - h(f, x_star) - data[pos + 1, 1] * d(f, data[pos + 1, 1]) + x_star * d(f, x_star)) /
      (d(f, x_star) - d(f, data[pos + 1, 1]))
    data[pos, 4] <- (h(f, x_star) - h(f, data[pos, 1]) - x_star * d(f, x_star) + data[pos, 1] * d(f, data[pos, 1])) /
      (d(f, data[pos, 1]) - d(f, x_star))
  }
  data <- rbind(data, newrow[1:4])
  # order the data
  data <- data[order(data[, 1]), ]
  return(data)
}

# Check if the distribution is log-concave or not ----------------------------------------
#' @importFrom stats nlm
check_dist<-function(formulas=dnorm, nsamples=100, min=-Inf, max=Inf, xinit=c(-1.5,-0.2, 0.5, 1)) {
  formulas <- Vectorize(formulas)
  h<-function(x) return(log(formulas(x)))
  xinit=sort(unique(xinit))
  xinit=xinit[(xinit>min)&(xinit<max)]
  #if ture, Non-log-Concative
  dd <-function(x) return(-fderiv(h,x, n=2, method = "central"))
  opt<-nlm(dd,tail(xinit,1))
  # Secondary derivative large than 0 (or 1e-4) and within the boundary, treat it as non-log-cancave
  if(((-opt$minimum)>1e-4)&(opt$estimate<max)&(opt$estimate>min)){
    print("This function is not a log-concave function, this function doesn't work!")
    return(FALSE)
  }

  if(!sum(!(abs(fderiv(h, xinit, n = 2, method = "central"))<=1e-4))){

    if(!sum(!(abs(fderiv(h, xinit, n = 1, method = "central"))<1e-4))){
      #uniform distribution
      print("This is a probaly uniform distribution!")
      if(formulas(xinit)<=0){
        print("Density function must be positive!")
        return(FALSE)
      }
      if((max==Inf) | (min==-Inf)){
        print("Please reset the min or max, they must be limited values")
        return(FALSE)
      }
      return(runif(nsamples, min ,max))

    }else{
      # exponential distribution
      lamuda=fderiv(h, tail(xinit,1), n = 1, method = "central")
      if((min==-Inf)|(lamuda>0&max==Inf)){
        print("min should be finite, and max should be finite for a increasing exponential function!")
        return(FALSE)
      }
      w=runif(nsamples)
      res=log(exp(lamuda*min)+w*(exp(lamuda*max)-exp(lamuda*min)))/lamuda
      print("This is probaly an exponential distribution!")
      return(res)
    }

  }

  if((!is.vector(xinit))|(length(xinit)<2)){
    print("The starting points should be a vector with at least two elements WITHIN the Boundary!")
    return(FALSE)
  }
  # if unbounded, the first derivative should be postive and the last should be negative.
  if(((min==-Inf)&(fderiv(h, xinit[1], n = 1, method = "central")<=0))|((max==Inf)&(fderiv(h, tail(xinit,1), n = 1, method = "central")>=0))){
    print("Please reset the starting pointing to make sure some are uphill and some are downhill! ")
    return(FALSE)
  }
  # next step, continued as below
  return(TRUE)
}

# perform the adaptive rejection sampling ----------------------------------------
#' Adaptive Rejection Sampling
#'
#' Perform the adaptive rejection sampling from the specified log-concave density
#'
#' @param f A function of the continuous density we want to sample from. Defaults to standard normal.
#' @param nsamples The number of samples desired. Defaults to 100
#' @param min The minimum of the domain of \code{f}. Defaults to -Inf.
#' @param max The maximum of the domain of \code{f}. Defaults to Inf.
#' @param xinit Vector of starting points for which \code{f} is defined (at least 2).
#' @return A vector (of length \code{nsamples}) of sampled points from the specified distribution.
#' @examples
#' a <- ars(dnorm, nsamples = 100, min = -Inf, max = Inf, xinit = c(-1.5, -0.2, 0.5, 1))
#' g <- function(x) {dchisq(x, df = 3)}
#' @importFrom pracma fderiv
#' @importFrom stats runif
#' @importFrom stats dnorm
#' @importFrom stats dchisq
#' @importFrom assertthat assert_that
#' @export
ars <- function(f = dnorm, nsamples = 100, min = -Inf, max = Inf, xinit = c(-1.5, -0.2, 0.5, 1)) {
  #check inputs nsamples, min, max, xinit
  assert_that((nsamples>=0) && (min<max) && (is.vector(xinit))
              &&(length(unique(xinit))>=2)  &&(is.function(f)))

  xinit=sort(unique(xinit))
  xinit=xinit[(xinit>min)&(xinit<max)]
  res=check_dist(formulas=f, nsamples=nsamples, min=min, max=max, xinit=xinit)
  if(is.numeric(res)){
    return(res)
  }
  if(res==FALSE){
    print("NO sampling result returned!")
    return(NULL)
  }
  #need check xinit for unbounded cases
  data <- initial(formulas = f, min = min, max = max, xinit = xinit)
  nsam = 0
  result = rep(0, nsamples)
  while (nsam < nsamples) {
    sam <- get_sample(f = f, y = data[, 1], hp = data[, 3], z = data[, 4], min = min)
    accept <- test(data, f = f, sam)
    if (accept[1]) {
      nsam = nsam + 1
      result[nsam] = sam
    }
    if (accept[2]) {
      data = update_step(data, f = f, x_star = sam, min = min, max = max)
    }
  }
  info <- paste0("In total, ", length(data[, 1]), " sampling points needed to be stored in the Matrix!")
  print(info)
  return(result)
}
