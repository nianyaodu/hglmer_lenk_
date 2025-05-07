rng_orderdiri <- function(n,shape,niter=100){
  K <- length(shape)
  result <- matrix(NA_real_,nrow=n,ncol=K)
  for(l in 1:n){
    # Initialize x
    x <- rgamma(K,shape)
    # run gibbs
    for(i in 1:niter){
      v <- x^(shape-1)*runif(K)
      lb <- v[1]^(1/(shape[1]-1)); u <- runif(1)
      x[1] <- lb - log(1-u)
      
      for(k in 2:K){
        lb <- max(v[k]^(1/(shape[k]-1)),x[k-1]); u <- runif(1)
        x[k] <- lb - log(1-u)
      }
    }
    x <- x/sum(x)
    result[l,] <- x
  }
  return(result)
}
