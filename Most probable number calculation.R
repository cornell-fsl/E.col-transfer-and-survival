The R script used to implement the method outlined by Cochran to calculate the most probable number of E. coli cells per sample.

# Author: David Kent (dk657@cornell.edu) and Daniel Weller (dlw263@cornell.edu)
# Date: December 23, 2016
#
# The following is an implementation of the method described in:
# Cochran, William G. "Estimation of bacterial densities by means of the
# 'most probable number'." Biometrics 6.2 (1950): 105-116.

# MPN Function to determine the most probable number of bacteria in a given
# sample given a specific pattern of positives tubes and replicates

mpn <- function(positives, tubes = 6, replicates = 2, greatestvolume = 10, dilution = 100) {
  
  combination <- function(n,k) factorial(n)/(factorial(k)*factorial(n-k)) #For convenience
  
  # Tubes is the number of dilutions used
  # Replicates is the number of replicates per dilution
  # Greatestvolume is the volume of sample per tube in the least diluted level
  # Dilution is the dilution factor (e.g. 10, 2, 5, etc.)
  # Positives is a (k x tubes) matrix containing the number of positive tubes at each dilution level
  
  n <- replicates
  
  
  v <- greatestvolume
  if (tubes > 1) {
    for(k in 1:(tubes-1)){
      v <- c(v,greatestvolume*dilution^(-k))
    }
  }
  
  resultVector <- vector("numeric",dim(positives)[1])
  
  for(k in 1:length(resultVector)) {
    cat(k," of ",length(resultVector),"\n")
    
    # Number of sterile tubes
    s <- n - positives[k,]
    
    # Probability that s out of n samples of size v are sterile for a given density
    prob <- function(density) {
      result <- 1
      
      for(i in 1:length(s)) {
        result <- result*combination(n,s[i])*exp(-v[i]*s[i]*density)*(1-exp(-v[i]*density))^(n-s[i])
      }
      
      return(as.numeric(result))
      
    }
    
    lowerBound <- 0 # Start with a small interval (0 to 1) and look for the maximum
    upperBound <- 1 # probability until it differs from the upper bound of the interval.
    MPN <- 1/5      # We are assuming that the function increases from 0 and has only one maximum.
    while(signif(MPN*5,digits=1)==signif(upperBound,digits=1)) {
      MPN <- optimize(f=prob,interval=c(lowerBound,upperBound),maximum=T)$maximum
      lowerBound <- upperBound
      upperBound <- upperBound*5
    }
    
    resultVector[k] <- signif(MPN,digits=2)
    
    
  }
  
  return(cbind(positives,MPN=resultVector))
}

