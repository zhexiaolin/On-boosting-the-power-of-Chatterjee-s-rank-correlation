# This function computes the revised Chatterjee’s rank correlation $\xi_{n,M}$.
# Input: xvec, yvec are two vectors, M is the number of right nearest neighbors.
# Output: the value of $\xi_{n,M}$.

XIMcalculate<-function(xvec, yvec, M){
  n <- length(xvec)
  xrank <- rank(xvec, ties.method = "random")
  yrank <- rank(yvec, ties.method = "random")
  ord <- order(xrank)
  yrank <- yrank[ord]
  coef.term <- function(m){
    return(sum(pmin(yrank[1:(n-m)], yrank[(m+1):n])) + sum(yrank[(n-m+1):n]))
  }
  coef.sum <- sapply(1:M, coef.term)
  coef.sum <- sum(coef.sum)
  coef.value <- -2+6*coef.sum/((n+1)*(n*M+M*(M+1)/4))
  return(coef.value)
}


# This function computes the test statistics $\xi_{n,M}^{\pm}$.
# Input: xvec, yvec are two vectors, M is the number of right nearest neighbors.
# Output: the value of $\xi_{n,M}^{\pm}$.
XIMstat<-function(xvec, yvec, M){
  n <- length(xvec)
  xrank <- rank(xvec, ties.method = "random")
  yrank <- rank(yvec, ties.method = "random")
  ord <- order(xrank)
  yrank <- yrank[ord]
  yrank.m <- n+1-yrank 
  coef.sum <- 0
  coef.sum.m <- 0
  coef.term <- function(m){
    return(sum(pmin(yrank[1:(n-m)], yrank[(m+1):n])) + sum(yrank[(n-m+1):n]))
  }
  coef.term.m <- function(m){
    return(sum(pmin(yrank.m[1:(n-m)], yrank.m[(m+1):n])) + sum(yrank.m[(n-m+1):n]))
  }
  coef.sum <- sapply(1:M, coef.term)
  coef.sum.m <- sapply(1:M, coef.term.m)
  coef.sum <- sum(coef.sum)
  coef.sum.m <- sum(coef.sum.m)
  coef.stat <- -2+6*max(coef.sum,coef.sum.m)/((n+1)*(n*M+M*(M+1)/4))
  return(coef.stat)
}


# This function computes the simulation statistics $\xi_{n,M}^{\pm(b)}$.
# Input: n is the number of sample, M is the number of right nearest neighbors, B is the number of simulation.
# Output: the vector of $\xi_{n,M}^{\pm(b)}$.
XIMsim<-function(n, M, B){
  XIMsim_single = function(){
    yrank <- sample(1:n,n)
    yrank.m <- n+1-yrank 
    coef.term <- function(m){
      return(sum(pmin(yrank[1:(n-m)], yrank[(m+1):n])) + sum(yrank[(n-m+1):n]))
    }
    coef.term.m <- function(m){
      return(sum(pmin(yrank.m[1:(n-m)], yrank.m[(m+1):n])) + sum(yrank.m[(n-m+1):n]))
    }
    coef.sum <- sapply(1:M, coef.term)
    coef.sum.m <- sapply(1:M, coef.term.m)
    coef.sum <- sum(coef.sum)
    coef.sum.m <- sum(coef.sum.m)
    coef.stat <- -2+6*max(coef.sum,coef.sum.m)/((n+1)*(n*M+M*(M+1)/4))
  }
  coef.sim<-replicate(B, XIMsim_single())
  return(coef.sim)
}


# This function performs the simulation based test using test statistics and simulation statistics.
# Input: XIMstat is the value of $\xi_{n,M}^{\pm}$, XIMsim is the vector of $\xi_{n,M}^{\pm(b)}$. alpha is the signiﬁcance level.
# Output: a logical value, TRUE reject, FALSE accept.
XIMtestT<-function(XIMstat, XIMsim, alpha = 0.05){
  B <- length(XIMsim)
  coef.test <- (1+sum(XIMstat<=XIMsim))/(1+B)
  return(coef.test <= alpha)
}


# This function performs the simulation based test directly.
# Input: xvec, yvec are two vectors, M is the number of right nearest neighbors. B is the number of simulation. alpha is the signiﬁcance level.
# Output: a logical value, TRUE reject, FALSE accept.
XIMtest<-function(xvec, yvec, M, B, alpha = 0.05){
  coef.sim <- XIMsim(n,M,B)
  coef.stat <- XIMstat(xvec,yvec,M)
  return(XIMtestT(coef.stat, coef.sim, alpha))
}





