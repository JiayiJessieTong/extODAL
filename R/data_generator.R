#' @title Data generation for extODAL
#'
#' @param N total number of observations
#' @param beta p*1 vector (coefficients of the covariates)
#' @param K total number of sites
#' @param n number of patients in each site
#'
#' @return data
#' @import matlib parallel survival
#' @importFrom stats binomial glm optim rbinom rnorm runif
#' @importFrom grDevices rgb
#' @importFrom graphics axis legend lines
#' @export
data_generator <-  function(N,beta,K,n)
{
  ######### Function 1: expit ########
  ## Input: x
  ## Output: expit(x)
  expit = function(x){exp(x)/(1+exp(x))}
  ##### ------------------------ #####

  X1 = rnorm(N)
  X2 = runif(N,0,1)
  X3 = rbinom(N,1,0.35)
  X4 = rbinom(N,1,0.6)
  X = cbind(1,X1,X2,X3,X4)
  meanY = expit(X%*%beta)
  Y = rbinom(N,1,meanY)
  data = data.frame(cbind(X1,X2,X3,X4,Y))
  # add hospital ID
  if (length(n) == 1){
    data$ID = rep(1:K, each = n)
  }else{
    data$ID = rep(1:K, n)
  }
  return(data)
}

