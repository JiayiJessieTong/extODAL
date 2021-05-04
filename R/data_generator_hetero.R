#' @title Heterogeneous data generation for extODAL
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
data_generator_hetero <-  function(N,beta,K,n)
{
  ######### Function 1: expit ########
  ## Input: x
  ## Output: expit(x)
  expit = function(x){exp(x)/(1+exp(x))}
  ##### ------------------------ #####

  # prevalence median
  pre_med = beta[1]
  pre_list = pre_med + runif(K,-2,0)

  meanY = X1 = X2 = X3 = X4 = rep(NA, N)
  for (i in 1:K){
    ind = seq((i-1)*n + 1, i*n, by = 1)
    X1[ind] = rnorm(n)
    X2[ind] = runif(n,0,1)
    X3[ind] = rbinom(n,1,0.35)
    X4[ind] = rbinom(n,1,0.6)
    X_tmp = cbind(1,X1[ind],X2[ind],X3[ind],X4[ind])
    meanY[ind] = expit(X_tmp%*%c(pre_list[i], beta[-1]))
  }

  Y = rbinom(N,1,meanY)
  data = data.frame(cbind(X1,X2,X3,X4,Y))
  # add hospital ID
  data$ID = rep(1:K, each = n)

  return(data)
}

