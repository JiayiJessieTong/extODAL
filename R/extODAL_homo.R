#' @title Run extODAL simulation -- homogeneous (i.e., replication of ODAL)
#'
#' @param Nsim an integer, total number of iterations
#' @param setting setting of the simulation, ("A" or "B")
#' @param parallel_run if the simulation run in parallel (default is "FALSE")
#' @param plotit if a plot will be made (default is "FALSE")
#'
#' @return MSE of three methods (Pooled, local, and ODAL)
#' @import matlib parallel survival
#' @importFrom stats binomial glm optim rbinom rnorm runif
#' @importFrom grDevices rgb
#' @importFrom graphics axis legend lines
#' @export

extODAL_homo <- function(Nsim, setting,
                         parallel_run = FALSE,
                         plotit = FALSE){
  set.seed(4321)
  beta_true = c(-1,1,-1,1,-1)

  # functions to be used
  data_generator <-  function(N,beta,K,n){

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

  ODAL_homo <-function(Data, beta_true){
    N = dim(Data)[1]
    p = dim(Data)[2]- 2
    K = length(unique(Data$ID))

    ########### METHOD 1: pooled analysis ##############
    # overall data
    X = as.matrix(Data[,-c(5,6)])
    Y = Data$Y

    # pooled beta and var
    fit_pooled = summary(glm(Y~X, family = "binomial"(link = "logit")))
    est_pooled = fit_pooled$coefficients[,1]
    sd_pooled = fit_pooled$coefficients[,2]
    var_pooled =sd_pooled^2
    ##### ------------------------------------ #####

    ########### METHOD 2: local analysis ##############
    # extract the local data. Treat ID = 1 as local site
    local_ind = which(Data$ID == 1)
    local_data = Data[local_ind,]
    Xlocal = as.matrix(local_data[,-c(5,6)])
    Ylocal = local_data$Y

    # Local beta and var
    fit_local = summary(glm(Ylocal~Xlocal, family = "binomial"(link = "logit")))
    est_local = fit_local$coefficients[,1]
    sd_local = fit_local$coefficients[,2]
    var_local = sd_local^2
    ##### ------------------------------------ #####

    ########### METHOD 3: ODAL method ##############

    ######### Function 1: expit ########
    ## Input: x
    ## Output: expit(x)
    expit = function(x){exp(x)/(1+exp(x))}
    ##### ------------------------ #####


    ######### Function 2: likelihood function ########
    ## Input: X, n*p matrix (n: total number of patients, p: dimension of covariates)
    ##        Y, n*1 binary vector (1 indicates case and 0 indicates control)
    ##        beta, p*1 vector (coefficients of the covariates)
    ## Output: log-likelihood value
    Lik = function(beta,X,Y){
      design = cbind(1,X)
      sum(Y*(design%*%t(t(beta)))-log(1+exp(design%*%t(t(beta)))))/length(Y)
    }
    ##### ------------------------------------ #####


    ######### Function 3: 1st order gradient of log-likelihood ########
    ## Input: X, n*p matrix (n: total number of patients. p: dimension of covariates)
    ##        Y, n*1 binary vector (1 indicates case and 0 indicates control)
    ##        beta, p*1 vector (coefficients of the covariates)
    ## Output: 1st order gradient of log-likelihood
    Lgradient = function(beta,X,Y){
      design = cbind(1,X)
      t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
    }
    ##### ------------------------------------ #####


    ######### Function 4: 2nd order gradient of log-likelihood ########
    ## Input: X, n*p matrix (n: total number of patients. p: dimension of covariates)
    ##        beta, p*1 vector (coefficients of the covariates)
    ## Output: 2nd order gradient of log-likelihood
    Lgradient2 =function(beta,X){
      design = cbind(1,X)
      Z=expit(design%*%beta)
      t(c(-1*Z*(1-Z))*design)%*%design/nrow(X)
    }
    ##### ------------------------------------ #####


    ######### Function 5: meat of the sandwich variance estimator ########
    ## Input: X, n*p matrix (n: total number of patients. p: dimension of covariates)
    ## Output: meat of the sandwich variance estimator
    Lgradient_meat = function(beta,X,Y){
      design = cbind(1,X)
      Z = expit(design%*%beta)
      t(c(Y-Z)*design)%*%(c(Y-Z)*design)/nrow(X)
    }
    ##### ------------------------------------ #####

    ######### Function 6: surrogate likelihood function ########
    ## Input: beta, p*1 vector (coefficients of the covariates)
    ##       Xlocal: local covariates
    ##       Ylocal: local outcome
    ## Output: surrogate likelihood value
    SL = function(beta){
      -Lik(beta,Xlocal,Ylocal) - L%*%beta
    }
    ##### ------------------------------------ #####

    ######### Function 7: sandwich var est ########
    ## Input:  beta, X, Y, N
    ## Output: sandwich variance estimate of est_ODAL
    Sandwich <- function(beta,X,Y,N){
      mat_L1 = Lgradient_meat(beta,X,Y)
      mat_L2 = Lgradient2(beta,X)
      inv_L2 = solve.default(mat_L2)
      out = inv_L2%*%mat_L1%*%inv_L2/N
      return(out)
    }
    ##### ------------------------------------ #####


    # use local point estimate as initial value
    L = Lgradient(est_local,X,Y)

    # ODAL point estimate
    est_ODAL = optim(est_local,SL,control = list(maxit = 10000,reltol = 1e-10))$par

    # var of the estimator
    cov_ODAL = Sandwich(est_ODAL,Xlocal,Ylocal,N)
    var_ODAL = diag(cov_ODAL)
    ##### ------------------------------------ #####

    bias_pooled = (est_pooled - beta_true)^2
    MSE_pooled = bias_pooled + var_pooled

    bias_local = (est_local - beta_true)^2
    MSE_local = bias_local + var_local

    bias_ODAL = (est_ODAL - beta_true)^2
    MSE_ODAL = bias_ODAL + var_ODAL


    return(list(MSE_pooled = MSE_pooled,
                MSE_local = MSE_local,
                MSE_ODAL = MSE_ODAL))
  }

  main_run_once <- function(N, beta_true, K, n){
    tryCatch(
      {# generate data
        data = data_generator(N, beta_true, K, n)

        # run ODAL
        out = ODAL_homo(data, beta_true)

        return(out)
      },error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
      })
  }

  Nsim = Nsim

  if (parallel_run == FALSE){
    if (setting == "A")
    {
      ### ------ Replication of Setting A ------ ###
      ### fix the number of sites ###
      ### sites have the same number of patients
      ### increase the number of patients in each site ###
      K_1 = 10 # total number of sites
      n_1_list = seq(100,1000,by=100)  # number in one site
      N_1_list = K_1 * n_1_list # overall patients across K_1 sites

      MSE_result_pooled = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
      for (i in 1:length(N_1_list)){
        result_all = replicate(n = Nsim,
                               main_run_once(N = N_1_list[i],
                                             beta_true = beta_true,
                                             K = K_1,
                                             n = n_1_list[i]))
        # result for Nsim iteration
        MSE_result_pooled[i] = mean(apply(matrix(unlist(result_all[1,]), ncol = 5, nrow = Nsim), 2, mean))
        MSE_result_local[i] = mean(apply(matrix(unlist(result_all[2,]), ncol = 5, nrow = Nsim), 2, mean))
        MSE_result_ODAL[i] = mean(apply(matrix(unlist(result_all[3,]), ncol = 5, nrow = Nsim), 2, mean))
      }

      result = as.data.frame(rbind(MSE_result_pooled, MSE_result_local, MSE_result_ODAL))


      if (plotit){
        plot(N_1_list, MSE_result_local,
             type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
             xlab="Total number of patients across 10 sites", ylab="Mean MSE of four covariates", ylim = c(0, max(MSE_result_local)),
             xaxt='n', main = "Replication A: fixed site number, increasing site size")
        axis(side = 1,at =N_1_list,
             labels=N_1_list,lwd.ticks = TRUE)
        lines(N_1_list, MSE_result_pooled, lty = 2,
              type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
        lines(N_1_list, MSE_result_ODAL, lty = 3,
              type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
        legend(6000, max(MSE_result_local), legend = c("Local", "Pooled", "ODAL"),
               lwd=2, lty = c(1,2,3), col=c(rgb(26/255,133/255,172/255,0.5),
                                            rgb(26/255,133/255,0/255,0.5),
                                            rgb(255/255,101/255,80/255,0.5)),
               pch=c(15, 17, 19), bty = "n")
      }
      return(result)
      ##### ------------------------------------ #####

    } else if (setting == "B")
    {

      ### ------ Replication of Setting B ------ ###
      ### fix patients in each site ###
      ### sites have the same number of patients
      ### increase the number of patients in each site ###
      K_1_list = c(2, 5, 10, 15, 30, 40, 50, 75, 90, 100) # total number of sites
      n_1 = 1000  # number in one site
      N_1_list = K_1_list * n_1 # overall patients across K_1 sites

      MSE_result_pooled = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
      for (i in 1:length(N_1_list)){
        result_all = replicate(n = Nsim,
                               main_run_once(N = N_1_list[i],
                                             beta_true = beta_true,
                                             K = K_1_list[i],
                                             n = n_1))
        # result for Nsim iteration
        MSE_result_pooled[i] = mean(apply(matrix(unlist(result_all[1,]), ncol = 5, nrow = Nsim), 2, mean))
        MSE_result_local[i] = mean(apply(matrix(unlist(result_all[2,]), ncol = 5, nrow = Nsim), 2, mean))
        MSE_result_ODAL[i] = mean(apply(matrix(unlist(result_all[3,]), ncol = 5, nrow = Nsim), 2, mean))
      }
      result = as.data.frame(rbind(MSE_result_pooled, MSE_result_local, MSE_result_ODAL))


      if (plotit){
        plot(K_1_list, MSE_result_local,
             type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
             xlab="Total number of sites", ylab="MSE", ylim = c(0, max(MSE_result_local)),
             xaxt='n', main = "Replication B: fixed site size, increasing site number")
        axis(side = 1,at =K_1_list,
             labels=K_1_list,lwd.ticks = TRUE)
        lines(K_1_list, MSE_result_pooled, lty = 2,
              type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
        lines(K_1_list, MSE_result_ODAL, lty = 3,
              type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
      }
      return(result)
      ##### ------------------------------------ #####

    }else{
      stop("Error: not valid setting input (setting is only A or B)")
    }
  }else if (parallel_run == TRUE){ #### parallelization

    if (setting == "A")
    {
      main_run_once_parallel_A <- function(N_n_list, beta_true, K){
        tryCatch(
          {# generate data
            data = data_generator(N_n_list[[1]], beta_true, K, N_n_list[[2]])

            # run ODAL
            out = ODAL_homo(data, beta_true)

            return(out)
          },error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
          })
      }
      ### ------ Replication of Setting A ------ ###
      ### fix the number of sites ###
      ### sites have the same number of patients
      ### increase the number of patients in each site ###
      K_1 = 10 # total number of sites
      n_1_list = seq(100,1000,by=100)  # number in one site
      N_1_list = K_1 * n_1_list # overall patients across K_1 sites

      N_n_list = c()
      for (i in 1:length(N_1_list)){
        N_n_list[[i]] = list(N_1_list[i], n_1_list[i])
      }
      out = mclapply(N_n_list,
                     main_run_once_parallel_A,
                     beta_true = beta_true,
                     K = K_1,
                     mc.cores = detectCores()/2)

      MSE_result_pooled = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
      for (i in 1:length(N_1_list)){
        MSE_result_pooled[i] = mean(out[[i]]$MSE_pooled)
        MSE_result_local[i] = mean(out[[i]]$MSE_local)
        MSE_result_ODAL[i] = mean(out[[i]]$MSE_ODAL)
      }

      result = as.data.frame(rbind(MSE_result_pooled, MSE_result_local, MSE_result_ODAL))


      if (plotit){
        plot(N_1_list, MSE_result_local,
             type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
             xlab="Total number of patients across 10 sites", ylab="Mean MSE of four covariates", ylim = c(0, max(MSE_result_local)),
             xaxt='n', main = "Replication A: fixed site number, increasing site size")
        axis(side = 1,at =N_1_list,
             labels=N_1_list,lwd.ticks = TRUE)
        lines(N_1_list, MSE_result_pooled, lty = 2,
              type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
        lines(N_1_list, MSE_result_ODAL, lty = 3,
              type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
        legend(6000, max(MSE_result_local), legend = c("Local", "Pooled", "ODAL"),
               lwd=2, lty = c(1,2,3), col=c(rgb(26/255,133/255,172/255,0.5),
                                            rgb(26/255,133/255,0/255,0.5),
                                            rgb(255/255,101/255,80/255,0.5)),
               pch=c(15, 17, 19), bty = "n")
      }
      return(result)
      ##### ------------------------------------ #####

    } else if (setting == "B")
    {
      main_run_once_parallel_B <- function(N_K_list, beta_true, n){
        tryCatch(
          {# generate data
            data = data_generator(N_K_list[[1]], beta_true, N_K_list[[2]], n)

            # run ODAL
            out = ODAL_homo(data, beta_true)

            return(out)
          },error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
          })
      }

      ### ------ Replication of Setting B ------ ###
      ### fix patients in each site ###
      ### sites have the same number of patients
      ### increase the number of patients in each site ###
      K_1_list = c(2, 5, 10, 15, 30, 40, 50, 75, 90, 100) # total number of sites
      n_1 = 1000  # number in one site
      N_1_list = K_1_list * n_1 # overall patients across K_1 sites

      N_K_list = c()
      for (i in 1:length(N_1_list)){
        N_K_list[[i]] = list(N_1_list[i], K_1_list[i])
      }
      out = mclapply(N_K_list,
                     main_run_once_parallel_B,
                     beta_true = beta_true,
                     n = n_1,
                     mc.cores = detectCores()/2)

      MSE_result_pooled = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
      for (i in 1:length(N_1_list)){
        MSE_result_pooled[i] = mean(out[[i]]$MSE_pooled)
        MSE_result_local[i] = mean(out[[i]]$MSE_local)
        MSE_result_ODAL[i] = mean(out[[i]]$MSE_ODAL)
      }

      result = as.data.frame(rbind(MSE_result_pooled, MSE_result_local, MSE_result_ODAL))


      if (plotit){
        plot(K_1_list, MSE_result_local,
             type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
             xlab="Total number of sites", ylab="MSE", ylim = c(0, max(MSE_result_local)),
             xaxt='n', main = "Replication B: fixed site size, increasing site number")
        axis(side = 1,at =K_1_list,
             labels=K_1_list,lwd.ticks = TRUE)
        lines(K_1_list, MSE_result_pooled, lty = 2,
              type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
        lines(K_1_list, MSE_result_ODAL, lty = 3,
              type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
      }
      return(result)
      ##### ------------------------------------ #####


    }else{
      stop("Error: not valid setting input (setting is only A or B)")
    }
  }
}
