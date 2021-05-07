#' @title Run extODAL simulation -- homogeneous (i.e., replication of ODAL)
#'
#' @param Nsim an integer, total number of iterations
#' @param setting setting of the extended ODAL simulation, ("1" or "2")
#' @param parallel_run if the simulation run in parallel (default is "FALSE")
#' @param plotit if a plot will be made (default is "TRUE")
#'
#' @return MSE of three methods (Pooled, local, and ODAL)
#' @import matlib parallel survival
#' @importFrom stats binomial glm optim rbinom rnorm runif
#' @importFrom grDevices rgb
#' @importFrom graphics axis legend lines
#' @export

extODAL_hetero <- function(Nsim, setting,
                           parallel_run = FALSE,
                           plotit = TRUE){
  set.seed(4321)
  beta_true = c(-1,1,-1,1,-1)

  # functions to be used
  Nsim = Nsim

  if (is.infinite(Nsim) | Nsim%%1!=0){
    stop("Error: Nsim should be a finite integer")
  }
  else
  {
    if (parallel_run == FALSE){
      if (setting == "1")
      {
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
        ### -------------- Extension 1 -------------- ###
        ### different size sites
        K_1 = 10 # total number of sites
        n_1_list = matrix(NA, ncol = K_1, nrow = 10)
        for (i in 1:K_1){
          tmp = sample(seq(200,1000,1),10)
          n_1_list[i,] = tmp # K_1 different size sites
        }
        N_1_list = apply(n_1_list, 1, sum) # overall patients across K_1 sites (10 settings)

        MSE_result_pooled = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
        for (i in 1:length(N_1_list)){
          result_all = replicate(n = Nsim,
                                 main_run_once(N = N_1_list[i],
                                               beta_true = beta_true,
                                               K = K_1,
                                               n = n_1_list[i,]))
          # result for Nsim iteration
          MSE_result_pooled[i] = mean(apply(matrix(unlist(result_all[1,]), ncol = 5, nrow = Nsim), 2, mean))
          MSE_result_local[i] = mean(apply(matrix(unlist(result_all[2,]), ncol = 5, nrow = Nsim), 2, mean))
          MSE_result_ODAL[i] = mean(apply(matrix(unlist(result_all[3,]), ncol = 5, nrow = Nsim), 2, mean))
        }

        result = as.data.frame(rbind(MSE_result_pooled, MSE_result_local, MSE_result_ODAL))


        if (plotit){
          order_index = order(N_1_list)
          plot(N_1_list[order_index], MSE_result_local[order_index],
               type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
               xlab="Total number of patients across 10 sites", ylab="Mean MSE of four covariates", ylim = c(0, max(MSE_result_local)),
               xaxt='n', main = "Extension 1: 10 sites with different sizes")
          axis(side = 1,at =N_1_list[order_index],
               labels=N_1_list[order_index],lwd.ticks = TRUE)
          lines(N_1_list[order_index], MSE_result_pooled[order_index], lty = 2,
                type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
          lines(N_1_list[order_index], MSE_result_ODAL[order_index], lty = 3,
                type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
          legend(6000, max(MSE_result_local), legend = c("Local", "Pooled", "ODAL"),
                 lwd=2, lty = c(1,2,3), col=c(rgb(26/255,133/255,172/255,0.5),
                                              rgb(26/255,133/255,0/255,0.5),
                                              rgb(255/255,101/255,80/255,0.5)),
                 pch=c(15, 17, 19), bty = "n")
        }
        return(result)
        ##### ------------------------------------ #####

      } else if (setting == "2")
      {

        main_run_once_hetero <- function(N, beta_true, K, n){
          tryCatch(
            {# generate data
              data = data_generator_hetero(N, beta_true, K, n)

              # run ODAL
              out = ODAL_hetero(data, beta_true)

              return(out)
            },error=function(e){
              cat("ERROR :",conditionMessage(e), "\n")
            })
        }
        ### -------------- Extension 2 -------------- ###
        ### heterogeneous covariates
        K_1 = 10 # total number of sites
        n_1_list = seq(100,1000,by=100)  # number in one site
        N_1_list = K_1 * n_1_list # overall patients across K_1 sites

        Bias_result_pooled = Bias_result_clogit = Bias_result_local = Bias_result_ODAL = rep(0, length(N_1_list))
        for (i in 1:length(N_1_list)){
          result_all = replicate(n = Nsim,
                                 main_run_once_hetero(N = N_1_list[i],
                                                      beta_true = beta_true,
                                                      K = K_1,
                                                      n = n_1_list[i]))
          # result for Nsim iteration
          Bias_result_pooled[i] = mean(apply(matrix(unlist(result_all[1,]), ncol = 5, nrow = Nsim), 2, mean))
          Bias_result_clogit[i] = mean(apply(matrix(unlist(result_all[2,]), ncol = 4, nrow = Nsim), 2, mean))
          Bias_result_local[i] = mean(apply(matrix(unlist(result_all[2,]), ncol = 5, nrow = Nsim), 2, mean))
          Bias_result_ODAL[i] = mean(apply(matrix(unlist(result_all[3,]), ncol = 5, nrow = Nsim), 2, mean))
        }

        result = as.data.frame(rbind(Bias_result_pooled, Bias_result_clogit, Bias_result_local, Bias_result_ODAL))


        if (plotit){
          plot(N_1_list, Bias_result_local,
               type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
               xlab="Total number of patients across 10 sites", ylab="Mean Bias of four covariates", ylim = c(0, max(Bias_result_pooled)),
               xaxt='n', main = "Extension 2: heterogenous prevalence across fixed K sites")
          axis(side = 1,at =N_1_list,
               labels=N_1_list,lwd.ticks = TRUE)
          lines(N_1_list, Bias_result_pooled, lty = 2,
                type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
          lines(N_1_list, Bias_result_ODAL, lty = 3,
                type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
          lines(N_1_list, Bias_result_clogit, lty = 4,
                type="b", col=rgb(172/255,85/255,255/255,0.5), lwd=2, pch=18)
          legend(7000, max(Bias_result_pooled), legend = c("Local", "Pooled","ODAL","clogit"),
                 lwd=2, lty = c(1,2,3,4),
                 col=c(rgb(26/255,133/255,172/255,0.5),
                       rgb(26/255,133/255,0/255,0.5),
                       rgb(255/255,101/255,80/255,0.5),
                       rgb(172/255,85/255,255/255,0.5)),
                 pch=c(15, 17, 19, 18), bty = "n")
        }
        return(result)
        ##### ------------------------------------ #####

      }else{
        stop("Error: not valid setting input (setting is only 1 or 2)")
      }
    }else if (parallel_run == TRUE){ #### parallelization

      if (setting == "1")
      {

        main_run_once_parallel_1 <- function(N_n_list, beta_true, K){
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
        ### -------------- Extension 1 -------------- ###
        ### different size sites
        K_1 = 10 # total number of sites
        n_1_list = matrix(NA, ncol = K_1, nrow = 10)
        for (i in 1:K_1){
          tmp = sample(seq(200,1000,1),10)
          n_1_list[i,] = tmp # K_1 different size sites
        }
        N_1_list = apply(n_1_list, 1, sum) # overall patients across K_1 sites (10 settings)

        N_n_list = c()
        for (i in 1:length(N_1_list)){
          N_n_list[[i]] = list(N_1_list[i], n_1_list[i,])
        }
        if (Sys.info()[1] == "Windows"){
          run_once <- function(input = NULL){
            cl = makeCluster(detectCores()/2)
            out = parLapply(cl, N_n_list,
                            main_run_once_parallel_1,
                            beta_true = beta_true,
                            K = K_1)
            stopCluster(cl)
          }

          cl = makeCluster(detectCores()/2)
          out = parLapply(cl, 1:Nsim, run_once)
          stopCluster(cl)

        } else {

          run_once <- function(input = NULL){
            out = mclapply(N_n_list,
                           main_run_once_parallel_1,
                           beta_true = beta_true,
                           K = K_1,
                           mc.cores = detectCores()/2)

            return(out)
          }

          out = mclapply(1:Nsim, run_once, mc.cores = detectCores()/2)

        }

        MSE_result_pooled = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
        for (i in 1:length(N_1_list)){
          for (iter in Nsim){
            MSE_result_pooled[i] = MSE_result_pooled[i] + mean(out[[iter]][[i]]$MSE_pooled)
            MSE_result_local[i] = MSE_result_local[i] + mean(out[[iter]][[i]]$MSE_local)
            MSE_result_ODAL[i] = MSE_result_ODAL[i] + mean(out[[iter]][[i]]$MSE_ODAL)
          }
        }

        result = as.data.frame(rbind(MSE_result_pooled, MSE_result_local, MSE_result_ODAL))

        if (plotit){
          order_index = order(N_1_list)
          plot(N_1_list[order_index], MSE_result_local[order_index],
               type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
               xlab="Total number of patients across 10 sites", ylab="Mean MSE of four covariates", ylim = c(0, max(MSE_result_local)),
               xaxt='n', main = "Extension 1: 10 sites with different sizes")
          axis(side = 1,at =N_1_list[order_index],
               labels=N_1_list[order_index],lwd.ticks = TRUE)
          lines(N_1_list[order_index], MSE_result_pooled[order_index], lty = 2,
                type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
          lines(N_1_list[order_index], MSE_result_ODAL[order_index], lty = 3,
                type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
          legend(6000, max(MSE_result_local), legend = c("Local", "Pooled", "ODAL"),
                 lwd=2, lty = c(1,2,3), col=c(rgb(26/255,133/255,172/255,0.5),
                                              rgb(26/255,133/255,0/255,0.5),
                                              rgb(255/255,101/255,80/255,0.5)),
                 pch=c(15, 17, 19), bty = "n")
        }
        return(result)

      } else if (setting == "2")
      {
        main_run_once_parallel_2 <- function(N_n_list, beta_true, K){
          print("here")
          tryCatch(
            {# generate data
              data = data_generator_hetero(N_n_list[[1]], beta_true, K, N_n_list[[2]])

              # run ODAL
              out = ODAL_hetero(data, beta_true)

              return(out)
            },error=function(e){
              cat("ERROR :",conditionMessage(e), "\n")
            })
        }

        ### -------------- Extension 2 -------------- ###
        ### heterogeneous covariates
        K_1 = 10 # total number of sites
        n_1_list = seq(100,1000,by=100)  # number in one site
        N_1_list = K_1 * n_1_list # overall patients across K_1 sites

        N_n_list = c()
        for (i in 1:length(N_1_list)){
          N_n_list[[i]] = list(N_1_list[i], n_1_list[i])
        }

        if (Sys.info()[1] == "Windows"){
          run_once <- function(input = NULL){
            cl = makeCluster(detectCores()/2)
            out = parLapply(cl, N_n_list,
                            main_run_once_parallel_2,
                            beta_true = beta_true,
                            K = K_1)
            stopCluster(cl)
          }

          cl = makeCluster(detectCores()/2)
          out = parLapply(cl, 1:Nsim, run_once)
          stopCluster(cl)

        } else {

          run_once <- function(input = NULL){
            out = mclapply(N_n_list,
                           main_run_once_parallel_2,
                           beta_true = beta_true,
                           K = K_1,
                           mc.cores = detectCores()/2)

            return(out)
          }

          out = mclapply(1:Nsim, run_once, mc.cores = detectCores()/2)

        }

        Bias_result_pooled = Bias_result_clogit  = Bias_result_local = Bias_result_ODAL = rep(0, length(N_1_list))
        for (i in 1:length(N_1_list)){
          for (iter in Nsim){
            Bias_result_pooled[i] = Bias_result_pooled[i] + mean(out[[iter]][[i]]$bias_pooled)
            Bias_result_clogit[i] = Bias_result_clogit[i] + mean(out[[iter]][[i]]$bias_clogit)
            Bias_result_local[i] = Bias_result_local[i] + mean(out[[iter]][[i]]$bias_local)
            Bias_result_ODAL[i] =  Bias_result_ODAL[i] + mean(out[[iter]][[i]]$bias_ODAL)
          }
        }


        result = as.data.frame(rbind(Bias_result_pooled, Bias_result_clogit,
                                     Bias_result_local, Bias_result_ODAL))


        if (plotit){
          plot(N_1_list, Bias_result_local,
               type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
               xlab="Total number of patients across 10 sites", ylab="Mean Bias of four covariates",
               ylim = c(0, max(Bias_result_pooled)),
               xaxt='n', main = "Extension 2: heterogenous prevalence across fixed K sites")
          axis(side = 1,at =N_1_list,
               labels=N_1_list,lwd.ticks = TRUE)
          lines(N_1_list, Bias_result_pooled, lty = 2,
                type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
          lines(N_1_list, Bias_result_ODAL, lty = 3,
                type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
          lines(N_1_list, Bias_result_clogit, lty = 4,
                type="b", col=rgb(172/255,85/255,255/255,0.5), lwd=2, pch=18)
          legend(7000, max(Bias_result_pooled), legend = c("Local", "Pooled","ODAL","clogit"),
                 lwd=2, lty = c(1,2,3,4),
                 col=c(rgb(26/255,133/255,172/255,0.5),
                       rgb(26/255,133/255,0/255,0.5),
                       rgb(255/255,101/255,80/255,0.5),
                       rgb(172/255,85/255,255/255,0.5)),
                 pch=c(15, 17, 19, 18), bty = "n")
        }
        return(result)
        ##### ------------------------------------ #####


      }else{
        stop("Error: not valid setting input (setting is only A or B)")
      }
    }
  }


}
