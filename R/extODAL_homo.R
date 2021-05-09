#' @title Run extODAL simulation -- homogeneous (i.e., replication of ODAL)
#'
#' @param Nsim an integer, total number of iterations
#' @param setting setting of the simulation, ("A" or "B")
#' @param parallel_run if the simulation run in parallel (default is "FALSE")
#' @param plotit if a plot will be made (default is "TRUE")
#'
#' @return mean MSEs of covariates for three methods (Pooled, local, and ODAL)
#' @import matlib parallel survival
#' @importFrom stats binomial glm optim rbinom rnorm runif
#' @importFrom grDevices rgb
#' @importFrom graphics axis legend lines
#' @export

extODAL_homo <- function(Nsim, setting,
                         parallel_run = FALSE,
                         plotit = TRUE){
  set.seed(4321)
  beta_true = c(-1,1,-1,1,-1)

  # functions to be used
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

  if (is.infinite(Nsim) | Nsim%%1!=0){
    stop("Error: Nsim should be a finite integer")
  }
  else
  {
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
          MSE_result_pooled[i] = mean(unlist(result_all[1,]))
          MSE_result_local[i] = mean(unlist(result_all[2,]))
          MSE_result_ODAL[i] = mean(unlist(result_all[3,]))
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
          MSE_result_pooled[i] = mean(unlist(result_all[1,]))
          MSE_result_local[i] = mean(unlist(result_all[2,]))
          MSE_result_ODAL[i] = mean(unlist(result_all[3,]))
        }
        result = as.data.frame(rbind(MSE_result_pooled, MSE_result_local, MSE_result_ODAL))


        if (plotit){
          plot(K_1_list, MSE_result_local,
               type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
               xlab="Total number of sites", ylab="Mean MSE of four covariates", ylim = c(0, max(MSE_result_local)),
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


        cat("Using parLapply to run in parallel -- Setting A")
        cl = makeCluster(detectCores()/2)
        out = parLapply(cl, rep(N_n_list, Nsim), main_run_once_parallel_A,
                        beta_true = beta_true,
                        K = K_1)
        stopCluster(cl)


        MSE_result_pooled = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
        for (i in 1:length(N_1_list)){
          for (iter in Nsim){
            MSE_result_pooled[i] = MSE_result_pooled[i] + out[[(iter-1) * length(N_1_list) + i]]$MSE_pooled/Nsim
            MSE_result_local[i] = MSE_result_local[i] + out[[(iter-1) * length(N_1_list) + i]]$MSE_local/Nsim
            MSE_result_ODAL[i] = MSE_result_ODAL[i] + out[[(iter-1) * length(N_1_list) + i]]$MSE_ODAL/Nsim
          }
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

        cat("Using parLapply to run in parallel -- Setting B")
        cl = makeCluster(detectCores()/2)
        out = parLapply(cl, rep(N_K_list, Nsim), main_run_once_parallel_B,
                        beta_true = beta_true,
                        n = n_1)
        stopCluster(cl)


        MSE_result_pooled = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
        for (i in 1:length(N_1_list)){
          for (iter in Nsim){
            MSE_result_pooled[i] = MSE_result_pooled[i] + out[[(iter-1) * length(N_1_list) + i]]$MSE_pooled/Nsim
            MSE_result_local[i] = MSE_result_local[i] + out[[(iter-1) * length(N_1_list) + i]]$MSE_local/Nsim
            MSE_result_ODAL[i] = MSE_result_ODAL[i] + out[[(iter-1) * length(N_1_list) + i]]$MSE_ODAL/Nsim
          }
        }

        result = as.data.frame(rbind(MSE_result_pooled, MSE_result_local, MSE_result_ODAL))


        if (plotit){
          plot(K_1_list, MSE_result_local,
               type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
               xlab="Total number of sites", ylab="Mean MSE of four covariates", ylim = c(0, max(MSE_result_local)),
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
}
