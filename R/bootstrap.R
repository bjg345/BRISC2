BRISC2_bootstrap <- function(BRISC2_Out, n_boot = 100, h = 1, n_omp = 1, init = "Initial", verbose = TRUE, nugget_status = 1){

  if(missing(BRISC2_Out)){stop("error: BRISC2_bootstrap expects BRISC2_Out\n")}

  if(nugget_status == 0){fix_nugget = 0}
  if(nugget_status == 1){fix_nugget = 1}
  X <- BRISC2_Out$X
  n.omp.threads <- as.integer(n_omp)
  n.neighbors <- BRISC2_Out$n.neighbors
  eps <- BRISC2_Out$eps
  cov.model <- BRISC2_Out$cov.model
  p  <- ncol(X)
  n <- nrow(X)


  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"
  storage.mode(n.neighbors) <- "integer"
  storage.mode(n.omp.threads) <- "integer"
  storage.mode(eps) <- "double"

  cov.model.names <- c("exponential","spherical","matern","gaussian")
  cov.model.indx <- which(cov.model == cov.model.names) - 1
  storage.mode(cov.model.indx) <- "integer"


  cov.model <- BRISC2_Out$cov.model
  norm.residual = BRISC2_Out$BRISC2_Object$norm.residual
  B =  BRISC2_Out$BRISC2_Object$B
  F = BRISC2_Out$BRISC2_Object$F
  Xbeta = BRISC2_Out$BRISC2_Object$Xbeta
  D = BRISC2_Out$BRISC2_Object$D
  d = BRISC2_Out$BRISC2_Object$d
  nnIndx = BRISC2_Out$BRISC2_Object$nnIndx
  nnIndxLU = BRISC2_Out$BRISC2_Object$nnIndxLU
  CIndx = BRISC2_Out$BRISC2_Object$CIndx
  Length.D = BRISC2_Out$BRISC2_Object$Length.D

  if(init == "Initial"){
    if(cov.model == "matern") {theta_boot_init <- c(BRISC2_Out$init[2]/BRISC2_Out$init[1], BRISC2_Out$init[3], BRISC2_Out$init[4])}
    else {theta_boot_init <- c(BRISC2_Out$init[2]/BRISC2_Out$init[1], BRISC2_Out$init[3])}
  }
  if(init == "Estimate"){
    if(cov.model == "matern") {theta_boot_init <- c(BRISC2_Out$Theta[2]/BRISC2_Out$Theta[1], BRISC2_Out$Theta[3], BRISC2_Out$Theta[4])}
    else {theta_boot_init <- c(BRISC2_Out$Theta[2]/BRISC2_Out$Theta[1], BRISC2_Out$Theta[3])}
  }

  theta_boot_init <- sqrt(theta_boot_init)

  p3 <- proc.time()

  if(h > 1){
    cl <- makeCluster(h)
    clusterExport(cl=cl, varlist=c("norm.residual", "X", "B", "F", "Xbeta", "D", "d", "nnIndx", "nnIndxLU",
                                   "CIndx", "n", "p", "n.neighbors", "theta_boot_init", "cov.model.indx", "Length.D",
                                   "n.omp.threads", "bootstrap_brisc", "eps", "fix_nugget"),envir=environment())
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tBootstrap Progress"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
      pboptions(type = "txt", char = "=")
      result <- pblapply(1:n_boot,bootstrap_brisc,norm.residual, X, B, F, Xbeta, D, d, nnIndx, nnIndxLU, CIndx, n, p, n.neighbors, theta_boot_init,
                                            cov.model.indx, Length.D, n.omp.threads, eps, fix_nugget, cl = cl)
      }
    if(verbose != TRUE){result <- parLapply(cl,1:n_boot,bootstrap_brisc,norm.residual, X, B, F, Xbeta, D, d, nnIndx, nnIndxLU, CIndx, n, p, n.neighbors, theta_boot_init,
                                            cov.model.indx, Length.D, n.omp.threads, eps, fix_nugget)}
    stopCluster(cl)
  }
  if(h == 1){
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tBootstrap Progress"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
      pboptions(type = "txt", char = "=")
      result <- pblapply(1:n_boot,bootstrap_brisc,norm.residual, X, B, F, Xbeta, D, d, nnIndx, nnIndxLU, CIndx, n, p, n.neighbors, theta_boot_init,
                       cov.model.indx, Length.D, n.omp.threads, eps, fix_nugget)
    }

    if(verbose != TRUE){
      result <- lapply(1:n_boot,bootstrap_brisc,norm.residual, X, B, F, Xbeta, D, d, nnIndx, nnIndxLU, CIndx, n, p, n.neighbors, theta_boot_init,
                       cov.model.indx, Length.D, n.omp.threads, eps, fix_nugget)
    }
  }

  p4 <- proc.time()

  result_table = arrange(result)

  estimate <- c(BRISC2_Out$Beta, BRISC2_Out$Theta)
  result_CI <- matrix(0,2,length(estimate))

  for(i in 1:length(estimate)){
    result_CI[,i] <- 2*estimate[i] - quantile(result_table[,i], c(.975,.025))
  }

  result_list <- list()


  result_list$boot.Theta <- result_table[,(length(BRISC2_Out$Beta) + 1):dim(result_table)[2]]
  if (cov.model != "matern") {colnames(result_list$boot.Theta) <- c("sigma.sq", "tau.sq", "phi")}
  if (cov.model == "matern") {colnames(result_list$boot.Theta) <- c("sigma.sq", "tau.sq", "phi", "nu")}
  result_list$boot.Beta <- as.matrix(result_table[,1:length(BRISC2_Out$Beta)])
  colnames(result_list$boot.Beta) <- rep(0, length(BRISC2_Out$Beta))
  for(i in 1:length(BRISC2_Out$Beta)){
    name_beta <- paste0("beta_",i)
    colnames(result_list$boot.Beta)[i] <- name_beta
  }
  result_list$confidence.interval <- cbind(result_CI[,1:length(BRISC2_Out$Beta)],pmax(result_CI[,(length(BRISC2_Out$Beta) + 1)
                                     :dim(result_table)[2]], 0*result_CI[,(length(BRISC2_Out$Beta) + 1):dim(result_table)[2]]))
  if (cov.model != "matern")  {colnames(result_list$confidence.interval)[(length(BRISC2_Out$Beta) + 1):dim(result_table)[2]] <-
    c("sigma.sq", "tau.sq", "phi")}
  if (cov.model == "matern")  {colnames(result_list$confidence.interval)[(length(BRISC2_Out$Beta) + 1):dim(result_table)[2]] <-
    c("sigma.sq", "tau.sq", "phi", "nu")}
  for(i in 1:length(BRISC2_Out$Beta)){
    name_beta <- paste0("beta_",i)
    colnames(result_list$confidence.interval)[i] <- name_beta
  }
  result_list$boot.time = p4 - p3
  result_list
}
