simulation<-function(runs_i){
  
  
  set.seed(runs_i)
  sample_i<-sample(1:1e5,2,replace = FALSE)
  library(locpol)
  library(MASS)
  
  library(grpreg)
  library(survcomp)
  library(fda)
  
  
  n1<-2000
  n2<-0+n1
  n<-0+2000
  
  use_n_all<-c(200,300,500,1000,2000)
  
  mfunctions<-300
  mfunctions_true<-2
  
  mfunctions_true_dev<-2
  q<-300
  q_true<-2
  
  
  rho<-0.8
  
  h<-calh<-0.01
  designpoints<-seq(0,1,calh)
  
  t<-1/h+1
  
  numberofbasis<-10
  fourierbasis10 <- create.fourier.basis(rangeval=c(0,1),nbasis=numberofbasis+1,period=1,dropind = 1)
  plot(fourierbasis10)
  
  
  fourierbasis11 <- create.fourier.basis(rangeval=c(0,1),nbasis=numberofbasis+1,period=1)
  plot(fourierbasis11)
  
  
  little_Sigma<-matrix(0,10,10)
  
  for (j in 1:5){
    little_Sigma[j,j]<- (j)^(-2)
  }
  for (j in 6:10){
    little_Sigma[j,j]<-(j-5)^(-2)
  }
  
  
  
  Sigma<-matrix(0,nrow(little_Sigma)*mfunctions,nrow(little_Sigma)*mfunctions)
  
  for(mfunctions_i in 1:mfunctions){
    Sigma[((mfunctions_i-1)*nrow(little_Sigma)+1):(mfunctions_i*nrow(little_Sigma)),((mfunctions_i-1)*nrow(little_Sigma)+1):(mfunctions_i*nrow(little_Sigma))]<-little_Sigma
    
  }
  
  truebeta_list<-list()
  
  truebeta_list[[1]]<-fd(c(0.5, rep(0,9)),fourierbasis10)
  truebeta_list[[2]]<-fd(c(0,-0.3,rep(0,8)),fourierbasis10)
  
  for (turebetalist_i in (mfunctions_true+1):mfunctions) {
    truebeta_list[[turebetalist_i]]<-fd(  rep(0,numberofbasis) ,fourierbasis10)
    
  }
  
  
  
  
  true_beta_flag<-c(rep(1,mfunctions_true), rep(0,mfunctions-mfunctions_true))
  
  gamma<-c(c(0.5, -0.3  ),rep(0,q-q_true))
  
  eta_dev<-c(0,0.005, -0.003 ,rep(0,mfunctions-1-mfunctions_true_dev))
  
  
  RMSE_true_gamma<-gamma
  RMSE_true_eta_dev<-eta_dev
  ISE_true <- matrix(0,nrow=length(designpoints),ncol = mfunctions)
  for(i in 1:mfunctions){
    ISE_true[,i] <- eval.fd(designpoints, truebeta_list[[i]])
  }
  
  
  
  correlationmatrix<-matrix(0,nrow(Sigma),q)
  
  
  SigmaZ<-rho^t(sapply(1:q, function(i, j) abs(i-j), 1:q))
  
  Bigsigma<-rbind(cbind(Sigma, correlationmatrix), cbind(t(correlationmatrix), SigmaZ))
  
  mu<-rep(0,nrow(Bigsigma))
  set.seed(sample_i[1])
  data<-mvrnorm(n,mu,Bigsigma)
  
  
  ximatrix<-list()
  for(ximatrix_i in 1:mfunctions){
    ximatrix[[ximatrix_i]]<-data[,((ximatrix_i-1)*nrow(little_Sigma)+1):(ximatrix_i*nrow(little_Sigma))]
  }
  
  Vxmatrix<-ximatrix
  for(vmatrix_i in 1:mfunctions){
    
    if((vmatrix_i>0)&(vmatrix_i<=(mfunctions-2))){
      Vxmatrix[[vmatrix_i]]<-ximatrix[[vmatrix_i]] +rho*(ximatrix[[vmatrix_i+1]]+ximatrix[[vmatrix_i+2]])
    }
    
    if(vmatrix_i ==(mfunctions-1)){
      Vxmatrix[[vmatrix_i]]<-ximatrix[[vmatrix_i]] +rho*(ximatrix[[vmatrix_i+1]]+ximatrix[[1]])
    }
    if(vmatrix_i ==mfunctions){
      Vxmatrix[[vmatrix_i]]<-ximatrix[[vmatrix_i]] +rho*(ximatrix[[1]]+ximatrix[[2]])
    }
  }
  
  
  
  rawX <- list()
  observeX <- list()
  X <- list()
  dev2 <- list()
  
  dev2_int_true <- matrix(0, n, mfunctions)
  dev2_int_obs <- matrix(0, n, mfunctions)
  
  lengthout <- length(designpoints)
  grid <- designpoints
  delta <- grid[2] - grid[1]
  
  
  
  
  
  trapezoid_integral <- function(x, y) {
    
    
    
    
    
    
    
    
    
    sum(    diff(x) * (y[-length(y)] + y[-1]) / 2  )
  }
  
  for (mfunctions_i in 1:mfunctions) {
    
    rawX[[mfunctions_i]]<- fd(t(Vxmatrix[[mfunctions_i]]),fourierbasis10)
    dev2[[mfunctions_i]]<-deriv.fd( rawX[[mfunctions_i]],2)
    
    for(dev2_i in 1:n){
      dev2fd_i<-fd(dev2[[mfunctions_i]]$coefs[,dev2_i],dev2[[mfunctions_i]]$basis)
      dev2_int_true[dev2_i,mfunctions_i] <-sqrt(inprod(dev2fd_i,dev2fd_i))
    }
    
    
    observeX[[mfunctions_i]]<- t(eval.fd(designpoints, rawX[[mfunctions_i]]))
    
    
    X[[mfunctions_i]] <- matrix( NA_real_,nrow = lengthout,ncol = n)
    
    for (smooth_i in 1:n) {
      
      dataframe <-data.frame( designpoints = designpoints,
                              observeX = observeX[[mfunctions_i]][smooth_i, ])
      
      lpfit <-locpol(observeX ~ designpoints, data = dataframe, xeval = designpoints)
      
      X[[mfunctions_i]][, smooth_i] <-lpfit$lpFit[, 2]
    }
    
    dev2_grid <- designpoints[2:(length(designpoints) - 1)]
    
    dev2_obs_fd <- diff(X[[mfunctions_i]], differences = 2) / delta^2
    
    for (dev_i in seq_len(n)) {
      dev2_int_obs[dev_i, mfunctions_i] <- sqrt(trapezoid_integral(x = dev2_grid, y = dev2_obs_fd[, dev_i]^2 ) )
    }
    
    cat("mfunctions_i:",mfunctions_i,";")
    
  }
  
  
  dim(X[[mfunctions_i]])
  
  
  
  Z<-data[,(nrow(Sigma)+1):(nrow(Sigma)+q)]
  
  
  
  yita<-rep(0,n)
  for(yita_i in 1:mfunctions_true){
    yita<-yita+ inprod(rawX[[yita_i]],truebeta_list[[yita_i]])
  }
  parameter<-exp(yita+Z%*%gamma+dev2_int_true%*%eta_dev)
  
  summary(yita)
  
  summary(Z%*%gamma)
  
  summary(dev2_int_true%*%eta_dev)
  
  
  
  failuretime<-rexp(n, rate = parameter)
  
  
  
  
  
  calculate_estimation_accuracy <- function(
    estimated_coef,
    true_scalar_coef,
    true_eta_dev_coef,
    true_functional_coef,
    fun_Eigfun,
    sn_chose,
    q,
    mfunctions) {
    
    
    estimated_coef <- as.numeric(estimated_coef)
    
    integration_grid <- seq( from = 0,to = 1,length.out = nrow(fun_Eigfun[[1]]))
    
    
    fpc_coef_index <- vector(mode = "list",length = mfunctions)
    
    current_index <- q + mfunctions
    
    for (j in seq_len(mfunctions)) {
      
      fpc_coef_index[[j]] <- current_index + seq_len(sn_chose[j])
      
      current_index <- current_index + sn_chose[j]
    }
    
    
    
    
    
    
    estimated_scalar_coef <-estimated_coef[seq_len(q)]
    
    scalar_RMSE <-sqrt(mean((estimated_scalar_coef -  true_scalar_coef)^2))
    
    
    estimated_dev_coef <-estimated_coef[q+seq_len(mfunctions)]
    dev_RMSE <-sqrt(mean((estimated_dev_coef -  true_eta_dev_coef)^2))
    
    
    
    
    
    functional_ISE <- rep(NA_real_, mfunctions)
    
    estimated_functional_coef <-  matrix(NA_real_, nrow = length(integration_grid),
                                         ncol = mfunctions)
    
    for (j in seq_len(mfunctions)) {
      
      
      estimated_fpc_coef <-estimated_coef[fpc_coef_index[[j]]]
      
      
      estimated_functional_coef[, j] <-as.numeric(fun_Eigfun[[j]] %*%estimated_fpc_coef  )
      
      
      squared_error <- (estimated_functional_coef[, j] -true_functional_coef[, j])^2
      
      
      functional_ISE[j] <- trapezoid_integral(x = integration_grid,y = squared_error)
    }
    
    
    
    
    
    mean_functional_ISE <-mean(functional_ISE)
    
    total_functional_ISE <- sum(functional_ISE)
    
    return(
      list(
        scalar_RMSE = scalar_RMSE,
        dev_RMSE = dev_RMSE,
        functional_ISE = functional_ISE,
        mean_functional_ISE = mean_functional_ISE,
        total_functional_ISE = total_functional_ISE,
        estimated_scalar_coef = estimated_scalar_coef,
        estimated_dev_coef = estimated_dev_coef,
        estimated_functional_coef =estimated_functional_coef
      )
    )
  }
  
  
  
  
  
  
  
  calculate_coxph_accuracy <- function(
    estimated_coef,
    true_scalar_coef,
    true_functional_coef,
    fun_Eigfun,
    sn_chose,
    q,
    mfunctions) {
    
    estimated_coef <- as.numeric(estimated_coef)
    
    
    expected_number_parameters <- q + sum(sn_chose)
    
    if (length(estimated_coef) != expected_number_parameters) {
      stop(
        paste0(
          "coxph coefficient length is ",
          length(estimated_coef),
          ", but expected ",
          expected_number_parameters,
          "."
        )
      )
    }
    
    
    
    
    
    integration_grid <- seq( from = 0,to = 1,length.out = nrow(fun_Eigfun[[1]]))
    cox_fpc_coef_index <- vector(mode = "list",length = mfunctions)
    
    current_index <- q
    
    for (j in seq_len(mfunctions)) {
      
      cox_fpc_coef_index[[j]] <-    current_index + seq_len(sn_chose[j])
      
      current_index <-    current_index + sn_chose[j]
    }
    
    
    
    
    
    
    estimated_scalar_coef <- estimated_coef[seq_len(q)]
    
    scalar_RMSE <-   sqrt(  mean( (  estimated_scalar_coef -   true_scalar_coef  )^2)   )
    
    
    
    
    
    functional_ISE <- rep( NA_real_,mfunctions)
    
    estimated_functional_coef <- matrix( NA_real_,  nrow = length(integration_grid),ncol = mfunctions)
    
    colnames(    estimated_functional_coef  ) <- paste0(    "functional_",    seq_len(mfunctions))
    
    for (j in seq_len(mfunctions)) {
      
      estimated_fpc_coef <-   estimated_coef[ cox_fpc_coef_index[[j]]      ]
      
      eigenfunction_j <- as.matrix(  fun_Eigfun[[j]]    )
      
      if (ncol(eigenfunction_j) != length(estimated_fpc_coef)) {
        stop(
          paste0(
            "Mismatch in functional predictor ",
            j,
            ": number of eigenfunctions and FPC coefficients differ."
          )
        )
      }
      
      if (nrow(eigenfunction_j) != nrow(true_functional_coef)) {
        stop(
          paste0(
            "Mismatch in functional predictor ",
            j,
            ": estimated and true functions are evaluated on different grids."
          )
        )
      }
      
      if (ncol(true_functional_coef) != mfunctions) {
        stop("true_functional_coef must have mfunctions columns.")
      }
      
      
      estimated_functional_coef[, j] <-  as.numeric(   eigenfunction_j %*%  estimated_fpc_coef      )
      
      true_functional_coef_j <-      true_functional_coef[, j]
      
      squared_error_j <-  ( estimated_functional_coef[, j] -   true_functional_coef_j    )^2
      
      functional_ISE[j] <-   trapezoid_integral(  x = integration_grid, y = squared_error_j    )
    }
    
    
    
    
    
    mean_functional_ISE <-  mean(    functional_ISE  )
    
    total_functional_ISE <- sum(  functional_ISE    )
    
    
    return(list(scalar_RMSE = scalar_RMSE,
                functional_ISE = functional_ISE,
                mean_functional_ISE = mean_functional_ISE,
                total_functional_ISE =  total_functional_ISE,
                estimated_scalar_coef =  estimated_scalar_coef,
                estimated_functional_coef = estimated_functional_coef))
  }
  
  
  
  
  
  
  censoringvalues_all<-c( 0.2)
  
  
  tau_all<-c(1.8 )
  
  tau<-0.3
  set.seed(sample_i[2])
  censoringtime<-runif(n,0,tau)
  event<-rep(0,n)
  event<-as.numeric(failuretime<censoringtime)
  sum(event)/n
  
  for (use_n in use_n_all) {
    
    
    
    
    p_group<-c()
    fun_Eigval<-list()
    fun_Eigfun<-list()
    eigenscore<-list()
    sn_chose<-rep(0,mfunctions)
    
    all_designmatrix<-cbind(Z[1:use_n,],dev2_int_obs[1:use_n,])
    
    for (pca_i in 1:mfunctions) {
      use_x_n <- X[[pca_i]][, 1:use_n, drop = FALSE]
      scaleX <- scale(  t(use_x_n),  center = TRUE,  scale = FALSE )
      svdresult <- svd(scaleX)
      Eigvec <- svdresult$v
      Singval <- svdresult$d
      h_fpca <-  1 / (dim(use_x_n)[1] - 1)
      n_subject <- nrow(scaleX)
      est_Eigval <-    Singval^2 * h_fpca / (n_subject - 1)
      total_variation <- sum(Singval^2)
      variance_contributions <- Singval^2 / total_variation
      cumulative_variance <-  cumsum(variance_contributions)
      sn_chose[pca_i] <-which(cumulative_variance >= 0.95)[1]
      fun_Eigfun[[pca_i]] <-Eigvec[,seq_len(sn_chose[pca_i]),drop = FALSE] / sqrt(h_fpca)
      eigenscore[[pca_i]] <-scaleX %*%fun_Eigfun[[pca_i]] *h_fpca
      fun_Eigval[[pca_i]] <-est_Eigval[seq_len(sn_chose[pca_i])]
      all_designmatrix <-cbind(all_designmatrix, eigenscore[[pca_i]])
    }
    
    
    
    
    
    
    
    retained_fpc_each <- sn_chose
    
    
    average_retained_fpc <- mean(retained_fpc_each)
    
    
    total_retained_fpc <- sum(retained_fpc_each)
    
    
    number_scalar_parameters <- q
    
    
    number_derivative_parameters <- mfunctions
    
    
    number_functional_parameters <- total_retained_fpc
    
    
    
    total_regression_parameters <-  number_scalar_parameters +  number_derivative_parameters +
      number_functional_parameters
    
    
    stopifnot(  ncol(all_designmatrix) == total_regression_parameters)
    
    cat(
      "\n----------------------------------------",
      "\nuse_n =", use_n,
      "\nRetained FPCs for each functional predictor =",
      paste(retained_fpc_each, collapse = ", "),
      "\nAverage retained FPCs =", average_retained_fpc,
      "\nTotal retained FPCs =", total_retained_fpc,
      "\nNumber of scalar parameters =", number_scalar_parameters,
      "\nNumber of derivative parameters =", number_derivative_parameters,
      "\nTotal number of regression parameters =",
      total_regression_parameters,
      "\n----------------------------------------\n"
    )
    
    
    
    
    
    p_group <-c(1:q,(1:mfunctions)+q,rep((1:mfunctions)+mfunctions +q, sn_chose ))
    m_fun<-c()
    for(m_i in 1:mfunctions){
      m_fun[m_i]<-sqrt(sum(fun_Eigval[[m_i]]))
    }
    
    m <- c(rep(1,q),rep(1,mfunctions),m_fun)
    
    designmatrix<-all_designmatrix
    for (tau_j in 1:length(tau_all)) {
      
      
      
      
      
      
      censoringvalues<-censoringvalues_all[tau_j]
      tau<- tau_all[tau_j]
      
      
      set.seed(sample_i[2])
      censoringtime<-runif(n,0,tau)
      
      
      event<-rep(0,n)
      event<-as.numeric(failuretime<censoringtime)
      sum(event)/n
      
      time<-failuretime*event+censoringtime*(rep(1,n)-event)
      
      
      
      
      
      
      number_events <- sum(event[1:use_n])
      
      
      event_rate <- number_events / use_n
      
      
      actual_censoring_rate <- 1 - event_rate
      
      
      EPV <- number_events / total_regression_parameters
      
      
      number_coxph_parameters <- q + total_retained_fpc
      EPV_coxph <- number_events / number_coxph_parameters
      
      
      
      cat(
        "\n========================================",
        "\nSimulation run =", runs_i,
        "\nSample size =", use_n,
        "\nTarget censoring value =", censoringvalues,
        "\nTau =", tau,
        "\nNumber of events =", number_events,
        "\nEvent rate =", event_rate,
        "\nActual censoring rate =", actual_censoring_rate,
        "\nTotal number of regression parameters =",
        total_regression_parameters,
        "\nEPV =", EPV,
        "\nEPV_coxph =",EPV_coxph,
        "\n========================================\n"
      )
      
      
      
      
      
      
      
      
      
      
      
      all_y <- data.frame(illt=time[1:n2], ills=event[1:n2])
      names(all_y) <- c("time", "status")
      
      
      y <- data.frame(illt=time[1:use_n], ills=event[1:use_n])
      names(y) <- c("time", "status")
      
      set.seed(  100000L * runs_i +  100L * match(use_n, use_n_all) +tau_j  )
      
      myfold <- sample(  rep(seq_len(5), length.out = use_n)  )
      
      
      
      use_method<-c("grLasso","grSCAD","grMCP")
      result<-matrix(NA_real_, nrow = 4, ncol =  4+4+4+3+mfunctions)
      result_sta<-matrix(0, nrow = 3*4 , ncol =  4+4+4+2)
      for(method_use_i in 1:length(use_method)){
        fit <- cv.grpsurv( X=designmatrix , y=y,group=p_group ,group.multiplier=m,  penalty=use_method[method_use_i] , fold = myfold  )
        
        re_fit <- grpsurv(lambda=fit$lambda.min, X=designmatrix , y=y,group=p_group ,group.multiplier=m,  penalty=use_method[method_use_i] , fold = myfold  )
        
        fit_coef<- fit$fit$beta[,fit$min]
        fit_coef<- re_fit$beta
        
        fit_coef_scalar<-fit_coef[1:q]
        tp_scalar<-sum(fit_coef_scalar[which(gamma!=0)]!=0)
        fp_scalar<-sum(fit_coef_scalar!=0)-tp_scalar
        tn_scalar<-sum(fit_coef_scalar[which(gamma==0)]==0)
        fn_scalar<-sum(fit_coef_scalar==0)-tn_scalar
        
        fit_coef_dev<-fit_coef[(q+1):(q+mfunctions)]
        tp_dev<-sum(fit_coef_dev[which(eta_dev!=0)]!=0)
        fp_dev<-sum(fit_coef_dev!=0)-tp_dev
        tn_dev<-sum(fit_coef_dev[which(eta_dev==0)]==0)
        fn_dev<-sum(fit_coef_dev==0)-tn_dev
        
        
        fit_coef_fun<-fit_coef[-c(1:(q+mfunctions))]
        fit_coef_fun_flag<-rep(0,mfunctions)
        for(subfit_i in 1:mfunctions){
          if(sum(abs(fit_coef_fun[which(as.numeric(p_group[-c(1:(q+mfunctions))]-q-mfunctions)==subfit_i)]))>0) fit_coef_fun_flag[subfit_i]<-1
        }
        tp_fun<-sum(fit_coef_fun_flag[which(true_beta_flag!=0)]!=0)
        fp_fun<-sum(fit_coef_fun_flag!=0)-tp_fun
        tn_fun<-sum(fit_coef_fun_flag[which(true_beta_flag==0)]==0)
        fn_fun<-sum(fit_coef_fun_flag==0)-tn_fun
        
        accuracy_RMSE_ISE <-calculate_estimation_accuracy(estimated_coef = fit_coef,
                                                          true_scalar_coef = RMSE_true_gamma,
                                                          true_eta_dev_coef = RMSE_true_eta_dev,
                                                          true_functional_coef =ISE_true,
                                                          fun_Eigfun =fun_Eigfun ,
                                                          sn_chose = sn_chose,
                                                          q = q,
                                                          mfunctions = mfunctions)
        accuracy_RMSE_ISE$scalar_RMSE
        accuracy_RMSE_ISE$dev_RMSE
        accuracy_RMSE_ISE$functional_ISE
        accuracy_RMSE_ISE$total_functional_ISE
        
        
        result[method_use_i,]<-c(tp_scalar,fp_scalar,tn_scalar,fn_scalar,
                                 tp_dev,fp_dev,tn_dev,fn_dev,
                                 tp_fun,fp_fun,tn_fun,fn_fun,
                                 accuracy_RMSE_ISE$scalar_RMSE,
                                 accuracy_RMSE_ISE$dev_RMSE,
                                 accuracy_RMSE_ISE$functional_ISE,
                                 accuracy_RMSE_ISE$total_functional_ISE)
        
        
        control_ev<-1
        thresholdvalues_lambda<-sqrt(0.8*(mfunctions+mfunctions+q) *control_ev)
        set.seed(9)
        
        
        
        fit <- cv.grpsurv(designmatrix,y,p_group,group.multiplier=m,penalty=use_method[method_use_i], fold = myfold)
        
        
        lam<-matrix(0,dim(fit$fit$beta)[2], q+mfunctions+mfunctions)
        for (lambdai in 1:(dim(fit$fit$beta)[2])) {
          fit_coef_fun<-fit$fit$beta[-c(1:(q+mfunctions)),lambdai]
          fit_coef_fun_flag<-rep(0,mfunctions)
          for(subfit_i in 1:mfunctions){
            if(sum(abs(fit_coef_fun[which(as.numeric(p_group)==(subfit_i+(q+mfunctions)))-(q+mfunctions)]))>0) fit_coef_fun_flag[subfit_i]<-1
          }
          lam[lambdai,]<- c(as.numeric(fit$fit$beta[c(1:(q+mfunctions)),lambdai]!=0),fit_coef_fun_flag)
        }
        
        sel_num<-apply(lam!=0, 1, sum)
        
        start_id<-1
        if(sel_num[1]==0) {
          start_id<-min(which(sel_num==0))
        }
        
        for (ss in (start_id+1):length(sel_num)) {
          if (sum(apply(lam[start_id:ss,]!=0, 2, sum)!=0) >=thresholdvalues_lambda)
            break
        }
        
        ss<-ss-1
        
        
        
        newlam<-fit$fit$lambda[start_id:ss]
        numLambda = length(newlam)
        
        numSample = 100
        nz = array(NA,c(numLambda,length(p_group),numSample))
        max_attempts <- 50L
        
        for (ii in 1:numSample) {
          
          cat(ii, fill = TRUE)
          
          attempt <- 0L
          fit_success <- FALSE
          usecoe <- NULL
          
          while (attempt < max_attempts && !fit_success) {
            
            attempt <- attempt + 1L
            
            sampleInds <- sample( 1:use_n,size = floor(use_n / 2), replace = FALSE    )
            
            aaa <- tryCatch(
              grpsurv(  designmatrix[sampleInds, , drop = FALSE],
                        y[sampleInds, , drop = FALSE], p_group,
                        group.multiplier = m,
                        penalty = use_method[method_use_i],
                        lambda = newlam),
              error = function(e) {
                NULL
              }
            )
            
            
            if (
              !is.null(aaa) &&
              !is.null(aaa$beta) &&
              is.matrix(aaa$beta) &&
              ncol(aaa$beta) == numLambda
            ) {
              usecoe <- aaa$beta
              fit_success <- TRUE
            }
          }
          
          
          if (!fit_success) {
            stop(
              paste0(
                "Stability selection failed for ii = ", ii,
                " after ", max_attempts, " attempts. ",
                "Unable to obtain all ", numLambda,
                " lambda solutions."
              )
            )
          }
          
          nz[, , ii] <- t(1 * (usecoe != 0))
        }
        pr_glasso=apply(nz,c(1,2),sum)/numSample
        
        
        maxPr = apply(pr_glasso,2,max);
        STAL_all<-c(0.9,0.8,0.7,0.6)
        
        for(STAL_i in 1:length(STAL_all)){
          STAL<-STAL_all[STAL_i]
          selVar = maxPr >= STAL;
          
          
          fit_coef_scalar<-as.numeric(selVar)[1:q]
          tp_scalar<-sum(fit_coef_scalar[which(gamma!=0)]!=0)
          fp_scalar<-sum(fit_coef_scalar!=0)-tp_scalar
          tn_scalar<-sum(fit_coef_scalar[which(gamma==0)]==0)
          fn_scalar<-sum(fit_coef_scalar==0)-tn_scalar
          
          fit_coef_dev<-as.numeric(selVar)[(q+1):(q+mfunctions)]
          tp_dev<-sum(fit_coef_dev[which(eta_dev!=0)]!=0)
          fp_dev<-sum(fit_coef_dev!=0)-tp_dev
          tn_dev<-sum(fit_coef_dev[which(eta_dev==0)]==0)
          fn_dev<-sum(fit_coef_dev==0)-tn_dev
          
          fit_coef_fun<-as.numeric(selVar)[-c(1:(q+mfunctions))]
          fit_coef_fun_flag<-rep(0,mfunctions)
          for(subfit_i in 1:mfunctions){
            if(sum(abs(fit_coef_fun[which(as.numeric(p_group[-c(1:(q+mfunctions))]-q-mfunctions)==subfit_i)]))>0) fit_coef_fun_flag[subfit_i]<-1
            
          }
          tp_fun<-sum(fit_coef_fun_flag[which(true_beta_flag!=0)]!=0)
          fp_fun<-sum(fit_coef_fun_flag!=0)-tp_fun
          tn_fun<-sum(fit_coef_fun_flag[which(true_beta_flag==0)]==0)
          fn_fun<-sum(fit_coef_fun_flag==0)-tn_fun
          
          
          Stability_grLasso<-c(tp_scalar,fp_scalar,tn_scalar,fn_scalar,
                               tp_dev,fp_dev,tn_dev,fn_dev,
                               tp_fun,fp_fun,tn_fun,fn_fun ,STAL,control_ev)
          
          result_sta[(method_use_i-1)*4+STAL_i,]<-Stability_grLasso
          
        }
        
        
      }
      
      
      
      result[4,1]<-1- sum(event[1:use_n]) / use_n
      
      result_all <- cbind(result,STAL = NA_real_,control_ev = NA_real_)
      
      result_sta_new <- cbind(
        result_sta[, 1:12, drop = FALSE],scalar_RMSE = rep( NA_real_, nrow(result_sta)),
        dev_RMSE = rep( NA_real_, nrow(result_sta)),
        functional_ISE = matrix( NA_real_, nrow(result_sta),mfunctions),
        total_functional_ISE = rep( NA_real_, nrow(result_sta)),
        result_sta[, 13:14, drop = FALSE])
      
      
      result_all<-rbind(result_all,result_sta_new)
      
      rownames(result_all)<-c( "result_grouplasso","result_groupSCAD","result_groupMCP","censoringvalues",rep("s_glasso",4),rep("s_glSCAD",4),rep("s_gMCP",4))
      colnames(result_all)<-c( "tp_scalar","fp_scalar","tn_scalar","fn_scalar","tp_dev","fp_dev","tn_dev","fn_dev","tp_fun","fp_fun","tn_fun","fn_fun",
                               "scalar_RMSE","dev_RMSE",paste0("functional_ISE_", seq_len(mfunctions)),"total_functional_ISE","STAL","control_ev" )
      result_all
      
      
      
      
      
      
      
      result_coxph <- setNames(  rep(  NA_real_, ncol(result_all) ),colnames(result_all))
      
      
      result_all <- rbind( result_all,   result_coxph = result_coxph      )
      result_all
      
      filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",use_n,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
      write.csv(result_all,filename)
      
      
      
      
      diagnostic_information <- data.frame(
        simulation_run = runs_i,
        rho = rho,
        sample_size = use_n,
        target_censoring = censoringvalues,
        actual_censoring = actual_censoring_rate,
        
        average_retained_fpc = average_retained_fpc,
        total_retained_fpc = total_retained_fpc,
        
        number_scalar_parameters = number_scalar_parameters,
        number_derivative_parameters = number_derivative_parameters,
        number_functional_parameters = number_functional_parameters,
        total_regression_parameters = total_regression_parameters,
        
        number_events = number_events,
        event_rate = event_rate,
        EPV = EPV,
        EPV_coxph = number_events / (q + total_retained_fpc)
      )
      
      filename_EPV=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",use_n,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,"EPV.csv",sep="")
      write.csv(diagnostic_information,filename_EPV)
      
      
    }
  }
  
}
