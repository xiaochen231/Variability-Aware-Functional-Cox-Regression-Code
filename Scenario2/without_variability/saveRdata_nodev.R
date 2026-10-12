nruns <- 100


methods <- 6


methods_dev <- 3

methods_MSE <- 6


method_names <- c(
  "lasso", "SCAD", "MCP",
  "lasso-w", "SCAD-w", "MCP-w"
)

censoringvalues_all<-c(0.2)
n1_all<-c(200,300,500,1000,2000)
q_all<-c(10)
q_true_all <- c(2)
rho_all<-c(0.3,0.5,0.8)


for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (q_i in 1:length(q_all)) {
      for (q_true_i in 1:length(q_true_all)) {
        for (rho_k in 1:length(rho_all)) {
          
          n1 <- n1_all[n_j]
          mfunctions <- q_all[q_i]
          mfunctions_true <- q_true_all[q_true_i]
          q <- q_all[q_i]
          q_true <- q_true_all[q_true_i]
          censoringvalues <- censoringvalues_all[cen_i]
          rho <- rho_all[rho_k]
          
          
          tp_scalar <- matrix(0, nruns, methods)
          fp_scalar <- matrix(0, nruns, methods)
          tn_scalar <- matrix(0, nruns, methods)
          fn_scalar <- matrix(0, nruns, methods)
          
          
          tp_dev <- matrix(0, nruns, methods_dev)
          fp_dev <- matrix(0, nruns, methods_dev)
          tn_dev <- matrix(0, nruns, methods_dev)
          fn_dev <- matrix(0, nruns, methods_dev)
          
          
          tp_fun <- matrix(0, nruns, methods)
          fp_fun <- matrix(0, nruns, methods)
          tn_fun <- matrix(0, nruns, methods)
          fn_fun <- matrix(0, nruns, methods)
          
          censoring <- matrix(0, nruns, 1)
          
          scalar_RMSE_all <- matrix(0, nruns, methods_MSE)
          
          
          dev_RMSE_all <- matrix(0, nruns, methods_dev)
          
          total_ISE_all <- matrix(0, nruns, methods_MSE)
          
          
          for (runs_i in 1:nruns) {
            
            filename <- paste(
              runs_i,
              "_rho", rho,
              "_censoringvalues", censoringvalues,
              "_n", n1,
              "_mfunctions", mfunctions,
              "_mfunctions_true", mfunctions_true,
              "_q", q,
              "_q_true", q_true,
              ".csv",
              sep = ""
            )
            
            result <- read.csv(
              filename,
              check.names = FALSE,
              stringsAsFactors = FALSE
            )
            
            
            censoring[runs_i] <- result[4, "tp_scalar"]
            
            
            
            scalar_RMSE_all[runs_i, ] <-
              as.numeric(result[c(1:3, 5:7), "scalar_RMSE"])
            
            total_ISE_all[runs_i, ] <-
              as.numeric(result[c(1:3, 5:7), "total_functional_ISE"])
            
            
            dev_RMSE_all[runs_i, ] <-
              as.numeric(result[1:3, "dev_RMSE"])
            
            
            
            result_select <- result[c(1:3, 5:7), ]
            
            tp_scalar[runs_i, ] <-
              as.numeric(result_select[, "tp_scalar"])
            
            fp_scalar[runs_i, ] <-
              as.numeric(result_select[, "fp_scalar"])
            
            tn_scalar[runs_i, ] <-
              as.numeric(result_select[, "tn_scalar"])
            
            fn_scalar[runs_i, ] <-
              as.numeric(result_select[, "fn_scalar"])
            
            
            tp_fun[runs_i, ] <-
              as.numeric(result_select[, "tp_fun"])
            
            fp_fun[runs_i, ] <-
              as.numeric(result_select[, "fp_fun"])
            
            tn_fun[runs_i, ] <-
              as.numeric(result_select[, "tn_fun"])
            
            fn_fun[runs_i, ] <-
              as.numeric(result_select[, "fn_fun"])
            
            
            
            tp_dev[runs_i, ] <-
              as.numeric(result[1:3, "tp_dev"])
            
            fp_dev[runs_i, ] <-
              as.numeric(result[1:3, "fp_dev"])
            
            tn_dev[runs_i, ] <-
              as.numeric(result[1:3, "tn_dev"])
            
            fn_dev[runs_i, ] <-
              as.numeric(result[1:3, "fn_dev"])
          }
          
          
          
          colnames(tp_scalar) <- method_names
          colnames(fp_scalar) <- method_names
          colnames(tn_scalar) <- method_names
          colnames(fn_scalar) <- method_names
          
          colnames(tp_fun) <- method_names
          colnames(fp_fun) <- method_names
          colnames(tn_fun) <- method_names
          colnames(fn_fun) <- method_names
          
          colnames(tp_dev) <- method_names[1:3]
          colnames(fp_dev) <- method_names[1:3]
          colnames(tn_dev) <- method_names[1:3]
          colnames(fn_dev) <- method_names[1:3]
          
          colnames(scalar_RMSE_all) <- method_names
          colnames(dev_RMSE_all) <- method_names[1:3]
          colnames(total_ISE_all) <- method_names
          
          
          tpr_scalar <- tp_scalar / (tp_scalar + fn_scalar)
          fpr_scalar <- fp_scalar / (fp_scalar + tn_scalar)
          
          tpr_fun <- tp_fun / (tp_fun + fn_fun)
          fpr_fun <- fp_fun / (fp_fun + tn_fun)
          
          
          tpr_dev <- tp_dev / (tp_dev + fn_dev)
          fpr_dev <- fp_dev / (fp_dev + tn_dev)
          
          
          fdr_scalar <- fp_scalar / (tp_scalar + fp_scalar)
          fnr_scalar <- fn_scalar / (tp_scalar + fn_scalar)
          sensitivity_scalar <- 1 - fnr_scalar
          specificity_scalar <- 1 - fp_scalar / (fp_scalar + tn_scalar)
          
          
          
          fdr_dev <- fp_dev / (tp_dev + fp_dev)
          fnr_dev <- fn_dev / (tp_dev + fn_dev)
          sensitivity_dev <- 1 - fnr_dev
          specificity_dev <- 1 - fp_dev / (fp_dev + tn_dev)
          
          
          fdr_fun <- fp_fun / (tp_fun + fp_fun)
          fnr_fun <- fn_fun / (tp_fun + fn_fun)
          sensitivity_fun <- 1 - fnr_fun
          specificity_fun <- 1 - fp_fun / (fp_fun + tn_fun)
          
          
          
          
          
          
          tp_all <- cbind(
            tp_scalar[, 1:3, drop = FALSE] +
              tp_fun[, 1:3, drop = FALSE] +
              tp_dev,
            
            tp_scalar[, 4:6, drop = FALSE] +
              tp_fun[, 4:6, drop = FALSE]
          )
          
          fp_all <- cbind(
            fp_scalar[, 1:3, drop = FALSE] +
              fp_fun[, 1:3, drop = FALSE] +
              fp_dev,
            
            fp_scalar[, 4:6, drop = FALSE] +
              fp_fun[, 4:6, drop = FALSE]
          )
          
          fn_all <- cbind(
            fn_scalar[, 1:3, drop = FALSE] +
              fn_fun[, 1:3, drop = FALSE] +
              fn_dev,
            
            fn_scalar[, 4:6, drop = FALSE] +
              fn_fun[, 4:6, drop = FALSE]
          )
          
          tn_all <- cbind(
            tn_scalar[, 1:3, drop = FALSE] +
              tn_fun[, 1:3, drop = FALSE] +
              tn_dev,
            
            tn_scalar[, 4:6, drop = FALSE] +
              tn_fun[, 4:6, drop = FALSE]
          )
          
          colnames(tp_all) <- method_names
          colnames(fp_all) <- method_names
          colnames(fn_all) <- method_names
          colnames(tn_all) <- method_names
          
          
          
          fdr_all <- fp_all / (tp_all + fp_all)
          fnr_all <- fn_all / (tp_all + fn_all)
          
          sensitivity_all <- 1 - fnr_all
          specificity_all <- 1 - fp_all / (fp_all + tn_all)
          
          
          apply(fdr_all, 2, mean, na.rm = TRUE)
          apply(fnr_all, 2, mean, na.rm = TRUE)
          apply(sensitivity_all, 2, mean, na.rm = TRUE)
          apply(specificity_all, 2, mean, na.rm = TRUE)
          
          apply(fp_all, 2, mean, na.rm = TRUE)
          apply(tp_all, 2, mean, na.rm = TRUE)
          
          apply(fn_all, 2, mean, na.rm = TRUE)
          apply(tn_all, 2, mean, na.rm = TRUE)
          
          apply(tpr_scalar, 2, mean, na.rm = TRUE)
          apply(fpr_scalar, 2, mean, na.rm = TRUE)
          apply(tpr_fun, 2, mean, na.rm = TRUE)
          apply(fpr_fun, 2, mean, na.rm = TRUE)
          
          apply(tpr_dev, 2, mean, na.rm = TRUE)
          apply(fpr_dev, 2, mean, na.rm = TRUE)
          
          summary(censoring)
          plot(censoring)
          
          
          
          save_file <- paste(
            "q", q,
            "true", q_true,
            "n", n1,
            "c", censoringvalues,
            "rho", rho,
            "nodev.RData",
            sep = ""
          )
          
          save(
            list = c(
              "censoring",
              "tp_all",
              "fp_all",
              "scalar_RMSE_all",
              "dev_RMSE_all",
              "total_ISE_all",
              "method_names"
            ),
            file = save_file
          )
        }
      }
    }
  }
}
