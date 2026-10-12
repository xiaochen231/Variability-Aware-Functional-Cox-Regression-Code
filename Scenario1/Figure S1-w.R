library(ggplot2)
library(plyr)
library(gridExtra)

censoringvalues_plot_all<-c( 0.1,0.2,0.3,0.5)
n1_plot_all<-c(200,300,500,1000,2000)
q_plot_all<-c(10)
q_true_plot_all<-c(2)
rho_plot_all<-c(0.5)
MAE_FP<-c()
MAE_TP<-c()
scalar_MSE<-c()
dev_MSE<-c()
fun_ISE<-c()


for (cen_plot_i in 1:length(censoringvalues_plot_all)) {
  for (n_plot_j in 1:length(n1_plot_all)) {
    for (q_plot_i in 1:length(q_plot_all)) {
      for (q_true_plot_i in 1:length(q_true_plot_all)) {
        for (rho_plot_k in 1:length(rho_plot_all)) {
          censoringvalues<-censoringvalues_plot_all[cen_plot_i]
          n1<-n1_plot_all[n_plot_j]
          q<-q_plot_all[q_plot_i]
          q_true<-q_true_plot_all[q_true_plot_i]
          rho<-rho_plot_all[rho_plot_k]
          save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,".RData",sep="")
          load(save_file)
          MAE_TP <- c(MAE_TP,c(tp_all[,c(1,2,3)]))
          MAE_FP <- c(MAE_FP,c(fp_all[,c(1,2,3)]))
          
          scalar_MSE<-c(scalar_MSE,c(scalar_RMSE_all[,c(1,2,3,4)]) )
          dev_MSE<-c(dev_MSE,c(dev_RMSE_all[,c(1,2,3)]))
          fun_ISE<-c(fun_ISE,c(total_ISE_all[,c(1,2,3,4)]))
          
          
          save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,"nodev.RData",sep="")
          load(save_file)
          
          scalar_MSE<-c(scalar_MSE,c(scalar_RMSE_all[,c(4,5,6)]) )
          fun_ISE<-c(fun_ISE,c(total_ISE_all[,c(4,5,6)]))
          
          print(save_file)
          print(mean(censoring))
        }
        
      }
    }
  }
}

method <-   rep(c(rep("lasso", 100), rep("SCAD", 100), rep("MCP", 100)),length(censoringvalues_plot_all)*length(n1_plot_all)*length(q_plot_all)*length(q_true_plot_all)*length(rho_plot_all))

method <- factor(method, levels = c("lasso", "SCAD", "MCP"))


grid_x <- rep(c(rep("n = 200", 300*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
                rep("n = 300", 300*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
                rep("n = 500", 300*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
                rep("n = 1000", 300*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
                rep("n = 2000", 300*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all))),
              length(censoringvalues_plot_all))
grid_x <- factor(grid_x, levels = c("n = 200", "n = 300", "n = 500", "n = 1000" , "n = 2000"))



grid_y <- c(rep("Censoring = 10%", 300*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
            rep("Censoring = 20%", 300*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
            rep("Censoring = 30%", 300*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
            rep("Censoring = 50%", 300*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)))
grid_y <- factor(grid_y, levels = c(  "Censoring = 10%", "Censoring = 20%",  "Censoring = 30%",  "Censoring = 50%"))


data_TP <- data.frame(MAE = MAE_TP, method = method, grid_x = grid_x, grid_y = grid_y)
data_FP <- data.frame(MAE = MAE_FP, method = method, grid_x = grid_x, grid_y = grid_y)



p_TP <- ggplot(data = data_TP, aes(x = method, y = MAE))  +
  geom_boxplot(aes(colour = method), notch = F) +
  xlab("") + ylab("TP")+
  theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold")) +
  facet_grid(grid_x~grid_y) +
  theme(axis.title.y = element_text(size = 22), axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 22)) +
  theme(strip.text.x = element_text(size = 20, colour = "purple")) +
  guides(fill = FALSE) +
  theme(strip.text.y = element_text(size = 15, colour = "purple")) +
  theme(legend.position = "none")+
  scale_colour_manual(values = c(
    "lasso" = "#E41A1C",
    "SCAD"  = "#4DAF4A",
    "MCP"   = "#377EB8",
    "FLCRM"   = "#984EA3"  ,
    "lasso-w" = "#FC9272",
    "SCAD-w"  = "#A1D99B",
    "MCP-w"   = "#9ECAE1"
  ))

p_FP <- ggplot(data = data_FP, aes(x = method, y = MAE))  +
  geom_boxplot(aes(colour = method), notch = F) +
  xlab("") + ylab("FP")+
  theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold")) +
  facet_grid(grid_x~grid_y) +
  theme(axis.title.y = element_text(size = 22), axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 22)) +
  theme(strip.text.x = element_text(size = 20, colour = "purple")) +
  guides(fill = FALSE) +
  theme(strip.text.y = element_text(size = 15, colour = "purple")) +
  theme(legend.position = "none")+
  scale_colour_manual(values = c(
    "lasso" = "#E41A1C",
    "SCAD"  = "#4DAF4A",
    "MCP"   = "#377EB8",
    "FLCRM"   = "#984EA3",
    "lasso-w" = "#FC9272",
    "SCAD-w"  = "#A1D99B",
    "MCP-w"   = "#9ECAE1"
  ))









grid_x2 <- rep(c(rep("n = 200", 700*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
                 rep("n = 300", 700*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
                 rep("n = 500", 700*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
                 rep("n = 1000", 700*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
                 rep("n = 2000", 700*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all))),
               length(censoringvalues_plot_all))
grid_x2 <- factor(grid_x2, levels = c("n = 200", "n = 300", "n = 500", "n = 1000" , "n = 2000"))



grid_y2 <- c(rep("Censoring = 10%", 700*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
             rep("Censoring = 20%", 700*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
             rep("Censoring = 30%", 700*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),
             rep("Censoring = 50%", 700*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)))
grid_y2 <- factor(grid_y2, levels = c(  "Censoring = 10%", "Censoring = 20%",  "Censoring = 30%",  "Censoring = 50%"))


method2 <-   rep(c(rep("lasso", 100), rep("SCAD", 100), rep("MCP", 100), rep("FLCRM", 100),rep("lasso-w", 100), rep("SCAD-w", 100), rep("MCP-w", 100)),length(censoringvalues_plot_all)*length(n1_plot_all)*length(q_plot_all)*length(q_true_plot_all)*length(rho_plot_all))
method2 <- factor(method2, levels = c("lasso", "lasso-w","SCAD","SCAD-w", "MCP",  "MCP-w","FLCRM"))


data_scalar_MSE <- data.frame(MAE = scalar_MSE, method = method2, grid_x = grid_x2, grid_y = grid_y2)

P_scalar_MSE <- ggplot(data = data_scalar_MSE, aes(x = method, y = MAE))  +
  geom_boxplot(aes(colour = method2), notch = F) +
  xlab("") + ylab(expression(RMSE[gamma]))+
  theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold")) +
  facet_grid(grid_x~grid_y) +
  theme(axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18, angle = 30, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 22)) +
  theme(strip.text.x = element_text(size = 20, colour = "purple")) +
  guides(fill = FALSE) +
  theme(strip.text.y = element_text(size = 15, colour = "purple")) +
  theme(legend.position = "none")+
  scale_colour_manual(values = c(
    "lasso" = "#E41A1C",
    "SCAD"  = "#4DAF4A",
    "MCP"   = "#377EB8",
    "FLCRM"   = "#984EA3"  ,
    "lasso-w" = "#FC9272",
    "SCAD-w"  = "#A1D99B",
    "MCP-w"   = "#9ECAE1"
  ))

data_fun_ISE <- data.frame(MAE = fun_ISE, method = method2, grid_x = grid_x2, grid_y = grid_y2)

P_fun_ISE <- ggplot(data = data_fun_ISE, aes(x = method, y = MAE))  +
  geom_boxplot(aes(colour = method2), notch = F) +
  xlab("") + ylab("ISE")+
  theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold")) +
  facet_grid(grid_x~grid_y, scales = "free_y") +
  theme(axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18, angle = 30, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 22)) +
  theme(strip.text.x = element_text(size = 20, colour = "purple")) +
  guides(fill = FALSE) +
  theme(strip.text.y = element_text(size = 15, colour = "purple")) +
  theme(legend.position = "none")+
  scale_colour_manual(values = c(
    "lasso" = "#E41A1C",
    "SCAD"  = "#4DAF4A",
    "MCP"   = "#377EB8",
    "FLCRM"   = "#984EA3"   ,
    "lasso-w" = "#FC9272",
    "SCAD-w"  = "#A1D99B",
    "MCP-w"   = "#9ECAE1"
  ))


data_dev_MSE <- data.frame(MAE = dev_MSE, method = method, grid_x = grid_x , grid_y = grid_y)

P_dev_MSE <- ggplot(data = data_dev_MSE, aes(x = method, y = MAE))  +
  geom_boxplot(aes(colour = method), notch = F) +
  xlab("") +  ylab(expression(RMSE[eta]))+
  theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold")) +
  facet_grid(grid_x~grid_y) +
  theme(axis.title.y = element_text(size = 22), axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 22)) +
  theme(strip.text.x = element_text(size = 20, colour = "purple")) +
  guides(fill = FALSE) +
  theme(strip.text.y = element_text(size = 15, colour = "purple")) +
  theme(legend.position = "none")+
  scale_colour_manual(values = c(
    "lasso" = "#E41A1C",
    "SCAD"  = "#4DAF4A",
    "MCP"   = "#377EB8",
    "FLCRM"   = "#984EA3" ,
    "lasso-w" = "#FC9272",
    "SCAD-w"  = "#A1D99B",
    "MCP-w"   = "#9ECAE1"
  ))

P_scalar_MSE
P_dev_MSE
P_fun_ISE
grid.arrange( P_scalar_MSE,P_dev_MSE,P_fun_ISE ,ncol=1,nrow=3)


combined_plot1 <- grid.arrange(p_TP, p_FP, ncol = 1, nrow = 2)
ggsave(filename = "FigureS1_1.pdf", plot = combined_plot1, width = 15, height = 15, units = "in", device = "pdf")

combined_plot2 <- grid.arrange(P_scalar_MSE, P_dev_MSE, P_fun_ISE, ncol = 1, nrow = 3)
ggsave(filename = "FigureS1_2-w.pdf", plot = combined_plot2, width = 18, height = 22, units = "in", device = "pdf")
