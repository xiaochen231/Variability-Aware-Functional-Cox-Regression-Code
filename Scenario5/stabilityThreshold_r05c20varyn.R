



library(ggplot2)
library(dplyr)
library(gridExtra)






data_dir <- "."






runs_all <- 1:100

rho <- 0.5
censoringvalues <- 0.2

n_all <- c(200, 300, 500, 1000, 2000)

mfunctions <- 300
mfunctions_true <- 2

q <- 300
q_true <- 2


pi_all <- c(0.9, 0.8, 0.7, 0.6)





plot_data <- data.frame(
  run       = integer(),
  n         = integer(),
  panel     = character(),
  pi        = numeric(),
  method    = character(),
  TP        = numeric(),
  FP        = numeric(),
  stringsAsFactors = FALSE
)






to_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}



calculate_tp_fp <- function(one_row) {
  
  tp_all <- to_numeric(one_row$tp_scalar) +
    to_numeric(one_row$tp_dev) +
    to_numeric(one_row$tp_fun)
  
  fp_all <- to_numeric(one_row$fp_scalar) +
    to_numeric(one_row$fp_dev) +
    to_numeric(one_row$fp_fun)
  
  return(
    data.frame(
      TP = tp_all,
      FP = fp_all
    )
  )
}



add_one_result <- function(current_data,
                           one_row,
                           run_value,
                           n_value,
                           panel_value,
                           pi_value,
                           method_value) {
  
  values <- calculate_tp_fp(one_row)
  
  new_result <- data.frame(
    run    = run_value,
    n      = n_value,
    panel  = panel_value,
    pi     = pi_value,
    method = method_value,
    TP     = values$TP,
    FP     = values$FP,
    stringsAsFactors = FALSE
  )
  
  bind_rows(current_data, new_result)
}





for (n1 in n_all) {
  
  for (runs_i in runs_all) {
    
    
    filename <- paste0(
      runs_i,
      "_rho", rho,
      "_censoringvalues", censoringvalues,
      "_n", n1,
      "_mfunctions", mfunctions,
      "_mfunctions_true", mfunctions_true,
      "_q", q,
      "_q_true", q_true,
      ".csv"
    )
    
    file_path <- file.path(data_dir, filename)
    
    
    
    
    
    if (!file.exists(file_path)) {
      warning(paste("文件不存在：", file_path))
      next
    }
    
    
    
    
    
    one_file <- tryCatch(
      {
        read.csv(
          file_path,
          header = TRUE,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      },
      error = function(e) {
        warning(
          paste(
            "无法读取文件：",
            file_path,
            "；错误信息：",
            e$message
          )
        )
        return(NULL)
      }
    )
    
    if (is.null(one_file)) {
      next
    }
    
    
    
    
    
    names(one_file)[1] <- "result_type"
    
    one_file$result_type <- trimws(
      as.character(one_file$result_type)
    )
    
    
    
    
    
    required_columns <- c(
      "result_type",
      "tp_scalar",
      "fp_scalar",
      "tp_dev",
      "fp_dev",
      "tp_fun",
      "fp_fun",
      "STAL"
    )
    
    missing_columns <- setdiff(
      required_columns,
      names(one_file)
    )
    
    if (length(missing_columns) > 0) {
      
      warning(
        paste0(
          "文件 ",
          filename,
          " 缺少变量：",
          paste(missing_columns, collapse = ", ")
        )
      )
      
      next
    }
    
    
    
    
    
    original_method_names <- c(
      "result_grouplasso",
      "result_groupSCAD",
      "result_groupMCP"
    )
    
    original_method_labels <- c(
      "lasso",
      "SCAD",
      "MCP"
    )
    
    for (method_i in seq_along(original_method_names)) {
      
      selected_row <- one_file[
        one_file$result_type == original_method_names[method_i],
        ,
        drop = FALSE
      ]
      
      if (nrow(selected_row) == 0) {
        
        warning(
          paste0(
            filename,
            " 中没有找到 ",
            original_method_names[method_i]
          )
        )
        
        next
      }
      
      
      selected_row <- selected_row[1, , drop = FALSE]
      
      plot_data <- add_one_result(
        current_data = plot_data,
        one_row      = selected_row,
        run_value    = runs_i,
        n_value      = n1,
        panel_value  = "Original",
        pi_value     = NA_real_,
        method_value = original_method_labels[method_i]
      )
    }
    
    
    
    
    
    s_method_names <- c(
      "s_glasso",
      "s_glSCAD",
      "s_gMCP"
    )
    
    s_method_labels <- c(
      "S-lasso",
      "S-SCAD",
      "S-MCP"
    )
    
    
    one_file$STAL_numeric <- to_numeric(one_file$STAL)
    
    for (pi_i in pi_all) {
      
      for (method_i in seq_along(s_method_names)) {
        
        selected_row <- one_file[
          one_file$result_type == s_method_names[method_i] &
            !is.na(one_file$STAL_numeric) &
            abs(one_file$STAL_numeric - pi_i) < 1e-8,
          ,
          drop = FALSE
        ]
        
        if (nrow(selected_row) == 0) {
          
          warning(
            paste0(
              filename,
              " 中没有找到方法 ",
              s_method_names[method_i],
              " 且 STAL = ",
              pi_i
            )
          )
          
          next
        }
        
        
        selected_row <- selected_row[1, , drop = FALSE]
        
        panel_name <- paste0("\u03c0 = ", pi_i)
        
        plot_data <- add_one_result(
          current_data = plot_data,
          one_row      = selected_row,
          run_value    = runs_i,
          n_value      = n1,
          panel_value  = panel_name,
          pi_value     = pi_i,
          method_value = s_method_labels[method_i]
        )
      }
    }
    
    cat("已经读取：", filename, "\n")
  }
}





cat("\n总共读取到", nrow(plot_data), "条方法结果。\n\n")

print(head(plot_data))


count_table <- plot_data %>%
  count(n, panel, method, name = "number_of_runs")

print(count_table)





library(dplyr)
library(ggplot2)





plot_data <- plot_data %>%
  mutate(
    
    
    grid_y = paste0("n = ", n),
    
    
    
    grid_x = case_when(
      panel == "Original"       ~ "phantom('Original')",
      abs(pi - 0.9) < 1e-8      ~ "pi[thr] == 0.9",
      abs(pi - 0.8) < 1e-8      ~ "pi[thr] == 0.8",
      abs(pi - 0.7) < 1e-8      ~ "pi[thr] == 0.7",
      abs(pi - 0.6) < 1e-8      ~ "pi[thr] == 0.6",
      TRUE                      ~ NA_character_
    )
  )







plot_data$grid_y <- factor(
  plot_data$grid_y,
  levels = c(
    "n = 200",
    "n = 300",
    "n = 500",
    "n = 1000",
    "n = 2000"
  )
)


plot_data$grid_x <- factor(
  plot_data$grid_x,
  levels = c(
    "phantom('Original')",
    "pi[thr] == 0.9",
    "pi[thr] == 0.8",
    "pi[thr] == 0.7",
    "pi[thr] == 0.6"
  )
)






plot_data$method <- factor(
  plot_data$method,
  levels = c(
    "lasso",
    "SCAD",
    "MCP",
    "S-lasso",
    "S-SCAD",
    "S-MCP"
  )
)






method_colours <- c(
  
  
  "lasso"   = "#E41A1C",
  "SCAD"    = "#4DAF4A",
  "MCP"     = "#377EB8",
  
  
  "S-lasso" = "#A50F15",
  "S-SCAD"  = "#238B45",
  "S-MCP"   = "#08519C"
)









method_labels <- c(
  "lasso"   = "lasso",
  "SCAD"    = "SCAD",
  "MCP"     = "MCP",
  "S-lasso" = "S-lasso",
  "S-SCAD"  = "S-SCAD",
  "S-MCP"   = "S-MCP"
)


p_TP <- ggplot(
  data = plot_data,
  aes(x = method, y = TP)
) +
  
  geom_boxplot(
    aes(colour = method),
    notch = FALSE,
    width = 0.65,
    linewidth = 0.7
  ) +
  
  xlab("") +
  ylab("TP") +
  
  facet_grid(
    grid_y ~ grid_x,
    scales = "free_x",
    labeller = labeller(
      grid_x = label_parsed
    )
  ) +
  
  scale_colour_manual(
    values = method_colours
  ) +
  
  scale_x_discrete(
    labels = method_labels
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      lineheight = 0.8,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 22
    ),
    
    axis.text.x = element_text(
      size = 14
    ),
    
    axis.text.y = element_text(
      size = 18
    ),
    
    legend.text = element_text(
      size = 22
    ),
    
    strip.text.x = element_text(
      size = 20,
      colour = "purple"
    ),
    
    strip.text.y = element_text(
      size = 15,
      colour = "purple"
    ),
    
    legend.position = "none"
  ) +
  
  guides(fill = "none")

print(p_TP)


p_FP <- ggplot(
  data = plot_data,
  aes(x = method, y = FP)
) +
  
  geom_boxplot(
    aes(colour = method),
    notch = FALSE,
    width = 0.65,
    linewidth = 0.7
  ) +
  
  xlab("") +
  ylab("FP") +
  
  facet_grid(
    grid_y ~ grid_x,
    scales = "free_x",
    labeller = labeller(
      grid_x = label_parsed
    )
  ) +
  
  scale_colour_manual(
    values = method_colours
  ) +
  
  scale_x_discrete(
    labels = method_labels
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      lineheight = 0.8,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 22
    ),
    
    axis.text.x = element_text(
      size = 14
    ),
    
    axis.text.y = element_text(
      size = 18
    ),
    
    legend.text = element_text(
      size = 22
    ),
    
    strip.text.x = element_text(
      size = 20,
      colour = "purple"
    ),
    
    strip.text.y = element_text(
      size = 15,
      colour = "purple"
    ),
    
    legend.position = "none"
  ) +
  
  guides(fill = "none")

print(p_FP)





combined_plot <- grid.arrange(
  p_TP,
  p_FP,
  ncol = 1,
  heights = c(1, 1)
)





combined_plot1 <- grid.arrange(p_TP, p_FP, ncol = 1, nrow = 2)
ggsave(filename = "q300stabilityThreshold_r05c20varyn.pdf", plot = combined_plot1, width = 15, height = 15, units = "in", device = "pdf")










