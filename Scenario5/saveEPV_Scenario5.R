library(dplyr)
library(purrr)
library(readr)
library(tidyr)


result_dir <- "."


rho_values <- 0.5
censoring_values <- c(0.1, 0.2, 0.3, 0.5)
sample_size_values <- c(200, 300, 500, 1000, 2000)

runs <- 1:100

q <- 300
q_true <- 2


mfunctions <- 300
mfunctions_true <- 2



file_information <- tidyr::crossing(
  simulation_run = runs,
  rho = rho_values,
  target_censoring = censoring_values,
  sample_size = sample_size_values
) %>%
  mutate(
    filename = paste0(
      simulation_run,
      "_rho", rho,
      "_censoringvalues", target_censoring,
      "_n", sample_size,
      "_mfunctions", mfunctions,
      "_mfunctions_true", mfunctions_true,
      "_q", q,
      "_q_true", q_true,
      "EPV.csv"
    ),
    filepath = file.path(result_dir, filename),
    file_exists = file.exists(filepath)
  )


nrow(file_information)
table(file_information$file_exists)
missing_files <- file_information %>%
  filter(!file_exists)

nrow(missing_files)
safe_read_EPV <- function(filepath,
                          expected_run,
                          expected_rho,
                          expected_censoring,
                          expected_n) {
  
  tryCatch({
    
    dat <- read.csv(
      filepath,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    
    if (ncol(dat) > 0) {
      
      first_name <- names(dat)[1]
      
      if (first_name %in% c("", "X", "...1")) {
        dat <- dat[, -1, drop = FALSE]
      }
    }
    
    
    if (nrow(dat) == 0) {
      warning("文件没有数据：", basename(filepath))
      return(NULL)
    }
    
    
    
    if (nrow(dat) > 1) {
      warning(
        "文件包含多行数据：",
        basename(filepath),
        "；行数 = ", nrow(dat)
      )
    }
    
    
    dat %>%
      mutate(
        source_file = basename(filepath),
        expected_run = expected_run,
        expected_rho = expected_rho,
        expected_censoring = expected_censoring,
        expected_sample_size = expected_n,
        read_success = TRUE
      )
    
  }, error = function(e) {
    
    message(
      "读取失败：", basename(filepath),
      "；错误信息：", conditionMessage(e)
    )
    
    return(NULL)
  })
}


existing_file_information <- file_information %>%
  filter(file_exists)


scenario1_raw <- purrr::pmap_dfr(
  list(
    filepath = existing_file_information$filepath,
    expected_run = existing_file_information$simulation_run,
    expected_rho = existing_file_information$rho,
    expected_censoring =
      existing_file_information$target_censoring,
    expected_n =
      existing_file_information$sample_size
  ),
  safe_read_EPV
)
nrow(scenario1_raw)
head(scenario1_raw)

numeric_variables <- c(
  "simulation_run",
  "rho",
  "sample_size",
  "target_censoring",
  "actual_censoring",
  "average_retained_fpc",
  "total_retained_fpc",
  "number_scalar_parameters",
  "number_derivative_parameters",
  "number_functional_parameters",
  "total_regression_parameters",
  "number_events",
  "event_rate",
  "EPV",
  "EPV_coxph"
)

scenario1_raw <- scenario1_raw %>%
  mutate(
    across(
      any_of(numeric_variables),
      as.numeric
    )
  )
parameter_check <- scenario1_raw %>%
  mutate(
    run_match =
      simulation_run == expected_run,
    
    rho_match =
      abs(rho - expected_rho) < 1e-10,
    
    censoring_match =
      abs(target_censoring - expected_censoring) < 1e-10,
    
    sample_size_match =
      sample_size == expected_sample_size
  )
parameter_mismatch <- parameter_check %>%
  filter(
    !run_match |
      !rho_match |
      !censoring_match |
      !sample_size_match
  )

nrow(parameter_mismatch)

scenario1_raw <- scenario1_raw %>%
  mutate(
    calculated_total_parameters =
      number_scalar_parameters +
      number_derivative_parameters +
      number_functional_parameters,
    
    total_parameter_match =
      total_regression_parameters ==
      calculated_total_parameters,
    
    calculated_EPV =
      number_events / total_regression_parameters,
    
    EPV_difference =
      EPV - calculated_EPV
  )
scenario1_raw %>%
  filter(
    !total_parameter_match |
      abs(EPV_difference) > 1e-10
  ) %>%
  select(
    source_file,
    total_regression_parameters,
    calculated_total_parameters,
    EPV,
    calculated_EPV,
    EPV_difference
  )
scenario1_run_check <- scenario1_raw %>%
  group_by(
    rho,
    sample_size,
    target_censoring
  ) %>%
  summarise(
    number_records = n(),
    number_unique_runs = n_distinct(simulation_run),
    minimum_run = min(simulation_run, na.rm = TRUE),
    maximum_run = max(simulation_run, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(sample_size, target_censoring)

scenario1_run_check



incomplete_settings <- scenario1_run_check %>%
  filter(number_unique_runs != 100)

incomplete_settings
scenario1_summary_numeric <- scenario1_raw %>%
  group_by(
    rho,
    sample_size,
    target_censoring
  ) %>%
  summarise(
    number_simulations = n(),
    
    actual_censoring_mean =
      mean(actual_censoring, na.rm = TRUE),
    actual_censoring_sd =
      sd(actual_censoring, na.rm = TRUE),
    
    average_retained_fpc_mean =
      mean(average_retained_fpc, na.rm = TRUE),
    average_retained_fpc_sd =
      sd(average_retained_fpc, na.rm = TRUE),
    
    total_retained_fpc_mean =
      mean(total_retained_fpc, na.rm = TRUE),
    total_retained_fpc_sd =
      sd(total_retained_fpc, na.rm = TRUE),
    
    total_regression_parameters_mean =
      mean(total_regression_parameters, na.rm = TRUE),
    total_regression_parameters_sd =
      sd(total_regression_parameters, na.rm = TRUE),
    
    number_events_mean =
      mean(number_events, na.rm = TRUE),
    number_events_sd =
      sd(number_events, na.rm = TRUE),
    
    EPV_mean =
      mean(EPV, na.rm = TRUE),
    EPV_sd =
      sd(EPV, na.rm = TRUE),
    
    EPV_min =
      min(EPV, na.rm = TRUE),
    EPV_max =
      max(EPV, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(sample_size, target_censoring)

scenario1_summary_numeric



scenario1_publication_table <- scenario1_summary_numeric %>%
  mutate(
    `Target censoring (%)` =
      sprintf("%.0f", 100 * target_censoring),
    
    `Average retained FPCs` =
      sprintf(
        "%.3f (%.3f)",
        average_retained_fpc_mean,
        average_retained_fpc_sd
      ),
    
    `Total retained FPCs` =
      sprintf(
        "%.3f (%.3f)",
        total_retained_fpc_mean,
        total_retained_fpc_sd
      ),
    
    `Total regression parameters` =
      sprintf(
        "%.3f (%.3f)",
        total_regression_parameters_mean,
        total_regression_parameters_sd
      ),
    
    `Number of events` =
      sprintf(
        "%.3f (%.3f)",
        number_events_mean,
        number_events_sd
      ),
    
    `EPV` =
      sprintf(
        "%.3f (%.3f)",
        EPV_mean,
        EPV_sd
      )
  ) %>%
  select(
    `Sample size` = sample_size,
    `Target censoring (%)`,
    `Average retained FPCs`,
    `Total retained FPCs`,
    `Total regression parameters`,
    `Number of events`,
    `EPV`
  )

scenario1_publication_table


write.csv(
  scenario1_publication_table,
  file.path(
    result_dir,
    "Scenario5_EPV_publication_table.csv"
  ),
  row.names = FALSE
)


