library(dplyr)
library(purrr)
library(readr)
library(tidyr)


result_dir <- "."


rho_values <- 0.5
censoring_values <- 0.2

sample_size_values <- c(  200, 300, 500, 1000, 2000)

q_values <- c(50,100,200,300)

runs <- 1:100


q_true <- 2
mfunctions_true <- 2






file_information <- tidyr::crossing(
  simulation_run = runs,
  rho = rho_values,
  target_censoring = censoring_values,
  sample_size = sample_size_values,
  q = q_values
) %>%
  mutate(
    
    mfunctions = q,
    
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
    
    filepath = file.path(
      result_dir,
      filename
    ),
    
    file_exists = file.exists(filepath)
  )



nrow(file_information)


table(file_information$file_exists)



missing_files <- file_information %>%
  filter(!file_exists)

nrow(missing_files)

missing_files



missing_file_summary <- missing_files %>%
  count(
    sample_size,
    q,
    mfunctions,
    name = "number_missing_files"
  ) %>%
  arrange(
    sample_size,
    q
  )

missing_file_summary






safe_read_EPV <- function(filepath,
                          expected_run,
                          expected_rho,
                          expected_censoring,
                          expected_n,
                          expected_q,
                          expected_mfunctions) {
  
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
      
      warning(
        "文件没有数据：",
        basename(filepath)
      )
      
      return(NULL)
    }
    
    
    
    
    if (nrow(dat) > 1) {
      
      warning(
        "文件包含多行数据：",
        basename(filepath),
        "；行数 = ",
        nrow(dat)
      )
    }
    
    
    
    dat %>%
      mutate(
        source_file = basename(filepath),
        expected_run = expected_run,
        expected_rho = expected_rho,
        expected_censoring = expected_censoring,
        expected_sample_size = expected_n,
        expected_q = expected_q,
        expected_mfunctions = expected_mfunctions,
        read_success = TRUE
      )
    
  }, error = function(e) {
    
    message(
      "读取失败：",
      basename(filepath),
      "；错误信息：",
      conditionMessage(e)
    )
    
    return(NULL)
  })
}






existing_file_information <- file_information %>%
  filter(file_exists)


scenario4_raw <- purrr::pmap_dfr(
  list(
    filepath =
      existing_file_information$filepath,
    
    expected_run =
      existing_file_information$simulation_run,
    
    expected_rho =
      existing_file_information$rho,
    
    expected_censoring =
      existing_file_information$target_censoring,
    
    expected_n =
      existing_file_information$sample_size,
    
    expected_q =
      existing_file_information$q,
    
    expected_mfunctions =
      existing_file_information$mfunctions
  ),
  safe_read_EPV
)



nrow(scenario4_raw)

head(scenario4_raw)






numeric_variables <- c(
  "simulation_run",
  "rho",
  "sample_size",
  "target_censoring",
  "actual_censoring",
  "q",
  "q_true",
  "mfunctions",
  "mfunctions_true",
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


scenario4_raw <- scenario4_raw %>%
  mutate(
    across(
      any_of(numeric_variables),
      as.numeric
    )
  )






parameter_check <- scenario4_raw %>%
  mutate(
    run_match =
      simulation_run == expected_run,
    
    rho_match =
      abs(rho - expected_rho) < 1e-10,
    
    censoring_match =
      abs(
        target_censoring -
          expected_censoring
      ) < 1e-10,
    
    sample_size_match =
      sample_size == expected_sample_size,
    
    
    q_match =
      if ("q" %in% names(scenario4_raw)) {
        q == expected_q
      } else {
        TRUE
      },
    
    
    mfunctions_match =
      if ("mfunctions" %in% names(scenario4_raw)) {
        mfunctions == expected_mfunctions
      } else {
        TRUE
      },
    
    
    q_mfunctions_match =
      expected_q == expected_mfunctions
  )


parameter_mismatch <- parameter_check %>%
  filter(
    !run_match |
      !rho_match |
      !censoring_match |
      !sample_size_match |
      !q_match |
      !mfunctions_match |
      !q_mfunctions_match
  )



nrow(parameter_mismatch)

parameter_mismatch






scenario4_raw <- scenario4_raw %>%
  mutate(
    calculated_total_parameters =
      number_scalar_parameters +
      number_derivative_parameters +
      number_functional_parameters,
    
    total_parameter_match =
      total_regression_parameters ==
      calculated_total_parameters,
    
    calculated_EPV =
      number_events /
      total_regression_parameters,
    
    EPV_difference =
      EPV - calculated_EPV
  )


parameter_EPV_check <- scenario4_raw %>%
  filter(
    !total_parameter_match |
      abs(EPV_difference) > 1e-10
  ) %>%
  select(
    source_file,
    expected_sample_size,
    expected_q,
    expected_mfunctions,
    total_regression_parameters,
    calculated_total_parameters,
    EPV,
    calculated_EPV,
    EPV_difference
  )


parameter_EPV_check






scenario4_run_check <- scenario4_raw %>%
  group_by(
    rho,
    sample_size,
    target_censoring,
    expected_q,
    expected_mfunctions
  ) %>%
  summarise(
    number_records = n(),
    
    number_unique_runs =
      n_distinct(simulation_run),
    
    minimum_run =
      min(simulation_run, na.rm = TRUE),
    
    maximum_run =
      max(simulation_run, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(
    sample_size,
    expected_q
  )


scenario4_run_check



incomplete_settings <- scenario4_run_check %>%
  filter(
    number_unique_runs != 100
  )


incomplete_settings






scenario4_summary_numeric <- scenario4_raw %>%
  group_by(
    rho,
    sample_size,
    target_censoring,
    q = expected_q,
    mfunctions = expected_mfunctions
  ) %>%
  summarise(
    number_simulations = n(),
    
    actual_censoring_mean =
      mean(
        actual_censoring,
        na.rm = TRUE
      ),
    
    actual_censoring_sd =
      sd(
        actual_censoring,
        na.rm = TRUE
      ),
    
    average_retained_fpc_mean =
      mean(
        average_retained_fpc,
        na.rm = TRUE
      ),
    
    average_retained_fpc_sd =
      sd(
        average_retained_fpc,
        na.rm = TRUE
      ),
    
    total_retained_fpc_mean =
      mean(
        total_retained_fpc,
        na.rm = TRUE
      ),
    
    total_retained_fpc_sd =
      sd(
        total_retained_fpc,
        na.rm = TRUE
      ),
    
    total_regression_parameters_mean =
      mean(
        total_regression_parameters,
        na.rm = TRUE
      ),
    
    total_regression_parameters_sd =
      sd(
        total_regression_parameters,
        na.rm = TRUE
      ),
    
    number_events_mean =
      mean(
        number_events,
        na.rm = TRUE
      ),
    
    number_events_sd =
      sd(
        number_events,
        na.rm = TRUE
      ),
    
    EPV_mean =
      mean(
        EPV,
        na.rm = TRUE
      ),
    
    EPV_sd =
      sd(
        EPV,
        na.rm = TRUE
      ),
    
    EPV_min =
      min(
        EPV,
        na.rm = TRUE
      ),
    
    EPV_max =
      max(
        EPV,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  arrange(
    sample_size,
    q
  )


scenario4_summary_numeric



scenario4_summary_numeric %>%
  filter(
    q != mfunctions
  )






scenario4_publication_table <- scenario4_summary_numeric %>%
  mutate(
    `Target censoring (%)` =
      sprintf(
        "%.0f",
        100 * target_censoring
      ),
    
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
    `p = q` = q,
    `Average retained FPCs`,
    `Total retained FPCs`,
    `Total regression parameters`,
    `Number of events`,
    `EPV`
  )


scenario4_publication_table






write.csv(
  scenario4_publication_table,
  file.path(
    result_dir,
    "Scenario8_EPV_publication_table.csv"
  ),
  row.names = FALSE
)





