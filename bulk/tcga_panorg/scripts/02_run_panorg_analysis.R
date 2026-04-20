options(stringsAsFactors = FALSE)

get_script_dir_bootstrap <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[1]),
      winslash = "/",
      mustWork = FALSE
    )))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(
      sys.frames()[[1]]$ofile,
      winslash = "/",
      mustWork = FALSE
    )))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir_bootstrap()
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(project_root)

source(file.path(project_root, "R", "tcga_user_style_utils.R"))
source(file.path(project_root, "R", "tcga_user_style_prepare.R"))
source(file.path(project_root, "R", "tcga_user_style_panorg.R"))

check_required_packages(unique(c(required_user_style_packages(), required_panorg_packages())))

ensure_dir(project_path("results"))

cohort_cfg <- readr::read_csv(project_path("config", "tcga_cohorts.csv"), show_col_types = FALSE)

for (i in seq_len(nrow(cohort_cfg))) {
  cohort_row <- as.list(cohort_cfg[i, , drop = FALSE])
  message("Running analysis for ", cohort_row$cohort_id, " ...")
  bundle_file <- project_path("results", cohort_row$cohort_id, "rds", "tcga_bundle.rds")
  if (!file.exists(bundle_file)) {
    stop("Missing bundle file: ", bundle_file, ". Please run scripts/01_download_tcga_data.R first.")
  }
  bundle <- readRDS(bundle_file)
  run_panorg_pipeline_user_style(bundle = bundle, cohort_row = cohort_row)
}

message("All TCGA PANoRG analyses finished.")
