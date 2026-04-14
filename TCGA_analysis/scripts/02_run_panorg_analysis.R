options(stringsAsFactors = FALSE)

source("R/tcga_user_style_utils.R")
source("R/tcga_user_style_prepare.R")
source("R/tcga_user_style_panorg.R")

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
