options(stringsAsFactors = FALSE)

source("R/tcga_user_style_utils.R")
source("R/tcga_user_style_prepare.R")

check_required_packages(required_user_style_packages())

ensure_dir(project_path("gdc_cache"))
ensure_dir(project_path("results"))

cohort_cfg <- readr::read_csv(project_path("config", "tcga_cohorts.csv"), show_col_types = FALSE)
candidate_cfg <- readr::read_csv(project_path("config", "molecular_entity_candidates.csv"), show_col_types = FALSE)
panorgs <- read_panorgs(project_path("config", "panorgs.txt"))

for (i in seq_len(nrow(cohort_cfg))) {
  cohort_row <- as.list(cohort_cfg[i, , drop = FALSE])
  message("========== ", cohort_row$cohort_id, " ==========")
  bundle <- prepare_one_tcga_cohort(cohort_row = cohort_row, panorgs = panorgs, candidate_cfg = candidate_cfg)
  write_tcga_bundle_outputs(bundle, project_path("results", cohort_row$cohort_id))
}

message("TCGA download and preprocessing finished.")
