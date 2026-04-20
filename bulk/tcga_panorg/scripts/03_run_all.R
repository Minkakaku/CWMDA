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

rscript_bin <- file.path(R.home("bin"), "Rscript")

run_step <- function(script_name) {
  output <- system2(
    rscript_bin,
    c(file.path(project_root, "scripts", script_name)),
    stdout = TRUE,
    stderr = TRUE
  )
  if (length(output) > 0) {
    cat(paste(output, collapse = "\n"), "\n")
  }
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop("Step failed: ", script_name)
  }
}

run_step("01_download_tcga_data.R")
run_step("02_run_panorg_analysis.R")

message("TCGA all-in-one run finished.")
