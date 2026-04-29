rm(list = ls())

cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(cmd_file) > 0) {
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", cmd_file[1]), fixed = TRUE), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
pipeline_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(pipeline_dir, "00_functions.R"))) {
  pipeline_dir <- normalizePath(script_dir, winslash = "/", mustWork = FALSE)
}

source(file.path(pipeline_dir, "00_functions.R"), local = FALSE)

step_subcluster_clustering()
