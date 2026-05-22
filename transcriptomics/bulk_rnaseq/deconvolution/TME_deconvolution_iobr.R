# TME deconvolution with IOBR
#
# Usage:
#   1) Edit config.yaml in the same directory
#   2) Run: Rscript TME_deconvolution_iobr.R
#   3) Or source this script in an R session

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  match <- grep(file_arg, cmd_args)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", cmd_args[match]))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  getwd()
}

resolve_path <- function(path, base_dir) {
  if (is.null(path) || path == "") {
    return(path)
  }
  if (grepl("^(/|[A-Za-z]:)", path)) {
    return(path)
  }
  file.path(base_dir, path)
}

load_config <- function(config_path = file.path(get_script_dir(), "config.yaml")) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("请先安装 yaml 包：install.packages('yaml')")
  }
  if (!file.exists(config_path)) {
    stop("未找到配置文件：", config_path)
  }
  yaml::read_yaml(config_path)
}

read_expression <- function(config) {
  input_path <- resolve_path(config$expression_file, get_script_dir())
  if (is.null(input_path) || !file.exists(input_path)) {
    stop("未找到表达矩阵：", input_path)
  }

  separator <- if (!is.null(config$separator)) config$separator else "\t"
  expr_df <- read.table(
    input_path,
    header = TRUE,
    sep = separator,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  gene_column <- if (!is.null(config$gene_column)) config$gene_column else 1
  if (is.character(gene_column)) {
    if (!gene_column %in% colnames(expr_df)) {
      stop("基因列不存在：", gene_column)
    }
    rownames(expr_df) <- expr_df[[gene_column]]
    expr_df[[gene_column]] <- NULL
  } else {
    rownames(expr_df) <- expr_df[[gene_column]]
    expr_df[[gene_column]] <- NULL
  }

  as.matrix(expr_df)
}

run_tme_deconvolution <- function(config = load_config()) {
  if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("请先安装 IOBR 包：\n  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n  BiocManager::install('IOBR')")
  }

  expr <- read_expression(config)

  output_dir <- resolve_path(
    if (!is.null(config$output_dir)) config$output_dir else "tme_deconvolution_results",
    get_script_dir()
  )
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  methods <- c(
    "cibersort",
    "cibersort_abs",
    "quanTIseq",
    "xcell",
    "mcp_counter",
    "epic",
    "timer",
    "estimate"
  )

  results <- list()
  for (method in methods) {
    message("Running deconvolution: ", method)
    deconv <- IOBR::deconvo_tme(
      eset = expr,
      method = method,
      arrays = FALSE
    )
    results[[method]] <- deconv

    output_path <- file.path(output_dir, paste0("tme_deconvolution_", method, ".tsv"))
    write.table(
      deconv,
      file = output_path,
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
  }

  saveRDS(results, file = file.path(output_dir, "tme_deconvolution_all_methods.rds"))

  results
}

results <- run_tme_deconvolution()
if (!is.null(results[[1]])) {
  print(head(results[[1]]))
}
