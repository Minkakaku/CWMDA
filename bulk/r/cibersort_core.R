suppressPackageStartupMessages({
  library(parallel)
  library(e1071)
  library(preprocessCore)
})

core_alg_cibersort <- function(signature_matrix, mixture_vector) {
  svn_itor <- 3

  fit_one <- function(index_value) {
    nu_value <- c(0.25, 0.5, 0.75)[[index_value]]
    e1071::svm(
      signature_matrix,
      mixture_vector,
      type = "nu-regression",
      kernel = "linear",
      nu = nu_value,
      scale = FALSE
    )
  }

  fit_list <- if (Sys.info()[["sysname"]] == "Windows") {
    parallel::mclapply(seq_len(svn_itor), fit_one, mc.cores = 1)
  } else {
    parallel::mclapply(seq_len(svn_itor), fit_one, mc.cores = svn_itor)
  }

  rmse_values <- numeric(svn_itor)
  cor_values <- numeric(svn_itor)

  for (index_value in seq_len(svn_itor)) {
    weights <- t(fit_list[[index_value]]$coefs) %*% fit_list[[index_value]]$SV
    weights[weights < 0] <- 0
    weights <- weights / sum(weights)
    weighted_signature <- sweep(signature_matrix, MARGIN = 2, weights, "*")
    predicted <- apply(weighted_signature, 1, sum)
    rmse_values[[index_value]] <- sqrt(mean((predicted - mixture_vector)^2))
    cor_values[[index_value]] <- cor(predicted, mixture_vector)
  }

  best_index <- which.min(rmse_values)
  best_model <- fit_list[[best_index]]

  raw_weights <- t(best_model$coefs) %*% best_model$SV
  raw_weights[raw_weights < 0] <- 0
  normalized_weights <- raw_weights / sum(raw_weights)

  list(
    w = normalized_weights,
    mix_rmse = rmse_values[[best_index]],
    mix_r = cor_values[[best_index]]
  )
}

run_cibersort_permutations <- function(perm, signature_matrix, mixture_matrix) {
  iter <- 1
  gene_values <- as.list(data.matrix(mixture_matrix))
  dist_values <- numeric()

  while (iter <= perm) {
    random_mixture <- as.numeric(gene_values[sample(length(gene_values), nrow(signature_matrix))])
    random_mixture <- (random_mixture - mean(random_mixture)) / stats::sd(random_mixture)
    result <- core_alg_cibersort(signature_matrix, random_mixture)
    dist_values <- c(dist_values, result$mix_r)
    iter <- iter + 1
  }

  sort(dist_values)
}

run_cibersort_matrix <- function(signature_matrix, mixture_matrix, perm = 0, qn = TRUE, output_file = NULL) {
  signature_matrix <- data.matrix(signature_matrix)
  mixture_matrix <- data.matrix(mixture_matrix)

  signature_matrix <- signature_matrix[order(rownames(signature_matrix)), , drop = FALSE]
  mixture_matrix <- mixture_matrix[order(rownames(mixture_matrix)), , drop = FALSE]

  if (max(mixture_matrix, na.rm = TRUE) < 50) {
    mixture_matrix <- 2^mixture_matrix
  }

  if (isTRUE(qn)) {
    sample_names <- colnames(mixture_matrix)
    gene_names <- rownames(mixture_matrix)
    mixture_matrix <- preprocessCore::normalize.quantiles(mixture_matrix)
    colnames(mixture_matrix) <- sample_names
    rownames(mixture_matrix) <- gene_names
  }

  keep_mixture <- rownames(mixture_matrix) %in% rownames(signature_matrix)
  keep_signature <- rownames(signature_matrix) %in% rownames(mixture_matrix)
  mixture_matrix <- mixture_matrix[keep_mixture, , drop = FALSE]
  signature_matrix <- signature_matrix[keep_signature, , drop = FALSE]

  signature_matrix <- (signature_matrix - mean(signature_matrix)) / stats::sd(as.vector(signature_matrix))

  null_dist <- NULL
  if (perm > 0) {
    null_dist <- run_cibersort_permutations(perm, signature_matrix, mixture_matrix)
  }

  header <- c("Mixture", colnames(signature_matrix), "P-value", "Correlation", "RMSE")
  output <- NULL

  for (iter in seq_len(ncol(mixture_matrix))) {
    y <- mixture_matrix[, iter]
    y <- (y - mean(y)) / stats::sd(y)
    result <- core_alg_cibersort(signature_matrix, y)

    p_value <- if (!is.null(null_dist)) {
      1 - (which.min(abs(null_dist - result$mix_r)) / length(null_dist))
    } else {
      9999
    }

    row_value <- c(colnames(mixture_matrix)[[iter]], result$w, p_value, result$mix_r, result$mix_rmse)
    output <- if (is.null(output)) row_value else rbind(output, row_value)
  }

  if (!is.null(output_file)) {
    write.table(
      rbind(header, output),
      file = output_file,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }

  result_matrix <- rbind(header, output)
  result_matrix <- result_matrix[, -1, drop = FALSE]
  result_matrix <- result_matrix[-1, , drop = FALSE]
  result_matrix <- matrix(as.numeric(unlist(result_matrix)), nrow = nrow(result_matrix))
  rownames(result_matrix) <- colnames(mixture_matrix)
  colnames(result_matrix) <- c(colnames(signature_matrix), "P-value", "Correlation", "RMSE")
  result_matrix
}
