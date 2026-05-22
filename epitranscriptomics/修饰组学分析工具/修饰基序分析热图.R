suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
})

# ===========================
# INPUT
# ===========================
fasta_file <- "proteome.fasta"          # 物种蛋白 FASTA（背景构建用）
fg_file    <- "motif.txt"    # 每行一个 13-mer
out_pdf    <- "琥珀酰化.pdf"


center_residue <- "K"   # K / S / T / Y

# 13-mer => ±6
k <- 6
L <- 2*k + 1
positions <- (-k):k
pos_labels <- ifelse(positions > 0, paste0("+", positions), as.character(positions))

AA <- c("A","C","D","E","F","G","H","I","K","L",
        "M","N","P","Q","R","S","T","V","W","Y")

# 显著性阈值（对应你 Methods 的 1e-06）
p_cut <- 0.05

# 色标截断，避免极端值把图拉平
lim <- 10

# ===========================
# helper: count matrix
# ===========================
calc_pos_count <- function(mat, AA, pos_labels){
  L <- ncol(mat)
  out <- matrix(0, nrow=length(AA), ncol=L,
                dimnames=list(AA, pos_labels))
  for (j in seq_len(L)){
    tab <- table(mat[,j])
    for (aa in AA){
      out[aa,j] <- ifelse(aa %in% names(tab), tab[[aa]], 0)
    }
  }
  out
}

# ===========================
# 1) read foreground 13-mers
# ===========================
fg <- readLines(fg_file)
fg <- fg[nzchar(fg)]
fg <- toupper(fg)

stopifnot(all(nchar(fg) == L))
stopifnot(all(substr(fg, k+1, k+1) == center_residue))

fg_mat <- do.call(rbind, strsplit(fg, ""))

# ===========================
# 2) build background from proteome
# ===========================
fa <- readAAStringSet(fasta_file)

extract_centered_windows <- function(seq, center_residue="K", k=6){
  s <- as.character(seq)
  n <- nchar(s)
  if (n < (2*k+1)) return(character(0))
  chars <- strsplit(s, "")[[1]]
  idx <- which(chars == center_residue)
  idx <- idx[idx > k & idx <= (n - k)]
  if (length(idx) == 0) return(character(0))
  
  vapply(idx, function(i){
    paste0(chars[(i-k):(i+k)], collapse = "")
  }, character(1))
}

bg <- unlist(lapply(fa, extract_centered_windows,
                    center_residue=center_residue, k=k),
             use.names = FALSE)

if (length(bg) == 0) stop("Background windows extracted = 0. Check fasta or center_residue.")

bg <- toupper(bg)
bg <- bg[nchar(bg) == L]
bg <- bg[substr(bg, k+1, k+1) == center_residue]
bg_mat <- do.call(rbind, strsplit(bg, ""))

# ===========================
# 3) counts
# ===========================
cnt_fg <- calc_pos_count(fg_mat, AA, pos_labels)
cnt_bg <- calc_pos_count(bg_mat, AA, pos_labels)

N_fg <- nrow(fg_mat)
N_bg <- nrow(bg_mat)

# ===========================
# 4) signed -log10(P) via Fisher test
# ===========================
score <- matrix(0, nrow=length(AA), ncol=length(pos_labels),
                dimnames=list(AA, pos_labels))

for (aa in AA) {
  for (pos in pos_labels) {
    
    a <- cnt_fg[aa, pos]
    c <- cnt_bg[aa, pos]
    b <- N_fg - a
    d <- N_bg - c
    
    mat2x2 <- matrix(c(a, b, c, d), nrow=2, byrow=TRUE)
    
    p <- fisher.test(mat2x2, alternative = "two.sided")$p.value
    
    # direction by frequency
    p_fg <- a / N_fg
    p_bg <- c / N_bg
    sgn <- ifelse(p_fg >= p_bg, 1, -1)
    
    # thresholding: non-significant -> 0 (white)
    if (p <= p_cut) {
      score[aa, pos] <- sgn * (-log10(p + 1e-300))
    } else {
      score[aa, pos] <- 0
    }
  }
}

# 0 位点（中心修饰位点）通常不解释：直接设为 0（白色）
score[, "0"] <- 0

# 截断到 [-lim, lim]，避免极端爆色
score[score > lim] <- lim
score[score < -lim] <- -lim

# ===========================
# 5) ggplot
# ===========================
df <- as.data.frame(score) %>%
  rownames_to_column("AA") %>%
  pivot_longer(-AA, names_to="Position", values_to="Score")

aa_order_topdown <- c("Y","W","V","T","S","R","Q","P","N","M",
                      "L","K","I","H","G","F","E","D","C","A")
df$AA <- factor(df$AA, levels = rev(aa_order_topdown))
df$Position <- factor(df$Position, levels = pos_labels)

p <- ggplot(df, aes(x=Position, y=AA, fill=Score)) +
  geom_tile(color="grey85", linewidth=0.35) +
  scale_fill_gradient2(
    low = "#2DB7A3",
    mid = "white",
    high = "#E94B4B",
    midpoint = 0,
    limits = c(-lim, lim),
    breaks = c(-lim, -6, 0, 6, lim),
    name = NULL
  ) +
  labs(x="Position", y="Amino acid") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm"),
    legend.key.width  = unit(0.6, "cm")
  )
p

ggsave(out_pdf, p, width=11, height=7)
cat("✅ Saved:", out_pdf, "\n")
cat("Foreground n =", N_fg, "Background n =", N_bg, "\n")
cat("P-cut =", p_cut, "  Color limit =", lim, "\n")

