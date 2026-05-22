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

step_findmarkers()


sc = read_qs2("04_qs2\\02.sc_global_clustered.qs2")

pdac_main_markers_mouse <- list(
  Epithelial = c("Epcam", "Krt8", "Krt18"),
  Ductal = c("Krt19", "Mmp7", "Sox9"),
  Acinar = c("Prss1", "Ctrb1", "Pnlip"),
  Endocrine = c("Chgb", "Ins1", "Gcg"),

  Fibroblast = c("Col1a1", "Dcn", "Lum"),
  Pericyte = c("Rgs5", "Pdgfrb", "Kcnj8"),
  Endothelial = c("Pecam1", "Cdh5", "Kdr"),

  Myeloid = c("Lyz2", "Cd68", "Tyrobp"),
  T_cell = c("Cd3d", "Cd3e", "Cd2"),
  CD4_T_cell = c("Cd4", "Il7r", "Ccr7"),
  CD8_T_cell = c("Cd8a", "Gzmk", "Cd8b1"),
  NK_cell = c("Nkg7", "Gzmb", "Klrb1c"),

  B_cell = c("Cd79a", "Ms4a1", "Cd19"),
  Plasma = c("Jchain", "Mzb1", "Xbp1"),
  Mast = c("Kit", "Cpa3", "Tpsab1"),
  Platelet = c("Ppbp", "Pf4"),
  Schwann = c("Sox10", "Plp1", "Mpz")
)
DotPlot(sc, features = pdac_main_markers_mouse) +
  theme_bw() +
  RotatedAxis()
ggsave("03_FindMarkers\\dot.png",width = 20,height = 6)
