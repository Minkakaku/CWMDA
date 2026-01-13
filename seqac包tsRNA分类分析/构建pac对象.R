str(count_list, max.level = 1)


Sample_ID <- colnames(count_list$counts)
pheno_fl <- data.frame(
  Sample_ID = c(
    "NC1", "NC2", "NC3",
    "TNBCexo1", "TNBCexo2", "TNBCexo3",
    "Her2exo1", "Her2exo2", "Her2exo3",
    "HRexo1", "HRexo2", "HRexo3"
  ),
  Group = c(
    rep("NC", 3),
    rep("TNBCexo", 3),
    rep("Her2exo", 3),
    rep("HRexo", 3)
  ),
  stringsAsFactors = FALSE
)


pheno <- make_pheno(
  pheno = pheno_fl,
  counts = count_list$counts
)

###--------------------------------------------------------------------- 
## Generate PAC object (default S4)
pac_master <- make_PAC(pheno=pheno, counts=count_list$counts)
pac_master
# 查看加载的对象
pac=pac_master

pac <- PAC_norm(pac)
pac
