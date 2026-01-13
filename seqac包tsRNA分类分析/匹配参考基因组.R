## tRNA analysis in Seqpac 

# First create an annotation blank PAC with group means
anno(pac) <- anno(pac)[,1, drop=FALSE]
pac_trna <- PAC_summary(pac, norm = "cpm", type = "means", 
                        pheno_target=list("Group"), merge_pac = TRUE)

# Then reannotate only tRNA using the PAC_mapper function
input  <- "hg38-tRNAs.fa"  
# or use the provided fasta

map_object <- PAC_mapper(
  pac_trna,
  input = input,
  mismatches = 2,
  threads = 16,
  N_up = "NNN",
  N_down = "NNN",
  report_string = TRUE,
  override = TRUE
)

###---------------------------------------------------------------------
## Generate range types using a ss-file

# Download ss object from GtRNAdb 
# ("http://gtrnadb.ucsc.edu/")
ss_file <- "hg38-tRNAs-detailed.ss"
# Classify fragments according to loop cleavage (small loops are omitted)       
map_object_ss <- map_rangetype(map_object, type="ss", 
                               ss=ss_file, min_loop_width=4)            

# Remove reference tRNAs with no hits
map_object_ss <-  map_object_ss[
  !unlist(lapply(map_object_ss, function(x){x[[2]][1,1] == "no_hits"}))]

# Now we have quite comprehensive tRNA loop annotations
names(map_object_ss)
map_object_ss[[1]][[2]]

# Don't forget to check ?map_rangetype to obtain more options 


###---------------------------------------------------------------------
## Function classifying 5'-tRF, 5'halves, i-tRF, 3'-tRF, 3'halves

# Does all tRNAs have 3 loops?
table(unlist(lapply(map_object_ss, function(x){unique(x$Alignments$n_ssloop)})))

# Set tolerance for classification as a terminal tRF
tolerance <- 5  # 2 nucleotides from start or end of full-length tRNA)

### Important:
# We set N_up and N_down to "NNN" in the PAC_mapper step. To make sure
# that we have a tolerance that include the original tRNA sequence
# we set terminal= 2+3 (5).  
## tRNA classifying function 
# Apply the tRNA_class function and make a tRNA type column
pac_trna <- tRNA_class(pac_trna, map=map_object_ss, terminal=tolerance)
## 
## -- Anno target was specified, will retain: 29 of 9131 seqs.
anno(pac_trna)$type <- paste0(anno(pac_trna)$decoder, anno(pac_trna)$acceptor)
anno(pac_trna)[1:5, 1:4] 
