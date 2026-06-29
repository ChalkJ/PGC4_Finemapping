library(data.table)
BASE       <- "C:/Users/c1977426/OneDrive - Cardiff University/Documents/PGC4"
gene_cache <- readRDS(file.path(BASE, "r_files/gene_cache_all.rds"))

# Check which loci in ALL analyses have a prioritised-gene flag
for (aname in c("sczvsbp","sczvsbp1","sczandbp")) {
  cat("\n===", aname, "===\n")
  dir  <- file.path(BASE, "r_files", aname)
  loci <- fread(file.path(dir, "loci_input.txt"), header=FALSE,
                select=1:4, col.names=c("loci_id","chr","start","end"))
  ss_file <- switch(aname,
    sczvsbp  = "daner_scz_vs_bip_1025_noduppos_hetpva.gz",
    sczvsbp1 = "daner_scz_vs_bip1_1025_noduppos_hetpva.gz",
    sczandbp = "daner_scz_bip_vs_allcontrols_1025_noduppos_hetpva.gz")
  ss <- fread(file.path(BASE, "sumstats", ss_file))
  ss <- ss[HetPVa > 0.05]

  for (i in seq_len(nrow(loci))) {
    li    <- loci[i]
    locus <- sprintf("%03d", li$loci_id)
    rds   <- file.path(dir, "susie", sprintf("%s.susie_l5.rds", locus))
    if (!file.exists(rds)) next
    fit <- readRDS(rds)
    if (is.null(fit$sets$cs) || length(fit$sets$cs)==0) next
    chr_genes <- gene_cache[[paste0("chr", li$chr)]]
    if (is.null(chr_genes)) next

    for (cs_i in seq_along(fit$sets$cs)) {
      cs_snps <- fit$snp_ids[fit$sets$cs[[cs_i]]]
      cs_bp   <- ss[SNP %in% cs_snps & CHR == li$chr, BP]
      if (length(cs_bp)==0) next
      for (j in seq_len(nrow(chr_genes))) {
        if (all(cs_bp >= chr_genes$start_position[j] & cs_bp <= chr_genes$end_position[j])) {
          cat(sprintf("  Locus %d CS%d -> %s (%d-%d)\n",
              li$loci_id, cs_i, chr_genes$external_gene_name[j],
              chr_genes$start_position[j], chr_genes$end_position[j]))
          break
        }
      }
    }
  }
}
