library(data.table)
library(openxlsx)

BASE       <- "C:/Users/c1977426/OneDrive - Cardiff University/Documents/PGC4"
GENE_CACHE <- readRDS(file.path(BASE, "r_files/gene_cache_v2.rds"))
OUT_FILE   <- file.path(BASE, "results/master_spreadsheet.xlsx")

analyses <- list(
  sczvsbp = list(
    label    = "SCZvsBP",
    dir      = file.path(BASE, "r_files/sczvsbp"),
    sumstats = file.path(BASE, "sumstats/daner_scz_vs_bip_1025_noduppos_hetpva.gz")
  ),
  sczvsbp1 = list(
    label    = "SCZvsBP1",
    dir      = file.path(BASE, "r_files/sczvsbp1"),
    sumstats = file.path(BASE, "sumstats/daner_scz_vs_bip1_1025_noduppos_hetpva.gz")
  ),
  sczandbp = list(
    label    = "SCZandBP",
    dir      = file.path(BASE, "r_files/sczandbp"),
    sumstats = file.path(BASE, "sumstats/daner_scz_bip_vs_allcontrols_1025_noduppos_hetpva.gz")
  )
)

# ── Helpers ──────────────────────────────────────────────────────────────────

load_vep <- function(vep_file) {
  if (!file.exists(vep_file)) return(list(gene = setNames(character(0), character(0)),
                                          csq  = setNames(character(0), character(0))))
  vep <- fread(vep_file)
  colnames(vep)[1] <- "SNP"
  rank_map <- c("HIGH"=1L, "MODERATE"=2L, "LOW"=3L, "MODIFIER"=4L, "-"=5L)
  vep[, rank := rank_map[IMPACT]]
  vep[is.na(rank), rank := 5L]
  best <- vep[order(SNP, rank)][!duplicated(SNP)]
  gene_col <- if ("SYMBOL" %in% names(best)) best$SYMBOL else rep("-", nrow(best))
  list(
    gene = setNames(gene_col,         best$SNP),
    csq  = setNames(best$Consequence, best$SNP)
  )
}

gene_by_pos <- function(snp_ids, chr_val, ss, chr_genes) {
  if (is.null(chr_genes) || nrow(chr_genes) == 0)
    return(rep(NA_character_, length(snp_ids)))
  bp_vec <- ss$BP[match(snp_ids, ss$SNP)]
  sapply(bp_vec, function(bp) {
    if (is.na(bp)) return(NA_character_)
    hits <- chr_genes[chr_genes$start_position <= bp & chr_genes$end_position >= bp, ]
    nms  <- hits$external_gene_name[!is.na(hits$external_gene_name) & hits$external_gene_name != ""]
    if (length(nms) == 0) return(NA_character_)
    paste(unique(nms), collapse = "; ")
  })
}

annotate_snps <- function(snp_ids, chr_val, ss, vep, chr_genes) {
  # Gene: VEP SYMBOL, fall back to positional
  gene <- ifelse(snp_ids %in% names(vep$gene), vep$gene[snp_ids], NA_character_)
  blank <- is.na(gene) | gene == "-" | gene == ""
  if (any(blank))
    gene[blank] <- gene_by_pos(snp_ids[blank], chr_val, ss, chr_genes)

  # Consequence
  csq <- ifelse(snp_ids %in% names(vep$csq), vep$csq[snp_ids], NA_character_)

  list(gene = gene, csq = csq)
}

best_cred_file <- function(locus_dir, locus) {
  creds <- sort(list.files(locus_dir, pattern = paste0("^", locus, "\\.cred[0-9]+$"),
                           full.names = TRUE))
  best <- NULL; best_pr <- -1
  for (cf in creds) {
    hdr <- readLines(cf, n = 2)
    pr  <- suppressWarnings(as.numeric(regmatches(hdr[1], regexpr("[0-9.eE+\\-]+$", hdr[1]))))
    if (!is.na(pr) && pr > best_pr) { best_pr <- pr; best <- cf }
  }
  best
}

# ── Build workbook ────────────────────────────────────────────────────────────

wb <- createWorkbook()

header_style <- createStyle(fontColour = "#FFFFFF", fgFill = "#2C3E50",
                             halign = "CENTER", textDecoration = "Bold",
                             border = "Bottom", borderColour = "#FFFFFF")

for (aname in names(analyses)) {
  a   <- analyses[[aname]]
  cat("Processing:", a$label, "\n")

  ss   <- fread(a$sumstats)
  ss   <- ss[HetPVa > 0.05]
  ss[, zscore := log(OR) / SE]

  loci <- fread(file.path(a$dir, "loci_input.txt"), header = FALSE,
                select = 1:4, col.names = c("locus_id", "chr", "start", "end"))

  vep      <- load_vep(file.path(a$dir, "vep_annotation.txt"))
  fm_dir   <- file.path(a$dir, "finemap")
  susie_dir <- file.path(a$dir, "susie")

  susie_rows <- list()
  fm_rows    <- list()

  for (i in seq_len(nrow(loci))) {
    locus_id  <- loci$locus_id[i]
    locus     <- sprintf("%03d", locus_id)
    chr_val   <- loci$chr[i]
    chr_genes <- GENE_CACHE[[paste0("chr", chr_val)]]

    # ── SuSiE ──────────────────────────────────────────────────────────────
    rds <- file.path(susie_dir, sprintf("%s.susie_l5.rds", locus))
    if (file.exists(rds)) {
      fit <- readRDS(rds)
      if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
        for (cs_i in seq_along(fit$sets$cs)) {
          idx     <- fit$sets$cs[[cs_i]]
          snps    <- fit$snp_ids[idx]
          pips    <- fit$pip[idx]
          ann     <- annotate_snps(snps, chr_val, ss, vep, chr_genes)
          ss_join <- ss[match(snps, SNP)]
          susie_rows[[length(susie_rows) + 1]] <- data.table(
            Locus       = locus_id,
            CS          = cs_i,
            PIP         = round(pips, 6),
            Gene        = ann$gene,
            Consequence = ann$csq,
            ss_join
          )
        }
      }
    }

    # ── FINEMAP ────────────────────────────────────────────────────────────
    locus_fm_dir <- file.path(fm_dir, locus)
    snp_file     <- file.path(locus_fm_dir, paste0(locus, ".snp"))
    cf           <- best_cred_file(locus_fm_dir, locus)

    if (file.exists(snp_file) && !is.null(cf)) {
      fm_snps <- fread(snp_file)
      clines  <- readLines(cf)
      skip    <- grep("^index", clines) - 1
      if (length(skip) > 0) {
        cdat     <- tryCatch(fread(cf, skip = skip[1]), error = function(e) NULL)
        cs_cols  <- if (!is.null(cdat)) grep("^cred[0-9]+$", names(cdat), value = TRUE) else character(0)

        for (cs_i in seq_along(cs_cols)) {
          snps <- cdat[[cs_cols[cs_i]]]
          snps <- snps[!is.na(snps) & snps != ""]
          if (length(snps) == 0) next
          pips    <- fm_snps$prob[match(snps, fm_snps$rsid)]
          ann     <- annotate_snps(snps, chr_val, ss, vep, chr_genes)
          ss_join <- ss[match(snps, SNP)]
          fm_rows[[length(fm_rows) + 1]] <- data.table(
            Locus       = locus_id,
            CS          = cs_i,
            PIP         = round(pips, 6),
            Gene        = ann$gene,
            Consequence = ann$csq,
            ss_join
          )
        }
      }
    }
  }

  # Write SuSiE sheet
  sheet_susie <- paste0(a$label, "_SuSiE")
  addWorksheet(wb, sheet_susie)
  if (length(susie_rows) > 0) {
    dt_susie <- rbindlist(susie_rows, fill = TRUE)
    dt_susie <- dt_susie[order(Locus, CS, -PIP)]
    writeData(wb, sheet_susie, dt_susie, headerStyle = header_style)
    setColWidths(wb, sheet_susie, cols = seq_len(ncol(dt_susie)), widths = "auto")
    freezePane(wb, sheet_susie, firstRow = TRUE)
  } else {
    writeData(wb, sheet_susie, data.frame(Note = "No credible sets found"))
  }
  cat("  SuSiE:", if (length(susie_rows) > 0) nrow(rbindlist(susie_rows, fill=TRUE)) else 0, "SNPs\n")

  # Write FINEMAP sheet
  sheet_fm <- paste0(a$label, "_FINEMAP")
  addWorksheet(wb, sheet_fm)
  if (length(fm_rows) > 0) {
    dt_fm <- rbindlist(fm_rows, fill = TRUE)
    dt_fm <- dt_fm[order(Locus, CS, -PIP)]
    writeData(wb, sheet_fm, dt_fm, headerStyle = header_style)
    setColWidths(wb, sheet_fm, cols = seq_len(ncol(dt_fm)), widths = "auto")
    freezePane(wb, sheet_fm, firstRow = TRUE)
  } else {
    writeData(wb, sheet_fm, data.frame(Note = "No credible sets found"))
  }
  cat("  FINEMAP:", if (length(fm_rows) > 0) nrow(rbindlist(fm_rows, fill=TRUE)) else 0, "SNPs\n")
}

saveWorkbook(wb, OUT_FILE, overwrite = TRUE)
cat("\nSaved:", OUT_FILE, "\n")
