# ============================================================================
# Phase 2: Test LD generation pipeline against tutorial's norcloz_chr2.ld
# Steps:
#   1. Run PLINK (via WSL) to generate pairwise LD from CLOZUK2 chr2
#   2. Run LDmerge (single-cohort) to produce square matrix
#   3. Compare against tutorial's pre-computed norcloz_chr2.ld
# ============================================================================

library(data.table)
library(readr)
library(dplyr)
library(purrr)
library(reshape2)
library(Matrix)

TESTDIR  <- "C:/Users/c1977426/OneDrive - Cardiff University/Documents/sczvsbp/met588_finemapping"
RESULTDIR <- file.path(TESTDIR, "results")
WORKDIR  <- file.path(TESTDIR, "test_phase2_work")
WSL_WORKDIR <- "/home/c1977426/sczvsbp/met588_finemapping/test_phase2_work"
WSL_TESTDIR <- "/home/c1977426/sczvsbp/met588_finemapping"
LDMERGE  <- "C:/Users/c1977426/OneDrive - Cardiff University/Documents/sczvsbp/new_scripts/LDmerge_v2.R"
RSCRIPT  <- "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"
LOCUS    <- "norcloz_chr2_test"

dir.create(WORKDIR, showWarnings = FALSE)

# ── 1. Run PLINK pipeline in WSL ──────────────────────────────────────────────
cat("=== 1. RUNNING PLINK PIPELINE (WSL) ===\n")
sh_wsl <- paste0(WSL_TESTDIR, "/test_phase2_setup.sh")
cmd <- sprintf("wsl bash %s", sh_wsl)
cat("  Command:", cmd, "\n")
t0 <- proc.time()
rc <- system(cmd, intern = FALSE)
cat(sprintf("  Exit code: %d | Time: %.1fs\n", rc, (proc.time() - t0)["elapsed"]))

if (rc != 0) {
  cat("  PLINK step failed — aborting\n")
  quit(status = 1)
}

# ── 2. Run LDmerge (single-cohort) ────────────────────────────────────────────
# We patch the stop() thresholds for single-cohort use via a temp modified copy
cat("\n=== 2. RUNNING LDmerge (single-cohort) ===\n")

ldmerge_src  <- readLines(LDMERGE)
ldmerge_text <- paste(ldmerge_src, collapse = "\n")

# Patch 1: lower "fewer than two" stop thresholds to allow single cohort
ldmerge_text <- gsub("length\\(corfile\\) < 2", "length(corfile) < 1", ldmerge_text)
ldmerge_text <- gsub("length\\(famfile\\) < 2", "length(famfile) < 1", ldmerge_text)

# Patch 2: replace the apply(LDframe[,c(-1,-2)], ...) block with a version that
# handles single cohort (1 correlation column) without the dim() error
old_apply <- "LDframe$wcor <-\n  apply(LDframe[, c(-1, -2)],\n        1,\n        weighted.mean,\n        w = as.numeric(neffvec),\n        na.rm = T)"
new_apply  <- "LDframe$wcor <- local({\n    cb <- LDframe[, c(-1, -2), drop = FALSE]\n    if (ncol(cb) == 1) cb[[1]] else apply(cb, 1, weighted.mean, w = as.numeric(neffvec), na.rm = TRUE)\n  })"
if (!grepl(old_apply, ldmerge_text, fixed = TRUE)) {
  cat("  WARNING: apply patch string not found in LDmerge — patch may have failed!\n")
  cat("  Trying alternative with leading whitespace...\n")
}
ldmerge_text <- gsub(old_apply, new_apply, ldmerge_text, fixed = TRUE)
cat(sprintf("  Patch 2 applied: %s\n",
    ifelse(grepl("if (ncol(cb) == 1)", ldmerge_text, fixed=TRUE), "YES", "NO - FAILED")))

ldmerge_tmp <- file.path(WORKDIR, "LDmerge_patched.R")
# Write with Unix line endings so Rscript reads it cleanly
con <- file(ldmerge_tmp, open = "wb"); writeLines(strsplit(ldmerge_text, "\n")[[1]], con, sep = "\n"); close(con)
cat("  Patched LDmerge written to:", ldmerge_tmp, "\n")

# Run LDmerge from the work directory using Windows Rscript (R not installed in WSL)
# setwd() propagates to child processes spawned by system()
RSCRIPT_WIN <- "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"
old_wd <- getwd()
setwd(WORKDIR)
ldmerge_cmd <- sprintf('"%s" --vanilla "%s" %s PLINK METAL',
                       RSCRIPT_WIN,
                       file.path(WORKDIR, "LDmerge_patched.R"),
                       LOCUS)
cat("  Command:", ldmerge_cmd, "\n")
t0 <- proc.time()
rc2 <- system(ldmerge_cmd, intern = FALSE)
setwd(old_wd)
cat(sprintf("  Exit code: %d | Time: %.1fs\n", rc2, (proc.time() - t0)["elapsed"]))

if (rc2 != 0) {
  cat("  LDmerge failed — check output above\n")
  quit(status = 1)
}

# ── 3. Compare matrices ───────────────────────────────────────────────────────
cat("\n=== 3. COMPARING LD MATRICES ===\n")

our_ld_file <- file.path(WORKDIR, paste0(LOCUS, ".ld"))
ref_ld_file <- file.path(RESULTDIR, "norcloz_chr2.ld")

if (!file.exists(our_ld_file)) {
  cat("  ERROR: LDmerge output not found:", our_ld_file, "\n")
  quit(status = 1)
}

cat("  Loading tutorial LD matrix...\n")
LD_ref <- as.matrix(fread(ref_ld_file, header = FALSE))
cat(sprintf("  Tutorial matrix: %d x %d\n", nrow(LD_ref), ncol(LD_ref)))

cat("  Loading pipeline LD matrix...\n")
LD_our <- as.matrix(fread(our_ld_file, header = FALSE))
cat(sprintf("  Pipeline matrix: %d x %d\n", nrow(LD_our), ncol(LD_our)))

# Load SNP lists for row/col labelling
snp_our <- readLines(file.path(WORKDIR, paste0(LOCUS, ".snp.log")))
cat(sprintf("  Pipeline SNP list: %d SNPs\n", length(snp_our)))

# Reference SNP order from .z file
z_dt   <- fread(file.path(RESULTDIR, "norcloz_chr2.z"))
snp_ref <- z_dt$rsid

# Check dimensions match
if (nrow(LD_ref) != nrow(LD_our)) {
  cat(sprintf("  DIMENSION MISMATCH: ref=%d, ours=%d\n", nrow(LD_ref), nrow(LD_our)))
  cat("  SNPs in ref but not pipeline:", sum(!snp_ref %in% snp_our), "\n")
  cat("  SNPs in pipeline but not ref:", sum(!snp_our %in% snp_ref), "\n")
  # Try to find common SNPs and subset
  common <- intersect(snp_ref, snp_our)
  cat(sprintf("  Common SNPs: %d — subsetting for comparison\n", length(common)))
  rownames(LD_ref) <- colnames(LD_ref) <- snp_ref
  rownames(LD_our) <- colnames(LD_our) <- snp_our
  LD_ref <- LD_ref[common, common]
  LD_our <- LD_our[common, common]
} else {
  rownames(LD_ref) <- colnames(LD_ref) <- snp_ref
  rownames(LD_our) <- colnames(LD_our) <- snp_our
}

# Numerical comparison (upper triangle only to avoid double-counting)
ut <- upper.tri(LD_ref)
r_ref <- LD_ref[ut]
r_our <- LD_our[ut]

diff   <- r_our - r_ref
pearson_r <- cor(r_ref, r_our)

cat(sprintf("\n  Pairs compared (upper triangle): %d\n", length(r_ref)))
cat(sprintf("  Pearson r (ref vs pipeline): %.6f\n", pearson_r))
cat(sprintf("  Mean |diff|: %.6f\n", mean(abs(diff))))
cat(sprintf("  Max  |diff|: %.6f\n", max(abs(diff))))
cat(sprintf("  RMSE:        %.6f\n", sqrt(mean(diff^2))))

# Sign discordance (pairs where sign differs)
sign_flip <- sum(sign(r_ref) != sign(r_our) & abs(r_ref) > 0.1)
cat(sprintf("  Sign flips (|r_ref|>0.1): %d\n", sign_flip))

# Spot-check: largest discrepancies
worst_idx <- order(abs(diff), decreasing = TRUE)[1:5]
worst_pos <- which(ut, arr.ind = TRUE)[worst_idx, ]
cat("\n  Top 5 worst pairs:\n")
cat(sprintf("  %-20s %-20s ref=%7.4f  ours=%7.4f  diff=%7.4f\n",
    "SNP_A", "SNP_B", 0, 0, 0))
for (k in seq_len(5)) {
  i <- worst_pos[k, 1]
  j <- worst_pos[k, 2]
  cat(sprintf("  %-20s %-20s ref=%7.4f  ours=%7.4f  diff=%7.4f\n",
      rownames(LD_ref)[i], colnames(LD_ref)[j],
      LD_ref[i, j], LD_our[i, j], diff[worst_idx[k]]))
}

cat("\n=== VERDICT ===\n")
if (pearson_r > 0.9999 && max(abs(diff)) < 0.01) {
  cat("PASS: matrices are essentially identical (r>0.9999, max|diff|<0.01)\n")
} else if (pearson_r > 0.999) {
  cat(sprintf("CLOSE: high correlation (r=%.5f) but some differences\n", pearson_r))
  cat("  Check sign flips and worst pairs above\n")
} else {
  cat(sprintf("FAIL: matrices differ substantially (r=%.4f)\n", pearson_r))
  cat("  Likely causes: allele orientation mismatch, SNP order, LD window\n")
}

cat("\n=== DONE ===\n")
