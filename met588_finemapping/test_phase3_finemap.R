# ============================================================================
# Phase 3: Run FINEMAP k=5 on pipeline-generated LD (test_phase2_work)
# Goal: verify the corrected PLINK pipeline resolves to 2 credible sets
# ============================================================================

library(data.table)

TESTDIR  <- "C:/Users/c1977426/OneDrive - Cardiff University/Documents/sczvsbp/met588_finemapping"
RESULTDIR  <- file.path(TESTDIR, "results")
WORKDIR  <- file.path(TESTDIR, "test_phase2_work")
OUTDIR   <- file.path(TESTDIR, "test_phase3_output")
WSL_OUTDIR <- "/home/c1977426/sczvsbp/met588_finemapping/test_phase3_output"
FINEMAP_BIN <- "/home/c1977426/bin/finemap"
N   <- 2989
K   <- 5

dir.create(OUTDIR, showWarnings = FALSE)

# ── 1. Load z-file and pipeline SNP list ─────────────────────────────────────
cat("=== 1. LOADING DATA ===\n")
z_dt     <- fread(file.path(RESULTDIR, "norcloz_chr2.z"))
z_dt[, zscore := beta / se]
pipe_snps <- readLines(file.path(WORKDIR, "norcloz_chr2_test.snp.log"))

cat(sprintf("  Tutorial z-file: %d SNPs\n", nrow(z_dt)))
cat(sprintf("  Pipeline LD SNPs: %d\n", length(pipe_snps)))

common <- intersect(z_dt$rsid, pipe_snps)
cat(sprintf("  Common SNPs: %d\n", length(common)))

if (length(common) < 10) {
  cat("  ERROR: too few common SNPs — aborting\n")
  quit(status = 1)
}

# Subset z to common SNPs in the pipeline order
z_sub <- z_dt[rsid %in% common]
# Reorder to match pipeline SNP order
z_sub <- z_sub[match(pipe_snps[pipe_snps %in% common], rsid)]
cat(sprintf("  z-score range: %.2f to %.2f | max |z|=%.2f (%s)\n",
    min(z_sub$zscore), max(z_sub$zscore),
    max(abs(z_sub$zscore)), z_sub$rsid[which.max(abs(z_sub$zscore))]))

# ── 2. Load pipeline LD, subset to common SNPs in same order ─────────────────
cat("\n=== 2. LOADING PIPELINE LD ===\n")
ld_file <- file.path(WORKDIR, "norcloz_chr2_test.ld")
cat(sprintf("  File size: %.1f MB\n", file.size(ld_file) / 1e6))
LD_full <- as.matrix(fread(ld_file, header = FALSE))
rownames(LD_full) <- colnames(LD_full) <- pipe_snps
cat(sprintf("  Full matrix: %d x %d\n", nrow(LD_full), ncol(LD_full)))

# Subset if needed
if (length(common) < length(pipe_snps)) {
  keep <- pipe_snps %in% common
  LD <- LD_full[keep, keep]
  cat(sprintf("  Subsetted to: %d x %d\n", nrow(LD), ncol(LD)))
} else {
  LD <- LD_full
  cat("  All SNPs common — no subsetting needed\n")
}
diag(LD) <- 1.0

# Quick health check
ev <- eigen(LD, only.values = TRUE)$values
cat(sprintf("  Min eigenvalue: %.4f | Negative count: %d\n",
    min(ev), sum(ev < 0)))

# ── 3. Compare LD orientation against z-scores (sign check) ──────────────────
cat("\n=== 3. LD SIGN CONSISTENCY CHECK ===\n")
# For each SNP, predicted direction of correlation with lead SNP
# should match the direction implied by z-scores
lead_idx <- which.max(abs(z_sub$zscore))
lead_snp <- z_sub$rsid[lead_idx]
cat(sprintf("  Lead SNP: %s (z=%.2f)\n", lead_snp, z_sub$zscore[lead_idx]))

ld_with_lead <- LD[lead_idx, ]
z_corr <- cor(z_sub$zscore, ld_with_lead)
cat(sprintf("  Pearson r(z-score vs LD-with-lead): %.4f\n", z_corr))
cat(sprintf("  [Expect r > 0.3 for a well-oriented LD matrix]\n"))

# ── 4. Write FINEMAP inputs ───────────────────────────────────────────────────
cat("\n=== 4. WRITING FINEMAP INPUTS ===\n")

# z file
fwrite(z_sub[, .(rsid, chromosome, position, allele1, allele2, maf, beta, se)],
       file.path(OUTDIR, "norcloz_chr2_pipe.z"),
       sep = " ", quote = FALSE, eol = "\n")

# LD file (exact 1s on diagonal)
fwrite(as.data.table(LD), file.path(OUTDIR, "norcloz_chr2_pipe.ld"),
       sep = " ", col.names = FALSE, quote = FALSE, eol = "\n")

# infile (Unix line endings for WSL)
infile_lines <- c(
  "z;ld;snp;config;cred;log;n_samples",
  paste(
    paste0(WSL_OUTDIR, "/norcloz_chr2_pipe.z"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_pipe.ld"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_pipe_k5.snp"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_pipe_k5.config"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_pipe_k5.cred"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_pipe_k5.log"),
    N, sep = ";"
  )
)
con <- file(file.path(OUTDIR, "norcloz_chr2_pipe.infile"), open = "wb")
writeLines(infile_lines, con, sep = "\n")
close(con)
cat(sprintf("  Written to: %s\n", OUTDIR))
cat(sprintf("  SNPs: %d | N: %d | k: %d\n", nrow(z_sub), N, K))

# ── 5. Run FINEMAP ────────────────────────────────────────────────────────────
cat("\n=== 5. RUNNING FINEMAP (--sss --n-causal-snps 5) ===\n")
infile_wsl <- paste0(WSL_OUTDIR, "/norcloz_chr2_pipe.infile")
cmd <- sprintf("wsl %s --sss --in-files %s --n-causal-snps %d --log",
               FINEMAP_BIN, infile_wsl, K)
cat("  Command:", cmd, "\n")
t0 <- proc.time()
rc <- system(cmd, intern = FALSE)
cat(sprintf("  Exit code: %d | Time: %.1fs\n", rc, (proc.time() - t0)["elapsed"]))

# ── 6. Parse and report results ───────────────────────────────────────────────
cat("\n=== 6. RESULTS ===\n")

log_file <- file.path(OUTDIR, "norcloz_chr2_pipe_k5.log")
snp_file <- file.path(OUTDIR, "norcloz_chr2_pipe_k5.snp")
cred_file <- file.path(OUTDIR, "norcloz_chr2_pipe_k5.cred")

if (!file.exists(log_file)) {
  cat("  ERROR: FINEMAP log not found — run may have failed\n")
  quit(status = 1)
}

# Parse log for Pr(k) posteriors
log_lines <- readLines(log_file)
cat("\n  --- FINEMAP posterior Pr(k) ---\n")
in_section <- FALSE
for (ln in log_lines) {
  if (grepl("Post-Pr\\(# of causal", ln)) in_section <- TRUE
  if (in_section) {
    cat(" ", ln, "\n")
    if (grepl("^\\s*5 ->", ln)) break
  }
}

# Expected k=5 probability (the old bug)
k5_line <- grep("^\\s*5 ->", log_lines, value = TRUE)
if (length(k5_line) > 0) {
  k5_prob <- as.numeric(trimws(gsub(".*-> ", "", k5_line)))
  if (!is.na(k5_prob) && k5_prob > 0.5) {
    cat(sprintf("\n  FAIL: Pr(k=5) = %.3f — this is the old bug!\n", k5_prob))
  } else {
    cat(sprintf("\n  PASS: Pr(k=5) = %.3f — not saturated\n", k5_prob))
  }
}

# Parse credible sets
if (file.exists(snp_file)) {
  snps <- fread(snp_file)
  cat(sprintf("\n  Top 5 SNPs by PIP:\n"))
  print(head(snps[order(-prob), .(rsid, prob, log10bf)], 5))
}

# Count credible sets from cred files
cred_files <- list.files(OUTDIR, pattern = "norcloz_chr2_pipe_k5\\.cred[0-9]+$", full.names = TRUE)
cat(sprintf("\n  Credible set files found: %d\n", length(cred_files)))
for (cf in sort(cred_files)) {
  lines <- readLines(cf)
  # First header line has Pr(k)
  pr_line <- lines[grep("^# Post-Pr", lines)[1]]
  cat(sprintf("  %s: %s\n", basename(cf), pr_line))
  # Count SNPs (non-header, non-NA rows for this CS)
  dat <- fread(cf, skip = grep("^index", lines) - 1)
  cs_num <- as.integer(gsub(".*cred([0-9]+)$", "\\1", basename(cf)))
  col <- paste0("cred", cs_num)
  if (col %in% names(dat)) {
    n_snps <- sum(!is.na(dat[[col]]))
    top_snp <- dat[[col]][!is.na(dat[[col]])][1]
    cat(sprintf("    -> %d SNPs | top: %s\n", n_snps, top_snp))
  }
}

# Compare to tutorial
cat("\n  --- Comparison to tutorial FINEMAP (results/) ---\n")
exp_log <- file.path(RESULTDIR, "norcloz_chr2_k5.log_sss")
if (file.exists(exp_log)) {
  exp_lines <- readLines(exp_log)
  cat("  Tutorial Pr(k):\n")
  in_sec <- FALSE
  for (ln in exp_lines) {
    if (grepl("Post-Pr\\(# of causal", ln)) in_sec <- TRUE
    if (in_sec) { cat(" ", ln, "\n"); if (grepl("5 ->", ln)) break }
  }
}

cat("\n=== VERDICT ===\n")
if (file.exists(snp_file) && length(cred_files) >= 2) {
  k5 <- if (length(k5_line) > 0) as.numeric(trimws(gsub(".*-> ", "", k5_line))) else NA
  if (!is.na(k5) && k5 < 0.1) {
    cat("PASS: Pipeline LD correctly resolves to multiple credible sets; Pr(k=5) is low\n")
  } else {
    cat("CHECK: Review Pr(k) and credible set structure above\n")
  }
} else {
  cat("INCOMPLETE: FINEMAP did not produce expected output files\n")
}

cat("\n=== DONE ===\n")
