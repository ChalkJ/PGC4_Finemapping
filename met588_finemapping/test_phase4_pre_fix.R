# ============================================================================
# Phase 4: Run FINEMAP k=5 on PRE-FIX (unoriented allele) LD
# Goal: reproduce the old bug — expect Pr(k=5)≈1 and degenerate CS2..CS5
# Contrast with Phase 3 (post-fix) results.
# ============================================================================

library(data.table)

TESTDIR     <- "C:/Users/c1977426/OneDrive - Cardiff University/Documents/sczvsbp/met588_finemapping"
RESULTDIR   <- file.path(TESTDIR, "results")
WORKDIR     <- file.path(TESTDIR, "test_phase2_work_pre_fix")
OUTDIR      <- file.path(TESTDIR, "test_phase4_output")
WSL_OUTDIR  <- "/home/c1977426/sczvsbp/met588_finemapping/test_phase4_output"
FINEMAP_BIN <- "/home/c1977426/bin/finemap"
LDMERGE     <- file.path(TESTDIR, "test_phase2_work", "LDmerge_patched.R")
RSCRIPT_WIN <- "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"
LOCUS       <- "norcloz_chr2_test"
N           <- 2989
K           <- 5

dir.create(WORKDIR, showWarnings = FALSE)
dir.create(OUTDIR,  showWarnings = FALSE)

# ── 1. Run pre-fix PLINK pipeline in WSL ─────────────────────────────────────
cat("=== 1. RUNNING PRE-FIX PLINK PIPELINE (WSL) ===\n")
sh_wsl <- "/home/c1977426/sczvsbp/met588_finemapping/test_phase2_setup_pre_fix.sh"
cmd <- sprintf("wsl bash %s", sh_wsl)
cat("  Command:", cmd, "\n")
t0 <- proc.time()
rc <- system(cmd, intern = FALSE)
cat(sprintf("  Exit code: %d | Time: %.1fs\n", rc, (proc.time() - t0)["elapsed"]))

if (rc != 0) {
  cat("  PLINK step failed — aborting\n")
  quit(status = 1)
}

# ── 2. Run LDmerge from pre-fix work dir (Windows Rscript; R not in WSL) ──────
cat("\n=== 2. RUNNING LDmerge (single-cohort, pre-fix LD) ===\n")

old_wd <- getwd()
setwd(WORKDIR)
ldmerge_cmd <- sprintf('"%s" --vanilla "%s" %s PLINK METAL', RSCRIPT_WIN, LDMERGE, LOCUS)
cat("  Command:", ldmerge_cmd, "\n")
t0 <- proc.time()
rc2 <- system(ldmerge_cmd, intern = FALSE)
setwd(old_wd)
cat(sprintf("  Exit code: %d | Time: %.1fs\n", rc2, (proc.time() - t0)["elapsed"]))

if (rc2 != 0) {
  cat("  LDmerge failed — check output above\n")
  quit(status = 1)
}

# ── 3. Load z-file and pre-fix LD ─────────────────────────────────────────────
cat("\n=== 3. LOADING DATA ===\n")
z_dt <- fread(file.path(RESULTDIR, "norcloz_chr2.z"))
z_dt[, zscore := beta / se]

pipe_snps <- readLines(file.path(WORKDIR, paste0(LOCUS, ".snp.log")))
cat(sprintf("  z-file: %d SNPs | Pre-fix LD SNPs: %d\n", nrow(z_dt), length(pipe_snps)))

common <- intersect(z_dt$rsid, pipe_snps)
cat(sprintf("  Common SNPs: %d\n", length(common)))

z_sub <- z_dt[rsid %in% common]
z_sub <- z_sub[match(pipe_snps[pipe_snps %in% common], rsid)]

ld_file <- file.path(WORKDIR, paste0(LOCUS, ".ld"))
cat(sprintf("  LD file size: %.1f MB\n", file.size(ld_file) / 1e6))
LD_full <- as.matrix(fread(ld_file, header = FALSE))
rownames(LD_full) <- colnames(LD_full) <- pipe_snps

if (length(common) < length(pipe_snps)) {
  keep <- pipe_snps %in% common
  LD <- LD_full[keep, keep]
  cat(sprintf("  Subsetted matrix: %d x %d\n", nrow(LD), ncol(LD)))
} else {
  LD <- LD_full
  cat(sprintf("  Full matrix: %d x %d\n", nrow(LD), ncol(LD)))
}
diag(LD) <- 1.0

# Sign consistency check — pre-fix should show ~0 correlation
lead_idx <- which.max(abs(z_sub$zscore))
lead_snp <- z_sub$rsid[lead_idx]
ld_with_lead <- LD[lead_idx, ]
z_corr <- cor(z_sub$zscore, ld_with_lead)
cat(sprintf("\n  Lead SNP: %s (z=%.2f)\n", lead_snp, z_sub$zscore[lead_idx]))
cat(sprintf("  Pearson r(z-score vs LD-with-lead): %.4f\n", z_corr))
cat(sprintf("  [Pre-fix: expect r ≈ 0; post-fix: r > 0.3]\n"))

# ── 4. Write FINEMAP inputs ───────────────────────────────────────────────────
cat("\n=== 4. WRITING FINEMAP INPUTS ===\n")

fwrite(z_sub[, .(rsid, chromosome, position, allele1, allele2, maf, beta, se)],
       file.path(OUTDIR, "norcloz_chr2_prefix.z"),
       sep = " ", quote = FALSE, eol = "\n")

fwrite(as.data.table(LD),
       file.path(OUTDIR, "norcloz_chr2_prefix.ld"),
       sep = " ", col.names = FALSE, quote = FALSE, eol = "\n")

infile_lines <- c(
  "z;ld;snp;config;cred;log;n_samples",
  paste(
    paste0(WSL_OUTDIR, "/norcloz_chr2_prefix.z"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_prefix.ld"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_prefix_k5.snp"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_prefix_k5.config"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_prefix_k5.cred"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_prefix_k5.log"),
    N, sep = ";"
  )
)
con <- file(file.path(OUTDIR, "norcloz_chr2_prefix.infile"), open = "wb")
writeLines(infile_lines, con, sep = "\n")
close(con)
cat(sprintf("  Written to: %s\n", OUTDIR))
cat(sprintf("  SNPs: %d | N: %d | k: %d\n", nrow(z_sub), N, K))

# ── 5. Run FINEMAP ────────────────────────────────────────────────────────────
cat("\n=== 5. RUNNING FINEMAP (--sss --n-causal-snps 5) ===\n")
infile_wsl <- paste0(WSL_OUTDIR, "/norcloz_chr2_prefix.infile")
cmd <- sprintf("wsl %s --sss --in-files %s --n-causal-snps %d --log",
               FINEMAP_BIN, infile_wsl, K)
cat("  Command:", cmd, "\n")
t0 <- proc.time()
rc3 <- system(cmd, intern = FALSE)
cat(sprintf("  Exit code: %d | Time: %.1fs\n", rc3, (proc.time() - t0)["elapsed"]))

# ── 6. Parse and report results ───────────────────────────────────────────────
cat("\n=== 6. RESULTS ===\n")

# FINEMAP writes .log_sss, not .log
log_file_sss <- file.path(OUTDIR, "norcloz_chr2_prefix_k5.log_sss")
snp_file     <- file.path(OUTDIR, "norcloz_chr2_prefix_k5.snp")

# Try .log_sss first, fall back to .log
if (file.exists(log_file_sss)) {
  log_lines <- readLines(log_file_sss)
  cat("  Reading from .log_sss\n")
} else {
  log_file <- file.path(OUTDIR, "norcloz_chr2_prefix_k5.log")
  if (!file.exists(log_file)) {
    cat("  ERROR: No FINEMAP log found — run may have failed\n")
    quit(status = 1)
  }
  log_lines <- readLines(log_file)
  cat("  Reading from .log\n")
}

cat("\n  --- Pre-fix FINEMAP posterior Pr(k) ---\n")
in_section <- FALSE
for (ln in log_lines) {
  if (grepl("Post-Pr\\(# of causal", ln)) in_section <- TRUE
  if (in_section) {
    cat(" ", ln, "\n")
    if (grepl("^\\s*5 ->", ln)) break
  }
}
# Expected k
exp_line <- grep("Post-expected # of causal", log_lines, value = TRUE)
if (length(exp_line) > 0) cat(" ", exp_line[1], "\n")

# Extract key values
k5_line   <- grep("^\\s*5 ->", log_lines, value = TRUE)
k2_line   <- grep("^\\s*2 ->", log_lines, value = TRUE)
k5_prob   <- if (length(k5_line) > 0) as.numeric(trimws(gsub(".*-> ", "", k5_line[1]))) else NA
k2_prob   <- if (length(k2_line) > 0) as.numeric(trimws(gsub(".*-> ", "", k2_line[1]))) else NA
exp_match <- regmatches(exp_line, regexpr("[0-9\\.]+$", exp_line))
exp_k     <- if (length(exp_match) > 0) as.numeric(exp_match[1]) else NA

# Parse credible sets
cred_files <- list.files(OUTDIR, pattern = "norcloz_chr2_prefix_k5\\.cred[0-9]+$", full.names = TRUE)
cat(sprintf("\n  Credible set files: %d\n", length(cred_files)))

cs_top <- character(0)
cs_logbf <- numeric(0)
cs_purity <- numeric(0)

for (cf in sort(cred_files)) {
  lines <- readLines(cf)
  pr_line  <- lines[grep("^# Post-Pr", lines)[1]]
  bf_line  <- lines[grep("^#log10bf", lines)[1]]
  pur_line <- lines[grep("^#min\\(\\|ld\\|\\)", lines)[1]]
  cs_num   <- as.integer(gsub(".*cred([0-9]+)$", "\\1", basename(cf)))

  # Extract log10bf for this CS
  bf_vals  <- as.numeric(strsplit(gsub("^#log10bf\\s+", "", bf_line), "\\s+")[[1]])
  pur_vals <- as.numeric(strsplit(gsub("^#min\\(\\|ld\\|\\)\\s+", "", pur_line), "\\s+")[[1]])

  skip_rows <- grep("^index", lines) - 1
  if (length(skip_rows) > 0) {
    dat   <- fread(cf, skip = skip_rows)
    col   <- paste0("cred", cs_num)
    if (col %in% names(dat)) {
      valid <- dat[[col]][!is.na(dat[[col]])]
      n_snp <- length(valid)
      top   <- if (n_snp > 0) valid[1] else NA_character_
    } else { n_snp <- NA; top <- NA_character_ }
  } else { n_snp <- NA; top <- NA_character_ }

  cs_top   <- c(cs_top, top)
  cs_logbf <- c(cs_logbf, if (!is.na(bf_vals[1])) bf_vals[1] else NA)
  cs_purity <- c(cs_purity, if (!is.na(pur_vals[1])) pur_vals[1] else NA)

  cat(sprintf("  CS%d: top=%s | n=%s | log10BF=%.3f | min|ld|=%.3f\n",
      cs_num, top, n_snp,
      if (!is.na(bf_vals[1])) bf_vals[1] else NA,
      if (!is.na(pur_vals[1])) pur_vals[1] else NA))
}

# ── 7. Three-column comparison table ─────────────────────────────────────────
cat("\n")
cat("=====================================================================\n")
cat(" COMPARISON TABLE: Tutorial LD vs Pipeline LD (post-fix vs pre-fix)\n")
cat("=====================================================================\n")

fmt_val <- function(x, fmt = "%.3f") if (is.na(x)) "N/A" else sprintf(fmt, x)

cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "", "Tutorial LD", "Pipeline (post-fix)", "Pipeline (pre-fix)"))
cat(paste(rep("-", 84), collapse = ""), "\n")
cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "Pr(k=2)", "0.677", "0.639", fmt_val(k2_prob)))
cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "Pr(k=5)", "0", "0", fmt_val(k5_prob)))
cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "Post-expected k", "2.02", "2.07", fmt_val(exp_k, "%.2f")))
cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "CS1 top SNP", "rs2011425", "rs2011425",
    if (length(cs_top) >= 1) cs_top[1] else "N/A"))
cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "CS2 top SNP", "rs13027283", "rs13027283",
    if (length(cs_top) >= 2) cs_top[2] else "N/A"))
cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "CS1 log10BF", "7.116", "7.116",
    if (length(cs_logbf) >= 1) fmt_val(cs_logbf[1]) else "N/A"))
cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "CS2 log10BF", "4.180", "4.180",
    if (length(cs_logbf) >= 2) fmt_val(cs_logbf[2]) else "N/A"))
cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "CS1 purity (min|ld|)", "0.901", "0.901",
    if (length(cs_purity) >= 1) fmt_val(cs_purity[1]) else "N/A"))
cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "CS2 purity (min|ld|)", "0 (independent)", "0 (independent)",
    if (length(cs_purity) >= 2) fmt_val(cs_purity[2]) else "N/A"))
cat(sprintf("%-26s  %-18s  %-18s  %-18s\n",
    "LD sign r (z vs lead)", ">0.3 (expected)", ">0.3 (verified)",
    fmt_val(z_corr)))
cat(paste(rep("-", 84), collapse = ""), "\n")

cat("\n=== VERDICT ===\n")
if (!is.na(k5_prob) && k5_prob > 0.5) {
  cat(sprintf("BUG CONFIRMED: Pr(k=5) = %.3f — pre-fix LD produces the degenerate result\n", k5_prob))
} else if (!is.na(k5_prob) && k5_prob < 0.1) {
  cat(sprintf("NOTE: Pr(k=5) = %.3f — pre-fix LD does NOT reproduce the bug (unexpected)\n", k5_prob))
} else {
  cat(sprintf("Pr(k=5) = %s — check full output above\n", fmt_val(k5_prob)))
}

cat("\n=== DONE ===\n")
