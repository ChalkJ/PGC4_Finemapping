# ============================================================================
# Phase 1: Test fine-mapping pipeline with tutorial norcloz_chr2 data
# Uses pre-computed .z and .ld files — tests SuSiE + FINEMAP end-to-end
# Expected: 2 well-resolved credible sets, k=5
# ============================================================================

library(data.table)
library(susieR)

TESTDIR    <- "C:/Users/c1977426/OneDrive - Cardiff University/Documents/sczvsbp/met588_finemapping"
RESULTDIR  <- file.path(TESTDIR, "results")
OUTDIR     <- file.path(TESTDIR, "test_phase1_output")
WSL_OUTDIR <- "/home/c1977426/sczvsbp/met588_finemapping/test_phase1_output"
FINEMAP_BIN <- "/home/c1977426/bin/finemap"

dir.create(OUTDIR, showWarnings = FALSE)

N <- 2989
L <- 5

# ── 1. Load sumstats ─────────────────────────────────────────────────────────
cat("=== 1. LOADING DATA ===\n")
z_dt <- fread(file.path(RESULTDIR, "norcloz_chr2.z"))
z_dt[, zscore := beta / se]
cat(sprintf("  %d SNPs | chr2:%d-%d\n",
    nrow(z_dt), min(z_dt$position), max(z_dt$position)))
cat(sprintf("  z-score range: %.2f to %.2f | max |z|: %.2f (%s)\n",
    min(z_dt$zscore), max(z_dt$zscore),
    max(abs(z_dt$zscore)), z_dt$rsid[which.max(abs(z_dt$zscore))]))

# ── 2. Load LD ───────────────────────────────────────────────────────────────
cat("\n=== 2. LOADING LD MATRIX ===\n")
cat(sprintf("  File size: %.1f MB\n",
    file.size(file.path(RESULTDIR, "norcloz_chr2.ld")) / 1e6))
LD <- as.matrix(fread(file.path(RESULTDIR, "norcloz_chr2.ld"), header = FALSE))
stopifnot(nrow(LD) == nrow(z_dt))
rownames(LD) <- colnames(LD) <- z_dt$rsid
diag(LD) <- 1.0
cat(sprintf("  %d x %d matrix\n", nrow(LD), ncol(LD)))

ev <- eigen(LD, only.values = TRUE)$values
cat(sprintf("  Min eigenvalue: %.4f | Negative count: %d\n",
    min(ev), sum(ev < 0)))

# ── 3. SuSiE-RSS ─────────────────────────────────────────────────────────────
cat("\n=== 3. SuSiE-RSS (L=5, lambda=0.1) ===\n")
t0 <- proc.time()
fit <- tryCatch(
  susie_rss(z       = z_dt$zscore,
            R       = LD,
            n       = N,
            L       = L,
            lambda  = 0.1,
            estimate_residual_variance = TRUE),
  error = function(e) { cat("  SuSiE ERROR:", conditionMessage(e), "\n"); NULL }
)
cat(sprintf("  Time: %.1fs\n", (proc.time() - t0)["elapsed"]))

susie_ok <- !is.null(fit)
if (susie_ok) {
  n_cs <- length(fit$sets$cs)
  cat(sprintf("  Credible sets: %d\n", n_cs))
  cat(sprintf("  Max PIP: %.4f (%s)\n",
      max(fit$pip), z_dt$rsid[which.max(fit$pip)]))
  for (i in seq_along(fit$sets$cs)) {
    idx  <- fit$sets$cs[[i]]
    best <- idx[which.max(fit$pip[idx])]
    cat(sprintf("  CS%d: %d SNPs | purity=%.3f | top: %s (PIP=%.4f, pos=%d)\n",
        i, length(idx),
        fit$sets$purity$min.abs.corr[i],
        z_dt$rsid[best], fit$pip[best], z_dt$position[best]))
  }
  saveRDS(fit, file.path(OUTDIR, "norcloz_chr2.susie.rds"))
} else {
  cat("  SuSiE returned NULL — check errors above\n")
}

# ── 4. FINEMAP input prep ─────────────────────────────────────────────────────
cat("\n=== 4. FINEMAP INPUT PREP ===\n")
fwrite(z_dt[, .(rsid, chromosome, position, allele1, allele2, maf, beta, se)],
       file.path(OUTDIR, "norcloz_chr2.z"), sep = " ", quote = FALSE, eol = "\n")

LD_fm       <- LD
diag(LD_fm) <- 1.0
fwrite(as.data.table(LD_fm), file.path(OUTDIR, "norcloz_chr2.ld"),
       sep = " ", col.names = FALSE, quote = FALSE, eol = "\n")

# FINEMAP runs in Linux (WSL) so the infile must use Unix line endings (\n not \r\n)
infile_lines <- c(
  "z;ld;snp;config;cred;log;n_samples",
  paste(
    paste0(WSL_OUTDIR, "/norcloz_chr2.z"),
    paste0(WSL_OUTDIR, "/norcloz_chr2.ld"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_k5.snp"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_k5.config"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_k5.cred"),
    paste0(WSL_OUTDIR, "/norcloz_chr2_k5.log"),
    N, sep = ";"
  )
)
con <- file(file.path(OUTDIR, "norcloz_chr2.infile"), open = "wb")
writeLines(infile_lines, con, sep = "\n")
close(con)
cat(sprintf("  Written to %s\n", OUTDIR))

# ── 5. Run FINEMAP ────────────────────────────────────────────────────────────
cat("\n=== 5. RUNNING FINEMAP (--sss --n-causal-snps 5) ===\n")
infile_wsl <- paste0(WSL_OUTDIR, "/norcloz_chr2.infile")
cmd <- sprintf('wsl %s --sss --in-files %s --n-causal-snps %d --log',
               FINEMAP_BIN, infile_wsl, L)
cat("  Command:", cmd, "\n")
t0 <- proc.time()
rc <- system(cmd, intern = FALSE)
elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("  Exit code: %d | Time: %.1fs\n", rc, elapsed))

# ── 6. Compare against expected output ───────────────────────────────────────
cat("\n=== 6. COMPARISON vs TUTORIAL EXPECTED OUTPUT ===\n")

exp_snp  <- fread(file.path(RESULTDIR, "norcloz_chr2_k5.snp"))
exp_cred <- fread(file.path(RESULTDIR, "norcloz_chr2_k5.cred"))

cat(sprintf("\nExpected (tutorial): %d SNPs, top PIP=%.4f (%s)\n",
    nrow(exp_snp), exp_snp$prob[1], exp_snp$rsid[1]))
cat(sprintf("Expected credible set: %d SNPs at 95%% confidence\n", nrow(exp_cred)))

# SuSiE comparison
if (susie_ok) {
  pip_dt <- data.table(rsid = z_dt$rsid, susie_pip = fit$pip)
  top_sus <- head(pip_dt[order(-susie_pip)], 10)
  top_exp <- head(exp_snp[order(-prob)][, .(rsid, finemap_prob = prob, log10bf)], 10)
  merged10 <- merge(top_sus, top_exp, by = "rsid", all = TRUE)
  cat("\nTop 10 SuSiE PIPs vs expected FINEMAP:\n")
  print(merged10)

  all_merged <- merge(pip_dt, exp_snp[, .(rsid, finemap_prob = prob)], by = "rsid")
  cat(sprintf("\nSuSiE vs FINEMAP PIP correlation (all %d overlapping SNPs): r=%.4f\n",
      nrow(all_merged), cor(all_merged$susie_pip, all_merged$finemap_prob)))
}

# Our FINEMAP vs expected FINEMAP
our_snp_file  <- file.path(OUTDIR, "norcloz_chr2_k5.snp")
our_cred_file <- file.path(OUTDIR, "norcloz_chr2_k5.cred")
if (file.exists(our_snp_file)) {
  our_snp  <- fread(our_snp_file)
  our_cred <- fread(our_cred_file)

  cat(sprintf("\nOur FINEMAP: %d SNPs, top PIP=%.4f (%s)\n",
      nrow(our_snp), our_snp$prob[1], our_snp$rsid[1]))
  cat(sprintf("Our credible set: %d SNPs at 95%% confidence\n", nrow(our_cred)))

  both <- merge(exp_snp[, .(rsid, exp_prob = prob)],
                our_snp[, .(rsid, our_prob = prob)], by = "rsid")
  cat(sprintf("Our vs Expected FINEMAP PIP correlation (n=%d): r=%.6f\n",
      nrow(both), cor(both$exp_prob, both$our_prob)))
  cat(sprintf("Max |PIP diff|: %.6f\n", max(abs(both$exp_prob - both$our_prob))))

  top_our <- head(our_snp[order(-prob)][, .(rsid, our_prob = prob, log10bf)], 5)
  top_exp2 <- head(exp_snp[order(-prob)][, .(rsid, exp_prob = prob)], 5)
  cat("\nTop 5 SNPs comparison:\n")
  print(merge(top_our, top_exp2, by = "rsid", all = TRUE))
} else {
  cat("\nFINEMAP .snp output not found — FINEMAP may have failed\n")
  cat("  Check FINEMAP log:", file.path(OUTDIR, "norcloz_chr2_k5.log"), "\n")
}

cat("\n=== DONE ===\n")
