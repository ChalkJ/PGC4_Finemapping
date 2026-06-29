#!/usr/bin/env Rscript
# =============================================================================
# Fine-Mapping Pipeline — SCZ vs BD1
# SuSiE-RSS (L=5) + FINEMAP-SSS (k=5) for 10 BD1 loci
# Inherits all fixes from sczvsbp_v2.r review:
#   - WSL_BASE corrected
#   - lambda = 0.1 for non-PSD LD matrices
#   - fit$snp_ids set for locus_zoom_bd1.Rmd compatibility
#   - FINEMAP Pr(k) reads posterior section only
#   - RDS named {locus}.susie_l5.rds to match locus_zoom pattern
# =============================================================================

library(data.table)
library(susieR)

# ── Paths ────────────────────────────────────────────────────────────────────
DATADIR     <- "C:/PGC4"
LDDIR       <- file.path(DATADIR, "ld/sczvsbp1")
SUSIEDIR    <- file.path(DATADIR, "r_files/sczvsbp1/susie")
FMDIR       <- file.path(DATADIR, "r_files/sczvsbp1/finemap")
SUMSTATS    <- file.path(DATADIR, "sumstats/daner_scz_vs_bip1_1025_noduppos_hetpva.gz")
FINEMAP_BIN <- "/home/c1977426/bin/finemap"

WSL_BASE <- "/home/c1977426/PGC4"

win_to_wsl <- function(path) {
  path <- gsub("\\\\", "/", path)
  rel  <- sub(gsub("\\\\", "/", DATADIR), "", path, fixed = TRUE)
  paste0(WSL_BASE, rel)
}

# ── Configuration ────────────────────────────────────────────────────────────
N_CAUSAL              <- 5
SKIP_EXISTING         <- TRUE
SUSIE_EST_RES_VAR     <- TRUE
K_POSTERIOR_THRESHOLD <- 0.5
N_CAUSAL_RERUN        <- 10
RERUN_MAX_ITER        <- 500

dir.create(SUSIEDIR, showWarnings = FALSE)
dir.create(FMDIR,    showWarnings = FALSE)

# Discover loci from ld_bd1/
snp_log_files <- list.files(LDDIR, pattern = "^[0-9]{3}\\.snp\\.log$", full.names = FALSE)
loci_all      <- sort(sub("\\.snp\\.log$", "", snp_log_files))
cat("Loci found in ld_bd1:", paste(loci_all, collapse = " "), "\n\n")

# =============================================================================
# 1. Load sumstats
# =============================================================================
cat("Loading BD1 sumstats...\n")
sumstats <- fread(SUMSTATS)
sumstats <- sumstats[HetPVa > 0.05]
sumstats[, zscore := log(OR) / SE]

N <- round(mean(4 / (1/sumstats$Nca + 1/sumstats$Nco), na.rm = TRUE))
cat("Effective sample size N =", N, "\n\n")

# =============================================================================
# 2. Per-locus: SuSiE + FINEMAP input prep
# =============================================================================
cat(strrep("=", 60), "\n")
cat("FINE-MAPPING (SuSiE + FINEMAP input prep)\n")
cat(strrep("=", 60), "\n\n")

fm_master_rows  <- list()
loci_processed  <- character()

for (locus in loci_all) {

  susie_out    <- file.path(SUSIEDIR, paste0(locus, ".susie_l", N_CAUSAL, ".rds"))
  fm_locus_dir <- file.path(FMDIR, locus)
  fm_z_file    <- file.path(fm_locus_dir, paste0(locus, ".z"))
  fm_ld_file   <- file.path(fm_locus_dir, paste0(locus, ".ld"))
  fm_snp_file  <- file.path(fm_locus_dir, paste0(locus, ".snp"))

  susie_done   <- SKIP_EXISTING && file.exists(susie_out)
  finemap_done <- SKIP_EXISTING && file.exists(fm_snp_file)

  if (susie_done && finemap_done) {
    cat("[", locus, "] Both complete, skipping\n")
    fm_master_rows[[locus]] <- data.table(
      z = fm_z_file, ld = fm_ld_file,
      snp    = fm_snp_file,
      config = file.path(fm_locus_dir, paste0(locus, ".config")),
      cred   = file.path(fm_locus_dir, paste0(locus, ".cred")),
      log    = file.path(fm_locus_dir, paste0(locus, ".log")),
      n_samples = N
    )
    loci_processed <- c(loci_processed, locus)
    next
  }

  snplog_file <- file.path(LDDIR, paste0(locus, ".snp.log"))
  ldfile      <- file.path(LDDIR, paste0(locus, ".ld.gz"))

  if (!file.exists(snplog_file) || !file.exists(ldfile)) {
    cat("[", locus, "] WARNING: LD files not found — skipping\n")
    next
  }

  ld_snps <- fread(snplog_file, header = FALSE)[[1]]
  LD      <- as.matrix(fread(ldfile, header = FALSE))
  rownames(LD) <- colnames(LD) <- ld_snps
  diag(LD) <- 1.0

  ss_locus    <- sumstats[SNP %in% ld_snps]
  common_snps <- intersect(ld_snps, ss_locus$SNP)

  if (length(common_snps) < 2) {
    cat("[", locus, "] WARNING: <2 overlapping SNPs — skipping\n")
    next
  }

  ss_locus <- ss_locus[match(common_snps, SNP)]
  LD       <- LD[common_snps, common_snps]

  cat("[", locus, "]", length(common_snps), "SNPs |")

  # ── SuSiE ──
  if (!susie_done) {
    z   <- ss_locus$zscore
    fit <- tryCatch(
      susie_rss(
        z = z, R = LD, n = N, L = N_CAUSAL,
        lambda = 0.1,
        estimate_residual_variance = SUSIE_EST_RES_VAR
      ),
      error = function(e) {
        cat(" SuSiE ERROR:", conditionMessage(e), " |")
        return(NULL)
      }
    )
    if (!is.null(fit)) {
      fit$snp      <- common_snps
      fit$snp_ids  <- common_snps
      fit$zscore   <- z
      fit$position <- ss_locus$BP
      fit$chr      <- ss_locus$CHR[1]
      saveRDS(fit, susie_out)
      cat(sprintf(" SuSiE: %d CS, max_pip=%.3f |", length(fit$sets$cs), max(fit$pip)))
    }
  } else {
    cat(" SuSiE: cached |")
  }

  # ── FINEMAP input prep ──
  if (!finemap_done) {
    dir.create(fm_locus_dir, showWarnings = FALSE)

    frq_col <- if ("FRQ_A_41065" %in% colnames(ss_locus)) "FRQ_A_41065" else
               grep("^FRQ_A_", colnames(ss_locus), value = TRUE)[1]

    z_dt <- data.table(
      rsid       = ss_locus$SNP,
      chromosome = ss_locus$CHR,
      position   = ss_locus$BP,
      allele1    = ss_locus$A1,
      allele2    = ss_locus$A2,
      maf        = pmin(ss_locus[[frq_col]], 1 - ss_locus[[frq_col]]),
      beta       = log(ss_locus$OR),
      se         = ss_locus$SE
    )
    fwrite(z_dt, fm_z_file, sep = " ", quote = FALSE, eol = "\n")

    LD_fm <- LD; diag(LD_fm) <- 1.0
    fwrite(as.data.table(LD_fm), fm_ld_file, sep = " ",
           col.names = FALSE, quote = FALSE, eol = "\n")

    cat(" FINEMAP inputs written |")
  } else {
    cat(" FINEMAP: cached |")
  }

  fm_master_rows[[locus]] <- data.table(
    z = fm_z_file, ld = fm_ld_file,
    snp    = fm_snp_file,
    config = file.path(fm_locus_dir, paste0(locus, ".config")),
    cred   = file.path(fm_locus_dir, paste0(locus, ".cred")),
    log    = file.path(fm_locus_dir, paste0(locus, ".log")),
    n_samples = N
  )
  loci_processed <- c(loci_processed, locus)
  cat("\n")
}

# ── Write FINEMAP master file ──
master_dt   <- rbindlist(fm_master_rows)
master_path <- file.path(FMDIR, "master")

master_wsl <- copy(master_dt)
for (col in c("z", "ld", "snp", "config", "cred", "log")) {
  master_wsl[[col]] <- win_to_wsl(master_wsl[[col]])
}
fwrite(master_wsl, master_path, sep = ";", quote = FALSE, eol = "\n")
cat("\nFINEMAP master written:", master_path, "(", nrow(master_dt), "loci)\n")

# =============================================================================
# 3. Run FINEMAP
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("RUNNING FINEMAP --sss (k =", N_CAUSAL, ")\n")
cat(strrep("=", 60), "\n\n")

finemap_cmd <- sprintf(
  "wsl %s --sss --in-files %s --n-causal-snps %d --log",
  FINEMAP_BIN, win_to_wsl(master_path), N_CAUSAL
)
cat("Command:", finemap_cmd, "\n\n")
exit_code <- system(finemap_cmd)
if (exit_code != 0) {
  cat("WARNING: FINEMAP exited with code", exit_code, "\n\n")
} else {
  cat("FINEMAP completed.\n\n")
}

# =============================================================================
# 3b. Adaptive re-run
# =============================================================================
load_locus_data <- function(loc) {
  snplog_file <- file.path(LDDIR, paste0(loc, ".snp.log"))
  ldfile      <- file.path(LDDIR, paste0(loc, ".ld.gz"))
  if (!file.exists(snplog_file) || !file.exists(ldfile)) return(NULL)

  ld_snps <- fread(snplog_file, header = FALSE)[[1]]
  LD      <- as.matrix(fread(ldfile, header = FALSE))
  rownames(LD) <- colnames(LD) <- ld_snps
  diag(LD) <- 1.0

  ss_locus    <- sumstats[SNP %in% ld_snps]
  common_snps <- intersect(ld_snps, ss_locus$SNP)
  if (length(common_snps) < 2) return(NULL)

  ss_locus <- ss_locus[match(common_snps, SNP)]
  LD       <- LD[common_snps, common_snps]
  list(ss_locus = ss_locus, LD = LD, common_snps = common_snps)
}

run_susie_locus <- function(loc, dat, L, max_iter = 100) {
  z   <- dat$ss_locus$zscore
  fit <- tryCatch(
    susie_rss(
      z = z, R = dat$LD, n = N, L = L,
      lambda = 0.1,
      estimate_residual_variance = SUSIE_EST_RES_VAR,
      max_iter = max_iter
    ),
    error = function(e) { cat(" SuSiE ERROR:", conditionMessage(e)); NULL }
  )
  if (!is.null(fit)) {
    fit$snp      <- dat$common_snps
    fit$snp_ids  <- dat$common_snps
    fit$zscore   <- z
    fit$position <- dat$ss_locus$BP
    fit$chr      <- dat$ss_locus$CHR[1]
    saveRDS(fit, file.path(SUSIEDIR, paste0(loc, ".susie_l", L, ".rds")))
  }
  fit
}

run_finemap_loci <- function(loci_to_run, k) {
  fm_rows <- lapply(loci_to_run, function(loc) {
    fm_locus_dir <- file.path(FMDIR, loc)
    data.table(
      z      = file.path(fm_locus_dir, paste0(loc, ".z")),
      ld     = file.path(fm_locus_dir, paste0(loc, ".ld")),
      snp    = file.path(fm_locus_dir, paste0(loc, ".snp")),
      config = file.path(fm_locus_dir, paste0(loc, ".config")),
      cred   = file.path(fm_locus_dir, paste0(loc, ".cred")),
      log    = file.path(fm_locus_dir, paste0(loc, ".log")),
      n_samples = N
    )
  })
  rerun_dt   <- rbindlist(fm_rows)
  rerun_path <- file.path(FMDIR, paste0("master_rerun_k", k))
  rerun_wsl  <- copy(rerun_dt)
  for (col in c("z", "ld", "snp", "config", "cred", "log")) {
    rerun_wsl[[col]] <- win_to_wsl(rerun_wsl[[col]])
  }
  fwrite(rerun_wsl, rerun_path, sep = ";", quote = FALSE, eol = "\n")
  cmd <- sprintf("wsl %s --sss --in-files %s --n-causal-snps %d --log",
                 FINEMAP_BIN, win_to_wsl(rerun_path), k)
  cat("  FINEMAP rerun:", cmd, "\n")
  system(cmd)
}

cat(strrep("=", 60), "\n")
cat("ADAPTIVE RE-RUN: PASS 1\n")
cat(strrep("=", 60), "\n\n")

pass1_loci <- character()
for (loc in loci_processed) {
  needs_rerun <- FALSE; reasons <- character()

  susie_file <- file.path(SUSIEDIR, paste0(loc, ".susie_l", N_CAUSAL, ".rds"))
  if (file.exists(susie_file)) {
    fit <- readRDS(susie_file)
    if (!fit$converged && length(fit$sets$cs) >= N_CAUSAL) {
      needs_rerun <- TRUE
      reasons <- c(reasons, sprintf("SuSiE: %d CS, not converged", length(fit$sets$cs)))
    }
  }

  log_file <- file.path(FMDIR, loc, paste0(loc, ".log_sss"))
  if (file.exists(log_file)) {
    log_lines  <- readLines(log_file)
    post_start <- grep("Post-Pr\\(# of causal", log_lines)
    post_lines <- if (length(post_start) > 0) log_lines[post_start[1]:length(log_lines)] else character(0)
    k_pattern  <- sprintf("^\\s+%d -> ", N_CAUSAL)
    k_line     <- grep(k_pattern, post_lines, value = TRUE)
    if (length(k_line) > 0) {
      k_prob <- as.numeric(sub(".*-> ", "", trimws(k_line[1])))
      if (!is.na(k_prob) && k_prob > K_POSTERIOR_THRESHOLD) {
        needs_rerun <- TRUE
        reasons <- c(reasons, sprintf("FINEMAP: Pr(k=%d) = %.3g", N_CAUSAL, k_prob))
      }
    }
  }

  if (needs_rerun) {
    cat(sprintf("[%s] SATURATED — %s\n", loc, paste(reasons, collapse = "; ")))
    pass1_loci <- c(pass1_loci, loc)
  }
}

if (length(pass1_loci) > 0) {
  cat(sprintf("\nRe-running %d loci with L/k = %d\n", length(pass1_loci), N_CAUSAL_RERUN))
  for (loc in pass1_loci) {
    dat <- load_locus_data(loc)
    if (is.null(dat)) next
    cat(sprintf("[%s] %d SNPs | L: %d → %d |", loc, length(dat$common_snps), N_CAUSAL, N_CAUSAL_RERUN))
    fit <- run_susie_locus(loc, dat, L = N_CAUSAL_RERUN)
    if (!is.null(fit)) cat(sprintf(" SuSiE: %d CS, conv=%s", length(fit$sets$cs), fit$converged))
    cat("\n")
  }
  run_finemap_loci(pass1_loci, k = N_CAUSAL_RERUN)
} else {
  cat("No saturated loci.\n\n")
}

cat(strrep("=", 60), "\n")
cat("ADAPTIVE RE-RUN: PASS 2\n")
cat(strrep("=", 60), "\n\n")

needs_more_iter <- character()
poorly_resolved <- character()

for (loc in pass1_loci) {
  susie_file <- file.path(SUSIEDIR, paste0(loc, ".susie_l", N_CAUSAL_RERUN, ".rds"))
  if (!file.exists(susie_file)) next
  fit <- readRDS(susie_file)
  if (fit$converged) {
    cat(sprintf("[%s] Converged at L=%d with %d CS\n", loc, N_CAUSAL_RERUN, length(fit$sets$cs)))
  } else if (length(fit$sets$cs) >= N_CAUSAL_RERUN) {
    cat(sprintf("[%s] HIT CEILING: %d CS at L=%d — poorly resolved\n",
                loc, length(fit$sets$cs), N_CAUSAL_RERUN))
    poorly_resolved <- c(poorly_resolved, loc)
  } else {
    cat(sprintf("[%s] NEEDS ITERATIONS: re-running with max_iter=%d\n", loc, RERUN_MAX_ITER))
    needs_more_iter <- c(needs_more_iter, loc)
  }
}

if (length(needs_more_iter) > 0) {
  for (loc in needs_more_iter) {
    dat <- load_locus_data(loc)
    if (is.null(dat)) next
    fit <- run_susie_locus(loc, dat, L = N_CAUSAL_RERUN, max_iter = RERUN_MAX_ITER)
    if (!is.null(fit)) cat(sprintf("[%s] %d CS, conv=%s\n", loc, length(fit$sets$cs), fit$converged))
  }
  run_finemap_loci(needs_more_iter, k = N_CAUSAL_RERUN)
}

if (length(poorly_resolved) > 0) {
  cat("\nPOORLY RESOLVED:", paste(poorly_resolved, collapse = ", "), "\n")
  cat("Complex LD — interpret results with caution.\n\n")
}

cat("Adaptive re-run complete.\n\n")

# =============================================================================
# 4. Diagnostics
# =============================================================================
cat(strrep("=", 60), "\n")
cat("DIAGNOSTICS\n")
cat(strrep("=", 60), "\n\n")

# Determine which L was actually used per locus (N_CAUSAL or N_CAUSAL_RERUN)
best_L <- setNames(rep(N_CAUSAL, length(loci_all)), loci_all)
for (loc in pass1_loci) {
  if (file.exists(file.path(SUSIEDIR, paste0(loc, ".susie_l", N_CAUSAL_RERUN, ".rds"))))
    best_L[loc] <- N_CAUSAL_RERUN
}

build_summary <- function() {
  rows <- lapply(loci_all, function(locus) {

    L_used     <- best_L[locus]
    susie_file <- file.path(SUSIEDIR, paste0(locus, ".susie_l", L_used, ".rds"))

    if (file.exists(susie_file)) {
      fit <- readRDS(susie_file)
      cs_lbf <- if (length(fit$sets$cs) > 0)
        sapply(seq_along(fit$sets$cs), function(i)
          round(max(fit$lbf_variable[i, fit$sets$cs[[i]]], na.rm = TRUE), 2))
        else numeric(0)
      cs_size   <- sapply(fit$sets$cs, length)
      cs_purity <- if (!is.null(fit$sets$purity) && nrow(fit$sets$purity) > 0)
        round(fit$sets$purity[, "min.abs.corr"], 3)
        else rep(NA_real_, length(fit$sets$cs))

      susie_dt <- data.table(
        s_L         = L_used,
        s_n_cs      = length(fit$sets$cs),
        s_max_pip   = round(max(fit$pip), 4),
        s_top_snp   = fit$snp[which.max(fit$pip)],
        s_converged = fit$converged,
        s_sigma2    = round(fit$sigma2, 4),
        s_cs_sizes  = paste(cs_size, collapse = ","),
        s_cs_lbf    = paste(cs_lbf, collapse = ","),
        s_cs_purity = paste(cs_purity, collapse = ",")
      )
    } else {
      susie_dt <- data.table(
        s_L = NA_integer_, s_n_cs = NA_integer_, s_max_pip = NA_real_,
        s_top_snp = NA_character_, s_converged = NA,
        s_sigma2 = NA_real_, s_cs_sizes = "", s_cs_lbf = "", s_cs_purity = ""
      )
    }

    fm_locus_dir      <- file.path(FMDIR, locus)
    fm_snp_file       <- file.path(fm_locus_dir, paste0(locus, ".snp"))
    fm_cred_files     <- list.files(fm_locus_dir, pattern = paste0("^", locus, "\\.cred\\d+$"),
                                    full.names = TRUE)

    if (file.exists(fm_snp_file)) {
      fm_snps <- fread(fm_snp_file)

      # Posterior k from .log_sss
      log_file   <- file.path(fm_locus_dir, paste0(locus, ".log_sss"))
      post_exp_k <- NA_real_; post_pr_k5 <- NA_real_
      if (file.exists(log_file)) {
        ll         <- readLines(log_file)
        post_start <- grep("Post-Pr\\(# of causal", ll)
        post_lines <- if (length(post_start) > 0) ll[post_start[1]:length(ll)] else character(0)
        exp_line   <- grep("Post-expected", ll, value = TRUE)
        if (length(exp_line) > 0)
          post_exp_k <- as.numeric(regmatches(exp_line[1], regexpr("[0-9\\.]+$", exp_line[1])))
        k5_line <- grep(sprintf("^\\s+%d -> ", N_CAUSAL), post_lines, value = TRUE)
        if (length(k5_line) > 0)
          post_pr_k5 <- as.numeric(sub(".*-> ", "", trimws(k5_line[1])))
      }

      # Credible set count from cred files
      n_cred <- length(fm_cred_files)

      # Best config
      fm_config_file <- file.path(fm_locus_dir, paste0(locus, ".config"))
      config_dt      <- if (file.exists(fm_config_file)) fread(fm_config_file) else NULL
      best_prob      <- if (!is.null(config_dt) && nrow(config_dt) > 0) config_dt[1, prob] else NA_real_

      fm_dt <- data.table(
        f_n_cred    = n_cred,
        f_max_pip   = round(max(fm_snps$prob, na.rm = TRUE), 4),
        f_top_snp   = fm_snps$rsid[which.max(fm_snps$prob)],
        f_exp_k     = round(post_exp_k, 2),
        f_pr_k5     = round(post_pr_k5, 4),
        f_best_prob = round(best_prob, 4)
      )
    } else {
      fm_dt <- data.table(
        f_n_cred = NA_integer_, f_max_pip = NA_real_, f_top_snp = NA_character_,
        f_exp_k = NA_real_, f_pr_k5 = NA_real_, f_best_prob = NA_real_
      )
    }

    cbind(data.table(locus = locus), susie_dt, fm_dt)
  })
  rbindlist(rows, fill = TRUE)
}

summary_dt <- build_summary()

cat("── SuSiE (L=5) summary ──\n")
print(summary_dt[, .(locus, s_L, s_n_cs, s_max_pip, s_top_snp, s_converged,
                      s_sigma2, s_cs_lbf, s_cs_purity)])

cat("\n── FINEMAP (k=5) summary ──\n")
print(summary_dt[, .(locus, f_n_cred, f_max_pip, f_top_snp, f_exp_k, f_pr_k5, f_best_prob)])

cat("\n── Top-SNP agreement ──\n")
both <- summary_dt[!is.na(s_max_pip) & !is.na(f_max_pip)]
if (nrow(both) > 0) {
  both[, agree := s_top_snp == f_top_snp]
  cat(sprintf("Agreement: %d / %d\n", sum(both$agree, na.rm = TRUE), nrow(both)))
  if (any(!both$agree, na.rm = TRUE))
    print(both[agree == FALSE, .(locus, s_top_snp, s_max_pip, f_top_snp, f_max_pip)])
}

out_file <- file.path(DATADIR, "finemapping_bd1_summary.tsv")
fwrite(summary_dt, out_file, sep = "\t")
cat("\nSummary written to:", out_file, "\nDone.\n")
