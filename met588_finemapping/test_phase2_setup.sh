#!/bin/bash
# Phase 2 setup: reproduce LD for norcloz_chr2 using the pipeline approach
# (PLINK2 extract + orient alleles, PLINK1.9 pairwise --r, then LDmerge)
# Runs inside WSL; called by test_phase2.R via system()

set -e

SCZVSBP="/home/c1977426/sczvsbp"
DATADIR="$SCZVSBP/met588_finemapping/data"
RESULTDIR="$SCZVSBP/met588_finemapping/results"
WORKDIR="$SCZVSBP/met588_finemapping/test_phase2_work"
PLINK1="$SCZVSBP/met588_finemapping/software/plink"
PLINK2="$SCZVSBP/met588_finemapping/software/plink2"

chmod +x "$PLINK1" "$PLINK2"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

LOCUS="norcloz_chr2_test"
COHORT="CLOZUK2"
BFILE="$DATADIR/CLOZUK2_chr2.qc_clozlevels"

echo "=== Phase 2: LD pipeline test for norcloz_chr2 ==="
echo "Working dir: $WORKDIR"
echo "PLINK1: $PLINK1"
echo "PLINK2: $PLINK2"

# ── Step 1: Build .snplist, .ref, and .a1 from .z file ──────────────────────
# .snplist = one rsid per line (for --extract)
# .ref     = "rsid allele2" per line (for PLINK2 --ref-allele: sets REF = non-effect allele)
# .a1      = "rsid allele1" per line (for PLINK1.9 --a1-allele: ensures correct LD sign)
echo ""
echo "--- Step 1: Build SNP list, ref file, and a1 file ---"
awk 'NR>1 {print $1}'      "$RESULTDIR/norcloz_chr2.z" > "${LOCUS}.snplist"
awk 'NR>1 {print $1, $5}'  "$RESULTDIR/norcloz_chr2.z" > "${LOCUS}.ref"
awk 'NR>1 {print $1, $4}'  "$RESULTDIR/norcloz_chr2.z" > "${LOCUS}.a1"
echo "  SNPs in list: $(wc -l < ${LOCUS}.snplist)"

# ── Step 2: PLINK2 — extract region, orient alleles, make BED ────────────────
echo ""
echo "--- Step 2: PLINK2 extract + allele orient ---"
$PLINK2 \
    --bfile "$BFILE" \
    --extract "${LOCUS}.snplist" \
    --ref-allele "${LOCUS}.ref" \
    --rm-dup force-first \
    --memory 4000 --threads 2 \
    --make-bed \
    --out "${COHORT}_${LOCUS}"

nvars=$(wc -l < "${COHORT}_${LOCUS}.bim")
echo "  Variants after filter: $nvars"

# ── Step 3: PLINK1.9 — pairwise signed LD (R) ────────────────────────────────
echo ""
echo "--- Step 3: PLINK1.9 pairwise --r ---"
$PLINK1 \
    --bfile "${COHORT}_${LOCUS}" \
    --a1-allele "${LOCUS}.a1" 2 1 \
    --r \
    --ld-window 999999 \
    --ld-window-kb 99999 \
    --ld-window-r2 0 \
    --memory 4000 --threads 2 \
    --out "${COHORT}_${LOCUS}"

npairs=$(wc -l < "${COHORT}_${LOCUS}.ld")
echo "  LD pairs: $npairs"

# ── Step 4: Compress and keep fam ────────────────────────────────────────────
echo ""
echo "--- Step 4: Compress LD, copy fam ---"
gzip -9 -f "${COHORT}_${LOCUS}.ld"
cp "$BFILE.fam" "${COHORT}_${LOCUS}.fam"
rm -f "${COHORT}_${LOCUS}".{bed,bim,log,nosex}

echo ""
echo "Done. Files in $WORKDIR:"
ls -lh "$WORKDIR"
