#!/bin/bash
# PRE-FIX version of test_phase2_setup.sh
# Replicates the OLD pipeline behaviour: no --ref-allele in PLINK2,
# no --a1-allele in PLINK1.9.  LD signs are therefore not oriented
# relative to the GWAS effect allele — reproducing the original bug.

set -e

SCZVSBP="/home/c1977426/sczvsbp"
DATADIR="$SCZVSBP/met588_finemapping/data"
RESULTDIR="$SCZVSBP/met588_finemapping/results"
WORKDIR="$SCZVSBP/met588_finemapping/test_phase2_work_pre_fix"
PLINK1="$SCZVSBP/met588_finemapping/software/plink"
PLINK2="$SCZVSBP/met588_finemapping/software/plink2"

chmod +x "$PLINK1" "$PLINK2"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

LOCUS="norcloz_chr2_test"
COHORT="CLOZUK2"
BFILE="$DATADIR/CLOZUK2_chr2.qc_clozlevels"

echo "=== PRE-FIX: LD pipeline (no allele orientation flags) ==="
echo "Working dir: $WORKDIR"

# ── Step 1: Build .snplist only (no .ref or .a1 used in PLINK calls) ─────────
echo ""
echo "--- Step 1: Build SNP list ---"
awk 'NR>1 {print $1}'      "$RESULTDIR/norcloz_chr2.z" > "${LOCUS}.snplist"
# Still build .ref and .a1 for LDmerge metadata but do NOT pass to PLINK
awk 'NR>1 {print $1, $5}'  "$RESULTDIR/norcloz_chr2.z" > "${LOCUS}.ref"
awk 'NR>1 {print $1, $4}'  "$RESULTDIR/norcloz_chr2.z" > "${LOCUS}.a1"
echo "  SNPs in list: $(wc -l < ${LOCUS}.snplist)"

# ── Step 2: PLINK2 — extract region, NO --ref-allele ─────────────────────────
echo ""
echo "--- Step 2: PLINK2 extract (NO --ref-allele) ---"
$PLINK2 \
    --bfile "$BFILE" \
    --extract "${LOCUS}.snplist" \
    --rm-dup force-first \
    --memory 4000 --threads 2 \
    --make-bed \
    --out "${COHORT}_${LOCUS}"

nvars=$(wc -l < "${COHORT}_${LOCUS}.bim")
echo "  Variants after filter: $nvars"

# ── Step 3: PLINK1.9 — pairwise --r, NO --a1-allele ─────────────────────────
echo ""
echo "--- Step 3: PLINK1.9 pairwise --r (NO --a1-allele) ---"
$PLINK1 \
    --bfile "${COHORT}_${LOCUS}" \
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
