#!/usr/bin/env bash
set -euo pipefail

# =========================
# EDIT THESE PATHS / OPTIONS
# =========================
VCF=/n/home08/dwhite/scripts/psmc/05_variant_calling/vcf_ilmn/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.ilmn.vcf.gz   
             # single- or multi-sample VCF (indexed)
SAMPLE=Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton              # sample ID to evaluate (for multi-sample VCF)
CALLABLE_BED=/n/home08/dwhite/scripts/genome-wide-Ho/callable.bed     # from 02_build_callable_bed.sh
OUT_TSV=heterozygosity_summary.tsv     # output summary row appended/created
# genotype quality / depth thresholds (optional; set to 0 to disable)
MIN_GQ=0
MIN_DP=0
# =========================

# --- sanity checks ---
for f in "$VCF" "$CALLABLE_BED"; do
  [[ -f "$f" ]] || { echo "[ERR] Missing file: $f" >&2; exit 2; }
done
[[ -f "${VCF}.tbi" || -f "${VCF}.csi" ]] || { echo "[ERR] VCF index (.tbi/.csi) missing for: $VCF" >&2; exit 3; }

# 1) helper: callable length
CALLABLE=$(awk '{s+=($3-$2)} END{print s}' "${CALLABLE_BED}")
[[ "$CALLABLE" -gt 0 ]] || { echo "[ERR] Callable length is zero"; exit 4; }

echo "[DBG] Sample names in VCF:"
bcftools query -l "$VCF" | tr '\n' ' ' && echo
if ! bcftools query -l "$VCF" | grep -qx "$SAMPLE"; then
  echo "[ERR] SAMPLE not found in VCF: $SAMPLE" >&2; exit 5
fi

echo "[DBG] Total PASS or unfiltered (.) biallelic SNPs in whole VCF:"
bcftools view -m2 -M2 -v snps -f .,PASS -H "$VCF" | wc -l

echo "[DBG] Total PASS or . biallelic SNPs within callable (using -T BED):"
bcftools view -T "$CALLABLE_BED" -m2 -M2 -v snps -f .,PASS -H "$VCF" | wc -l | sed 's/^/  /'


# 2) subset to sample, callable regions, PASS, biallelic SNPs
#    apply optional per-genotype filters via bcftools query -i expression
bcftools view -s "${SAMPLE}" -T "${CALLABLE_BED}" -m2 -M2 -v snps -f .,PASS -Ou "${VCF}" > _tmp.snps.bcf

expr_base='GT="het"'
expr_hom='GT=="1/1" || GT=="1|1"'
# Add GQ/DP thresholds if requested
if [[ "$MIN_GQ" -gt 0 ]]; then
  expr_base="${expr_base} && FMT/GQ>=${MIN_GQ}"
  expr_hom="${expr_hom} && FMT/GQ>=${MIN_GQ}"
fi
if [[ "$MIN_DP" -gt 0 ]]; then
  expr_base="${expr_base} && FMT/DP>=${MIN_DP}"
  expr_hom="${expr_hom} && FMT/DP>=${MIN_DP}"
fi

# 3) counts
HETS=$(bcftools query -i "${expr_base}" -f '%CHROM\t%POS\n' _tmp.snps.bcf | wc -l | awk '{print $1}')
HOMS=$(bcftools query -i "${expr_hom}"  -f '%CHROM\t%POS\n' _tmp.snps.bcf | wc -l | awk '{print $1}')
TOTAL_VAR=$((HETS + HOMS))

# 4) heterozygosity metrics
python3 - "$CALLABLE" "$HETS" "$HOMS" <<'PY'
import sys
callable=int(sys.argv[1]); hets=int(sys.argv[2]); homs=int(sys.argv[3])
h = hets/callable if callable else float('nan')
print(f"Callable bases           : {callable:,}")
print(f"Heterozygous SNPs        : {hets:,}")
print(f"Homozygous-ALT SNPs      : {homs:,}")
print(f"Heterozygosity per site  : {h:.6e}")
print(f"Heterozygosity per kb    : {h*1e3:.3f}")
print(f"Heterozygosity per Mb    : {h*1e6:.1f}")
PY

# 5) write/append a summary row
if [[ ! -f "${OUT_TSV}" ]]; then
  echo -e "sample\tvcf\tcallable_bases\thets\thoms\ttotal_var\thet_per_site" > "${OUT_TSV}"
fi
HET_PER_SITE=$(python3 - <<PY
callable=${CALLABLE}; hets=${HETS}
print(f"{(hets/callable):.6e}" if callable else "nan")
PY
)
echo -e "${SAMPLE}\t${VCF}\t${CALLABLE}\t${HETS}\t${HOMS}\t${TOTAL_VAR}\t${HET_PER_SITE}" >> "${OUT_TSV}"

echo "[OK] Appended summary to ${OUT_TSV}"
rm -f _tmp.snps.bcf

