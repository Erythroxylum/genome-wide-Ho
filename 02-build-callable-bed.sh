#!/usr/bin/env bash
set -euo pipefail

# =========================
# EDIT THESE PATHS
# =========================
REF=/n/netscratch/davis_lab/Everyone/dwhite/psmc/06_consensus/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.iupac.fa                       # reference used for mapping/calling
COMBINED_MASK=/n/home08/dwhite/scripts/genome-wide-Ho/mask.combined.bed            # output from 01_add-N-to-mask-bed.sh
KEEP_CHROMS=                                        # optional: file with contigs to KEEP (one per line) or leave empty
CALLABLE_BED=callable.bed                           # output BED
# =========================

# --- sanity checks ---
for f in "$REF" "$COMBINED_MASK"; do
  [[ -f "$f" ]] || { echo "[ERR] Missing file: $f" >&2; exit 2; }
done
[[ -z "${KEEP_CHROMS}" || -s "${KEEP_CHROMS}" ]] || { echo "[ERR] KEEP_CHROMS set but empty: ${KEEP_CHROMS}" >&2; exit 3; }

# 1) reference genome lengths and order
[[ -f "${REF}.fai" ]] || samtools faidx "${REF}"
cut -f1,2 "${REF}.fai" > ref.genome

# 2) sort/merge the mask to match ref.genome order
bedtools sort  -g ref.genome -i "${COMBINED_MASK}" > _mask.sorted.bed
bedtools merge -i _mask.sorted.bed                > _mask.sorted.merged.bed

# 3) build "keep" intervals across full contigs (either selected contigs or all)
if [[ -n "${KEEP_CHROMS}" ]]; then
  awk 'NR==FNR{ok[$1]=1; next} ok[$1]{print $1,0,$2}' OFS='\t' "${KEEP_CHROMS}" ref.genome > keep.bed
else
  awk '{print $1,0,$2}' OFS='\t' ref.genome > keep.bed
fi

# 4) callable = kept contigs minus mask, then limit back to kept contigs
bedtools complement -i _mask.sorted.merged.bed -g ref.genome \
  | bedtools intersect -a - -b keep.bed -wa \
  > "${CALLABLE_BED}"

CALLABLE=$(awk '{s+=($3-$2)} END{print s}' "${CALLABLE_BED}")
echo "[OK] Wrote ${CALLABLE_BED}"
echo "[OK] Callable bases: ${CALLABLE}"

# optional diagnostics: contigs in mask but not in reference
comm -23 <(cut -f1 _mask.sorted.bed | sort -u) <(cut -f1 ref.genome | sort -u) | head -n 5 | sed 's/^/[WARN] Mask contig not in ref: /' || true

rm -f _mask.sorted.bed _mask.sorted.merged.bed
