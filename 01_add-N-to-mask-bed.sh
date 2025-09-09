# ---- inputs (edit) ----
REF=/n/netscratch/davis_lab/Everyone/dwhite/psmc/06_consensus/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.iupac.fa
CONSFA=/n/netscratch/davis_lab/Everyone/dwhite/psmc/06_consensus/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.iupac.fa
MASK_STEP6=/n/home08/dwhite/scripts/psmc/06_consensus/coverage_masks_Pfurb2/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.mask.bed   # from home, step 6
KEEP_CHROMS=                 # optional; one contig per line
SAMPLE=Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton
VCF=/n/home08/dwhite/scripts/psmc/05_variant_calling/vcf_ilmn/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.ilmn.vcf.gz          # or vcf_hifi/...

# 0) index the reference if needed
[ -f "${REF}.fai" ] || samtools faidx "${REF}"

# 1) extract contiguous N-runs from the consensus as BED (0-based, half-open)
awk 'BEGIN{OFS="\t"}
  /^>/{
    if (seq!="") { # flush trailing N-run
      if (inN) print name, start, pos
    }
    name=substr($0,2); seq=""; pos=0; inN=0
    next
  }
  {
    line=toupper($0)
    for (i=1; i<=length(line); i++) {
      c=substr(line,i,1)
      if (c=="N") {
        if (!inN) { start=pos; inN=1 }
      } else {
        if (inN) { print name, start, pos; inN=0 }
      }
      pos++
    }
  }
  END{
    if (inN) print name, start, pos
  }' "${CONSFA}" > n_runs_from_consensus.bed

# 2) (optional) restrict to chosen contigs up front
if [[ -s "$KEEP_CHROMS" ]]; then
  awk 'NR==FNR{ok[$1]=1;next} ok[$1]' "$KEEP_CHROMS" n_runs_from_consensus.bed > n_runs.keep.bed
  mv n_runs.keep.bed n_runs_from_consensus.bed
fi

# 3) combine the N-runs with your existing mask, then sort+merge
cat "${MASK_STEP6}" n_runs_from_consensus.bed \
  | bedtools sort -i - \
  | bedtools merge -i - \
  > mask.combined.bed

