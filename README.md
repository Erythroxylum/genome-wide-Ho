# genome-wide-Ho
This repository contains shell scripts to calculate genome-wide heterozygosity (Ho) from whole-genome sequencing data. The pipeline extracts regions of Ns from a consensus fasta sequence and combines those regions with other coverage masks (such as those generated in the github.com/Erythroxylum/psmc-pipeline) and then uses BEDTools and BCFtools to mask non-callable regions (N's) from a fasta, extract SNPs, and compute heterozygosity metrics.

## Requirements

-BCFtools ≥1.20

-BEDTools ≥2.31.1

-SAMtools

-Python ≥3.7

-Bash shell (tested on Linux HPC environments)

Make sure these tools are available in your PATH.

## Scripts
-01_add-N-to-mask-bed.sh – extracts regions of Ns from the consensus assembly and combines with existing coverage masks.

-02_build_callable_bed.sh – creates a BED file of callable sites (genome minus masked regions).

-03_count-variants-Ho.sh – counts heterozygous and homozygous-ALT SNPs and computes heterozygosity statistics.

## How to Run
Each script defines input/output paths at the top. Edit these variables before running:

### Example from 03_count-variants-Ho.sh
````
VCF=/path/to/your/sample.vcf.gz         # input VCF (bgzipped + indexed)
SAMPLE=SampleName                       # sample name in VCF
CALLABLE_BED=/path/to/callable.bed      # output from step 02
OUT_TSV=heterozygosity_summary.tsv      # results file
````

Once paths are set, run the script directly:
````
bash 01_add-N-to-mask-bed.sh
bash 02_build_callable_bed.sh
bash 03_count-variants-Ho.sh
````

The final output will be a tab-delimited file (heterozygosity_summary.tsv) with summary statistics per sample.

## Notes
-Always check that your VCF is indexed (.tbi or .csi).
-The pipeline expects bgzipped VCFs. Use bcftools view -Oz -o out.vcf.gz in.vcf to compress.

You can adjust genotype quality (GQ) and depth (DP) thresholds in 03_count-variants-Ho.sh (default = no filter).

For multiple samples, run step 03 separately for each sample name.
