#!/usr/bin/env bash
# Usage: ./01_qc_annotate.sh input.vcf.gz
# Requires: bcftools, tabix, Java 21+, SnpEff (https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip)
# Run in a temporary environment and delete all data afterwards.
set -euo pipefail
VCF="$1"

tabix -f -p vcf "$VCF"

# 1) Quality filter
bcftools view -f PASS -i 'FMT/DP>=10 && FMT/GQ>=20 && GT="alt"' "$VCF" -Oz -o qc.vcf.gz
tabix -f -p vcf qc.vcf.gz
bcftools stats qc.vcf.gz | grep '^SN'

# 2) Annotation (strip Java warning lines that start with '[')
java -Xmx8g -jar snpEff/snpEff.jar download GRCh38.mane.1.5.ensembl
java -Xmx8g -jar snpEff/snpEff.jar eff -noStats GRCh38.mane.1.5.ensembl qc.vcf.gz 2>/dev/null \
  | grep -v '^\[' | bgzip > qc_ann.vcf.gz
tabix -f -p vcf qc_ann.vcf.gz

# 3) ClinVar (GRCh38)
wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz -O clinvar.vcf.gz
wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi -O clinvar.vcf.gz.tbi
bcftools annotate -a clinvar.vcf.gz -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT \
  qc_ann.vcf.gz -Oz -o qc_ann_clinvar.vcf.gz
tabix -f -p vcf qc_ann_clinvar.vcf.gz
