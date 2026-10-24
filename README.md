# MVA Hackathon 2026: Track 1 pipeline

Team: Rocio Estefani Vasquez Villalba, Luana Beltran, Cesar Perez

Phenotype-driven prioritisation of variants from a single-sample WGS VCF.

**This repository contains no patient data**: no VCF/BAM/FASTQ, no clinical document and no
variant-level tables.

## Steps
1. `01_qc_annotate.sh input.vcf.gz`: quality filter, SnpEff annotation
   (GRCh38.mane.1.5.ensembl), ClinVar annotation.
2. `02_candidates.py`: HIGH/MODERATE variants, panel genes (`mva_panel.txt`),
   ClinVar pathogenic hits.
3. `03_gnomad.py`: gnomAD v4 allele counts for panel candidates.
4. Manual step: missense effect via Ensembl VEP REST (AlphaMissense, CADD); ranking of
   candidates by population rarity, ClinVar status and phenotype fit.

Requirements: bcftools, tabix, Java 21+, SnpEff, Python 3 with pandas and requests.

## Report
See the Track 1 report submitted with the predictions. AI assistance is disclosed there.

## Data handling
All data were processed in a temporary environment and will be deleted after the
hackathon, as required by the data use terms.
