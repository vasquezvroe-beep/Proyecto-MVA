"""Extract HIGH/MODERATE variants and ClinVar pathogenic hits into local CSVs.
Output files contain variant-level data: keep them local (they are gitignored)."""
import re, subprocess
import pandas as pd

VCF = "qc_ann_clinvar.vcf.gz"
FMT = "%CHROM\t%POS\t%REF\t%ALT\t%ID\t[%GT\t%DP\t%GQ]\t%INFO/CLNSIG\t%INFO/ANN\n"
out = subprocess.run(["bcftools", "query", "-f", FMT, VCF], capture_output=True, text=True, check=True).stdout

rows = []
for line in out.splitlines():
    c = line.split("\t")
    if len(c) < 10 or not re.search(r"\|(HIGH|MODERATE)\|", c[9]):
        continue
    for a in c[9].split(","):
        x = a.split("|")
        if len(x) > 10 and x[2] in ("HIGH", "MODERATE"):
            rows.append(dict(chrom=c[0], pos=int(c[1]), ref=c[2], alt=c[3], rsID=c[4],
                             gen=x[3], efecto=x[1], impacto=x[2], cambio_proteina=x[10],
                             genotipo=c[5], DP=c[6], GQ=c[7], CLNSIG=c[8]))
            break

df = pd.DataFrame(rows)
df["zigosidad"] = df.genotipo.map(lambda g: "homocigoto" if g in ("1/1", "2/2") else "heterocigoto")
df.to_csv("candidates_all.csv", index=False)
print("HIGH/MODERATE variants:", len(df))

panel = [g.strip() for g in open("mva_panel.txt") if g.strip()]
df[df.gen.isin(panel)].to_csv("candidates_panel.csv", index=False)
print("Variants in panel genes:", int(df.gen.isin(panel).sum()))
print("ClinVar pathogenic/likely pathogenic:",
      int(df.CLNSIG.fillna("").str.contains("athogenic") .sum()))
