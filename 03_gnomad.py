"""Query gnomAD v4 (GRCh38) allele counts for the panel candidates."""
import time, requests, pandas as pd

Q = """query($v:String!){variant(variantId:$v,dataset:gnomad_r4){
        exome{ac an homozygote_count} genome{ac an homozygote_count}}}"""

def gnomad(vid, tries=5):
    for i in range(tries):
        try:
            r = requests.post("https://gnomad.broadinstitute.org/api",
                              json={"query": Q, "variables": {"v": vid}},
                              headers={"User-Agent": "Mozilla/5.0"}, timeout=60)
            if r.status_code == 200 and r.text.strip().startswith("{"):
                j = r.json()
                return None if "errors" in j and not j.get("data") else j["data"]["variant"]
            time.sleep(15 * (i + 1))
        except Exception:
            time.sleep(10)
    return None

df = pd.read_csv("candidates_panel.csv")
res = []
for _, x in df.iterrows():
    v = gnomad(f"{x.chrom}-{x.pos}-{x.ref}-{x.alt}")
    ex, ge = (v or {}).get("exome") or {}, (v or {}).get("genome") or {}
    ac = (ex.get("ac") or 0) + (ge.get("ac") or 0)
    an = (ex.get("an") or 0) + (ge.get("an") or 0)
    res.append({"gen": x.gen, "cambio": x.cambio_proteina, "en_gnomAD": v is not None,
                "AC": ac, "AF": ac / an if an else None})
    time.sleep(3)
pd.DataFrame(res).sort_values("AF", na_position="first").to_csv("gnomad_panel.csv", index=False)
print(pd.DataFrame(res).sort_values("AF", na_position="first").to_string())
