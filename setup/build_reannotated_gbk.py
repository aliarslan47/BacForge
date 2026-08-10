#!/usr/bin/env python3
"""
build_reannotated_gbk.py — Pharokka GBK'yi reanotasyon tablosuyla zenginlestirip
Proksee/CGView icin ETIKET-DOSTU bir GBK uretir.

NE YAPAR:
  1) reannotation_tablo.tsv'deki (locus_tag -> fonksiyon) atamalarini GBK'ya isler
     (ilgili CDS'lerin /product alanini gunceller).
  2) Fonksiyonu BILINEN her CDS'e kisa, okunur bir /gene etiketi ekler.
     -> Proksee gene'i locus_tag'e tercih ettiginden haritada FONKSIYON ADI gorunur.
  3) Hipotetik CDS'lere DOKUNMAZ (/gene yok) -> haritada sade kalir.

Kullanim:
  conda run -n ali-pharokka python3 setup/build_reannotated_gbk.py \
      --gbk  <contig.gbk> \
      --reann <reannotation_tablo.tsv>   (opsiyonel) \
      --out  <contig_REANNOTATED.gbk>
"""
import argparse, csv, re, sys
from pathlib import Path
from Bio import SeqIO

HYPO = "hypothetical protein"

# Uzun urun adlarini harita etiketi icin kisaltma kurallari (sira onemli).
SHORTEN = [
    (r"crossover junction endodeoxyribonuclease", "resolvase"),
    (r"\bfamily protein\b", ""),
    (r"\bfamily\b", ""),
    (r"\bputative\b", ""),
    (r"single strand DNA binding protein", "SSB"),
    (r"\bhomologue\b", "homolog"),
    (r"DnaD-like helicase loader", "DnaD loader"),
    (r"metal-dependent hydrolase", "metallohydrolase"),
    (r"Panton-Valentine leukoci\w+", "PVL"),
]

def short_label(product: str) -> str:
    s = product.strip()
    for pat, rep in SHORTEN:
        s = re.sub(pat, rep, s, flags=re.IGNORECASE)
    s = re.sub(r"\s{2,}", " ", s).strip(" -")
    # cok uzun kalirsa kelime bazinda kirp
    if len(s) > 28:
        s = " ".join(s.split()[:3])
    return s

def load_reann(path):
    m = {}
    if not path:
        return m
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            lt = (row.get("protein_id") or "").strip()
            fn = (row.get("atanan_fonksiyon") or "").strip()
            if lt and fn:
                m[lt] = fn
    return m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gbk", required=True)
    ap.add_argument("--reann", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    reann = load_reann(a.reann)
    recs = list(SeqIO.parse(a.gbk, "genbank"))
    n_reann = n_gene = n_hypo = n_cds = n_vfdb = 0

    for rec in recs:
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            n_cds += 1
            lt = (feat.qualifiers.get("locus_tag", [""])[0]).strip()
            # 1) reanotasyon atamasi
            if lt in reann:
                feat.qualifiers["product"] = [reann[lt]]
                n_reann += 1
            product = (feat.qualifiers.get("product", [""])[0]).strip()
            # 2a) VFDB hit varsa etiket VIRULANS faktoru olsun (Pharokka'nin
            #     genel adina ONCELIKLI) -> ornek: "staphylokinase (sak)"
            sn = (feat.qualifiers.get("vfdb_short_name", [""])[0]).strip()
            desc = (feat.qualifiers.get("vfdb_description", [""])[0]).strip()
            if sn:
                base = short_label(desc) if desc else sn
                feat.qualifiers["gene"] = [f"{base} ({sn})"]
                n_vfdb += 1
                n_gene += 1
                continue
            # 2b) fonksiyonel ise /gene; hipotetik ise dokunma
            if product and product.lower() != HYPO:
                feat.qualifiers["gene"] = [short_label(product)]
                n_gene += 1
            else:
                feat.qualifiers.pop("gene", None)
                n_hypo += 1

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(recs, a.out, "genbank")
    print(f"CDS toplam      : {n_cds}")
    print(f"reanote edilen  : {n_reann}")
    print(f"VFDB virulans   : {n_vfdb}")
    print(f"etiketli (gene) : {n_gene}")
    print(f"hipotetik (sade): {n_hypo}")
    print(f"-> {a.out}")

if __name__ == "__main__":
    main()
