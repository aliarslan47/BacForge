#!/usr/bin/env python3
"""GBK'dan 2 sutunlu CDS->isim anotasyon tablosu uretir (VFDB onceligi ile)."""
import sys
from Bio import SeqIO

gbk, out = sys.argv[1], sys.argv[2]
rows = []
for rec in SeqIO.parse(gbk, "genbank"):
    for f in rec.features:
        if f.type != "CDS":
            continue
        lt = f.qualifiers.get("locus_tag", [""])[0]
        vfdb = f.qualifiers.get("vfdb_description", [""])[0].strip()
        prod = f.qualifiers.get("product", [""])[0].strip()
        rows.append((lt, vfdb if vfdb else prod))

with open(out, "w") as fh:
    fh.write("CDS\tisim\n")
    for lt, name in rows:
        fh.write(f"{lt}\t{name}\n")

print(f"{'CDS':<20} isim")
print("-" * 60)
for lt, name in rows:
    print(f"{lt:<20} {name}")
print("-" * 60)
print(f"Toplam {len(rows)} CDS -> {out}")
