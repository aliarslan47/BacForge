#!/usr/bin/env python3
"""GENEL FİLOGENİ (genom-geneli) — mash k-mer mesafesi NJ ağacı.
Marker-bağımsız: 33 fajimiz + 19 referans, TAMAMI. TerL ML ağacını tamamlar.
(Referans-çıpalı; asla sadece birbiriyle değil.)"""
import subprocess, itertools
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INP = ROOT/"runs"/"_viptree_input"   # 33 QUERY + 19 REF .fna (önceden hazır)
SBU = Path("/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj")
OUT = SBU/"YAYIN_PAKETI_TUM"/"09_Genel_Filogeni_MASH"; OUT.mkdir(parents=True, exist_ok=True)
ENV = "ali-pharokka"  # mash + biopython burada

def cr(cmd, t=900):
    return subprocess.run(["conda","run","-n",ENV]+cmd, capture_output=True, text=True, timeout=t)

files = sorted(INP.glob("*.fna"))
print(f"genom: {len(files)} (QUERY+REF)")
# 1) mash sketch (hepsi tek msh)
msh = OUT/"all.msh"
cr(["mash","sketch","-s","20000","-k","17","-o",str(OUT/"all")] + [str(f) for f in files])
# 2) pairwise dist
d = cr(["mash","dist",str(msh),str(msh)])
(OUT/"mash_dist.tsv").write_text(d.stdout)
# 3) mesafe matrisi + NJ ağaç (biopython)
labels=[f.stem for f in files]
idx={f.stem:i for i,f in enumerate(files)}
import re
n=len(labels); M=[[0.0]*n for _ in range(n)]
for ln in d.stdout.splitlines():
    c=ln.split("\t")
    if len(c)>=3:
        a=Path(c[0]).stem; b=Path(c[1]).stem
        if a in idx and b in idx: M[idx[a]][idx[b]]=float(c[2])
# alt-üçgen
lower=[[M[i][j] for j in range(i+1)] for i in range(n)]
NJ = '''
import sys, json
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
data=json.load(open(sys.argv[1]))
dm=DistanceMatrix(names=data["labels"], matrix=data["lower"])
tree=DistanceTreeConstructor().nj(dm)
tree.ladderize()
Phylo.write(tree, sys.argv[2], "newick")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig=plt.figure(figsize=(13,max(4,len(data["labels"])*0.4))); ax=fig.add_subplot(111)
def lab(c):
    if not c.is_terminal(): return ""
    nm=(c.name or "")
    return ("▶ "+nm.replace("QUERY_","")) if nm.startswith("QUERY_") else nm.replace("REF_","")
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: lab(c).replace("_"," "))
for t in ax.texts:
    if t.get_text().startswith("▶"): t.set_color("crimson"); t.set_fontweight("bold")
ax.set_title("GENEL FILOGENI (genom-geneli, mash) | ▶=bizim fajlar, digerleri REFERANS", fontsize=12, fontweight="bold")
ax.set_xlabel("mash genom mesafesi (k-mer) | Neighbor-Joining")
for s in ["top","right","left"]: ax.spines[s].set_visible(False)
ax.get_yaxis().set_visible(False)
plt.tight_layout(); plt.savefig(sys.argv[3], dpi=200, bbox_inches="tight"); print("PNG OK")
'''
import json
(OUT/"_njdata.json").write_text(json.dumps({"labels":labels,"lower":lower}))
Path("/tmp/nj_tree.py").write_text(NJ)
r=cr(["python","/tmp/nj_tree.py",str(OUT/"_njdata.json"),str(OUT/"mash_nj.treefile"),str(OUT/"genel_filogeni_mash.png")],300)
print("NJ+render:", r.stdout.strip() or r.stderr[-300:])
# ascii
if (OUT/"mash_nj.treefile").exists():
    a=cr(["python","-c",f"from Bio import Phylo;t=Phylo.read('{OUT/'mash_nj.treefile'}','newick');Phylo.draw_ascii(t)"],120)
    (OUT/"genel_filogeni_mash_ascii.txt").write_text(a.stdout or a.stderr)
print("MASH GENEL FİLOGENİ BİTTİ:", OUT)
