#!/usr/bin/env python3
"""Cinssiz/novel query'ler için KONAK-bağlam ICTV ağacı (host fajları [Title] + en yakın hit + kardeşler)."""
import sys, subprocess, json, shutil
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_per_query_trees import genus_members, fetch_fasta, latest, clean, ENVP

runs = Path(__file__).resolve().parents[1]/"runs"
Path("/tmp/nj_per_query.py").write_text((Path(__file__).resolve().parents[1]/"setup"/"_nj_template.py").read_text()
    if (Path(__file__).resolve().parents[1]/"setup"/"_nj_template.py").exists() else "")

NJ='''
import sys, json
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
d=json.load(open(sys.argv[1])); focal=d["focal"]
dm=DistanceMatrix(names=d["labels"], matrix=d["lower"])
tree=DistanceTreeConstructor().nj(dm); tree.ladderize()
fig=plt.figure(figsize=(12,max(3,len(d["labels"])*0.45))); ax=fig.add_subplot(111)
def lab(c):
    if not c.is_terminal(): return ""
    nm=c.name or ""
    if nm==focal: return "★ "+nm.replace("QUERY_","")
    if nm.startswith("QUERY_"): return "▶ "+nm.replace("QUERY_","")
    return nm.replace("REF_","")
Phylo.draw(tree,axes=ax,do_show=False,label_func=lambda c: lab(c).replace("_"," "))
for t in ax.texts:
    s=t.get_text()
    if s.startswith("★"): t.set_color("crimson"); t.set_fontweight("bold"); t.set_fontsize(11)
    elif s.startswith("▶"): t.set_color("darkorange"); t.set_fontweight("bold")
ax.set_title(d["title"], fontsize=11, fontweight="bold")
ax.set_xlabel("mash genom mesafesi (k-mer) NJ | ★=bu query, ▶=ayni grup, digerleri host referans")
for sp in ["top","right","left"]: ax.spines[sp].set_visible(False)
ax.get_yaxis().set_visible(False)
plt.tight_layout(); plt.savefig(sys.argv[2],dpi=190,bbox_inches="tight"); print("OK")
'''
Path("/tmp/nj_novel.py").write_text(NJ)

def build_tree(genomes, focal, title, out_png, workdir):
    workdir.mkdir(parents=True, exist_ok=True); fnas=[]
    for lab,fa in genomes:
        p=workdir/f"{lab}.fna"
        seq="".join(l for l in Path(fa).read_text().splitlines() if not l.startswith(">"))
        p.write_text(f">{lab}\n{seq}\n"); fnas.append(p)
    subprocess.run(["conda","run","-n",ENVP,"mash","sketch","-s","20000","-k","17","-o",str(workdir/"all")]+[str(f) for f in fnas],capture_output=True,timeout=600)
    d=subprocess.run(["conda","run","-n",ENVP,"mash","dist",str(workdir/"all.msh"),str(workdir/"all.msh")],capture_output=True,text=True,timeout=600)
    labels=[f.stem for f in fnas]; idx={l:i for i,l in enumerate(labels)}; n=len(labels)
    M=[[0.0]*n for _ in range(n)]
    for ln in d.stdout.splitlines():
        c=ln.split("\t")
        if len(c)>=3:
            a=Path(c[0]).stem; b=Path(c[1]).stem
            if a in idx and b in idx: M[idx[a]][idx[b]]=float(c[2])
    lower=[[M[i][j] for j in range(i+1)] for i in range(n)]
    jd=workdir/"_d.json"; jd.write_text(json.dumps({"labels":labels,"lower":lower,"focal":focal,"title":title}))
    subprocess.run(["conda","run","-n",ENVP,"python","/tmp/nj_novel.py",str(jd),str(out_png)],capture_output=True,timeout=300)
    shutil.rmtree(workdir, ignore_errors=True)

JOBS=[("2858","contig_1","Enterococcus","OP072375.1"),
      ("21857478","contig_228","Pseudomonas","PZ179922.1"),
      ("21857478","contig_273","Enterococcus","NC_112203.1"),
      ("21857478","contig_274","Enterococcus","OP072375.1"),
      ("21663260","contig_186","Enterococcus","OP072375.1"),
      ("4188mrsa","contig_88","Enterococcus","OP172804.1")]
mag=[("2858","contig_1"),("21857478","contig_274"),("21663260","contig_186")]
for s,c,host,acc in JOBS:
    d=latest(s)/"18_Final_Report"/"yayina_uygun"/c
    cfa=d/"01_Genome"/f"{c}.fasta"
    if not cfa.exists(): cfa=d/f"{c}.fasta"
    focal=f"QUERY_{s}_{c}_novel"; genomes=[(focal,cfa)]
    if (s,c) in mag:
        for s2,c2 in mag:
            if (s2,c2)!=(s,c):
                f2=latest(s2)/"18_Final_Report"/"yayina_uygun"/c2/"01_Genome"/f"{c2}.fasta"
                if f2.exists(): genomes.append((f"QUERY_{s2}_{c2}_novel",f2))
    rf=fetch_fasta(acc)
    if rf:
        t=open(rf).readline()[1:].strip(); genomes.append((f"REF_{clean(t)}",rf))
    mem=genus_members(f'"{host} phage"[Title] AND complete genome AND viruses[filter]', 12)
    seen={acc}
    for a,t,f in mem:
        if a not in seen: genomes.append((f"REF_{clean(t)}",f)); seen.add(a)
    if len(genomes)<3:
        print(f"  {s}/{c}: yetersiz ({len(genomes)})"); continue
    o=d/"09_Phylogeny"; o.mkdir(parents=True,exist_ok=True)
    build_tree(genomes, focal, f"{s}/{c} ICTV yerlesim (cinssiz/novel) — host {host} baglami",
               o/f"{c}_ICTV_agac.png", runs/"_pq_novel"/f"{s}_{c}")
    print(f"  {s}/{c}: {len(genomes)} genom (host {host})", flush=True)
print("NOVEL host-baglam agaclari BITTI")
