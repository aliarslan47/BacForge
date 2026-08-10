#!/usr/bin/env python3
"""GENEL FİLOGENİ — bizim 33 fajımız + REFERANS fajlar (terminaz/TerL marker).
ASLA contig'leri sadece birbirleriyle değil; referanslarla yerleştirir.
Adımlar: QUERY TerL (pharokka) + REF TerL (genbank /translation) -> mafft -> IQ-TREE -> PNG."""
import csv, subprocess, sys, shutil, time, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_pub_advanced import get_terL_seq, latest

ROOT = Path(__file__).resolve().parents[1]
runs = ROOT/"runs"
CACHE = ROOT/"databases"/"ref_genomes_cache"
SBU = Path("/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj")
SAMPLES = ["200225319","2858","21857478","21663260","4188mrsa","19576470psa_001"]
OUT = SBU/"YAYIN_PAKETI_TUM"/"09_Genel_Filogeni_TerL"; OUT.mkdir(parents=True, exist_ok=True)

def cr(cmd, env, t=1800):
    return subprocess.run(["conda","run","-n",env]+cmd, capture_output=True, text=True, timeout=t)

# --- 1) REF terminaz çıkar (genbank /translation, product=terminase, large tercih) ---
EXTRACT = r'''
import sys
from Bio import SeqIO
gb, acc, sp = sys.argv[1], sys.argv[2], sys.argv[3]
best=None; best_large=None
for rec in SeqIO.parse(gb,"genbank"):
    for f in rec.features:
        if f.type=="CDS":
            prod=" ".join(f.qualifiers.get("product",[""])).lower()
            tr=f.qualifiers.get("translation",[""])[0]
            if "terminase" in prod and tr:
                if "large" in prod and not best_large: best_large=tr
                elif not best: best=tr
seq=best_large or best
if seq: print(f">REF_{sp}_{acc}\n{seq}")
'''
Path("/tmp/extract_terL.py").write_text(EXTRACT)

# accession -> tür adı (etiket)
acc2sp={}
for s in SAMPLES:
    run=latest(s); yu=run/"18_Final_Report"/"yayina_uygun"
    if not yu.exists(): continue
    for d in yu.iterdir():
        if not d.is_dir(): continue
        raw=d/f"{d.name}_blast_NCBI_nt.tsv"
        if raw.exists() and raw.stat().st_size>0:
            cc=raw.read_text().splitlines()[0].split("\t"); acc2sp[cc[1]]=cc[8].split(",")[0]
import re
def clean(x): return re.sub(r'[^A-Za-z0-9]+','_',x)[:32].strip('_')

seqs=[]; nref=0
for gb in sorted(CACHE.glob("*.gb")):
    acc=gb.stem; sp=clean(acc2sp.get(acc,acc))
    r=subprocess.run(["conda","run","-n","ali-pharokka","python","/tmp/extract_terL.py",str(gb),acc,sp],
                     capture_output=True,text=True,timeout=120)
    if r.stdout.strip().startswith(">"): seqs.append(r.stdout.strip()); nref+=1

# --- 2) QUERY terminaz (bizim 33) ---
tax={}
tf=SBU/"TUM_TAKSONOMI_RANK.tsv"
if tf.exists():
    rd=list(csv.reader(open(tf),delimiter="\t")); H={h:i for i,h in enumerate(rd[0])}
    for r in rd[1:]: tax[(r[H["ornek"]],r[H["contig"]])]=r[H["genus"]] or "novel"
nq=0
for s in SAMPLES:
    run=latest(s); yu=run/"18_Final_Report"/"yayina_uygun"
    if not yu.exists(): continue
    for d in sorted([x for x in yu.iterdir() if x.is_dir()],key=lambda x:int(x.name.split("_")[1])):
        body=get_terL_seq(d/"pharokka")
        if body:
            g=clean(tax.get((s,d.name),"novel"))
            seqs.append(f">QUERY_{s}_{d.name}_{g}\n{body}"); nq+=1
print(f"QUERY TerL: {nq}  |  REF TerL: {nref}  |  toplam: {len(seqs)}", flush=True)

fa=OUT/"terL_query_plus_ref.faa"; fa.write_text("\n".join(seqs)+"\n")
aln=OUT/"terL_aln.fasta"
aln.write_text(cr(["mafft","--auto",str(fa)],"annotation_prokka",1200).stdout)
print("mafft hizalama bitti.", flush=True)
cr(["iqtree","-s",str(aln),"-m","LG+G","-bb","1000","-nt","4","-redo","-pre",str(OUT/"genel_filogeni")],"ali-iqtree",2400)
print("IQ-TREE bitti.", flush=True)

# --- 3) PNG render (QUERY kırmızı, REF siyah) ---
RENDER=r'''
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo
import sys
tf,out=sys.argv[1],sys.argv[2]
t=Phylo.read(tf,"newick"); t.ladderize()
n=t.count_terminals()
fig=plt.figure(figsize=(13,max(4,n*0.42))); ax=fig.add_subplot(111)
def lab(c):
    if not c.is_terminal(): return ""
    return (c.name or "").replace("QUERY_","▶ ").replace("REF_","").replace("_"," ")
Phylo.draw(t,axes=ax,do_show=False,label_func=lab)
for txt in ax.texts:
    if txt.get_text().startswith("▶"): txt.set_color("crimson"); txt.set_fontweight("bold")
ax.set_title("GENEL FİLOGENİ — TerL (terminaz) ML ağacı | ▶ = bizim fajlar, digerleri REFERANS",fontsize=12,fontweight="bold")
ax.set_xlabel("terminaz protein uzakligi | IQ-TREE LG+G, 1000 bootstrap")
for s in ["top","right","left"]: ax.spines[s].set_visible(False)
ax.get_yaxis().set_visible(False)
plt.tight_layout(); plt.savefig(out,dpi=200,bbox_inches="tight"); print("PNG OK")
'''
Path("/tmp/render_general.py").write_text(RENDER)
tfile=OUT/"genel_filogeni.treefile"
if tfile.exists():
    r=subprocess.run(["conda","run","-n","ali-pharokka","python","/tmp/render_general.py",str(tfile),str(OUT/"genel_filogeni.png")],capture_output=True,text=True,timeout=300)
    print("render:", r.stdout.strip() or r.stderr[-200:], flush=True)
    rr=subprocess.run(["conda","run","-n","ali-pharokka","python","-c",
        f"from Bio import Phylo;t=Phylo.read('{tfile}','newick');Phylo.draw_ascii(t)"],capture_output=True,text=True,timeout=120)
    (OUT/"genel_filogeni_ascii.txt").write_text(rr.stdout or rr.stderr)
print("GENEL FİLOGENİ BİTTİ:", OUT)
