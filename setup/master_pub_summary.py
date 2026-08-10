#!/usr/bin/env python3
"""Çapraz-örnek YAYIN master özeti: tüm örneklerin YAYIN_TABLOSU birleşimi +
33 contig'in çapraz-örnek TerL ağacı (aynı faj farklı örnekte kümelenir)."""
import csv, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_pub_advanced import get_terL_seq, latest

ROOT = Path(__file__).resolve().parents[1]
SBU = Path("/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj")
SAMPLES = ["200225319","2858","21857478","21663260","4188mrsa","19576470psa_001"]
OUT = SBU/"YAYIN_PAKETI_TUM"; OUT.mkdir(exist_ok=True)

def cr(cmd, env, t=1200):
    return subprocess.run(["conda","run","-n",env]+cmd, capture_output=True, text=True, timeout=t)

# 1) birleşik tablo
master=[["ornek","contig","uzunluk","tamlik%","genus","yasam_tarzi","terapotik","intergenomik%","ICTV"]]
for s in SAMPLES:
    run=latest(s)
    if not run: continue
    yt=run/"18_Final_Report"/"YAYIN_PAKETI"/"YAYIN_TABLOSU.tsv"
    if yt.exists():
        for r in list(csv.reader(open(yt),delimiter="\t"))[1:]:
            master.append([s]+r)
with open(OUT/"TUM_YAYIN_TABLOSU.tsv","w",newline="") as fh:
    csv.writer(fh,delimiter="\t").writerows(master)
print(f"TUM_YAYIN_TABLOSU.tsv: {len(master)-1} contig")

# 2) çapraz-örnek TerL ağacı
seqs=[]
# genus map
tax={}
tf=SBU/"TUM_TAKSONOMI_RANK.tsv"
if tf.exists():
    rd=list(csv.reader(open(tf),delimiter="\t")); H={h:i for i,h in enumerate(rd[0])}
    for r in rd[1:]: tax[(r[H["ornek"]],r[H["contig"]])]=r[H["genus"]] or "?"
for s in SAMPLES:
    run=latest(s)
    if not run: continue
    yu=run/"18_Final_Report"/"yayina_uygun"
    if not yu.exists(): continue
    for d in sorted([x for x in yu.iterdir() if x.is_dir()], key=lambda x:int(x.name.split("_")[1])):
        body=get_terL_seq(d/"pharokka")
        if body:
            g=tax.get((s,d.name),"?")
            seqs.append(f">{s}__{d.name}__{g}\n{body}")
print(f"TerL toplandı: {len(seqs)}/33")
if len(seqs)>=3:
    fa=OUT/"terL_all_33.faa"; fa.write_text("\n".join(seqs)+"\n")
    aln=OUT/"terL_all_33_aln.fasta"
    aln.write_text(cr(["mafft","--auto",str(fa)],"annotation_prokka",900).stdout)
    cr(["iqtree","-s",str(aln),"-m","LG","-bb","1000","-nt","3","-redo","-pre",str(OUT/"terL_all_33_tree")],"ali-iqtree",1800)
    tfile=OUT/"terL_all_33_tree.treefile"
    if tfile.exists():
        r=cr(["python","-c",f"from Bio import Phylo;t=Phylo.read('{tfile}','newick');Phylo.draw_ascii(t)"],"ali-pharokka",120)
        (OUT/"terL_all_33_tree_ascii.txt").write_text(r.stdout or r.stderr)
        print("çapraz-örnek TerL ağacı çizildi.")
print("MASTER ÖZET BİTTİ:", OUT)
