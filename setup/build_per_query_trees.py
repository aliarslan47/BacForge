#!/usr/bin/env python3
"""PER-QUERY ICTV AĞAÇLARI — her query'i KENDİ cinsinin NCBI türleri + komşu cins (aile)
üyeleriyle yerleştirir. Query'i ASLA sadece birbiriyle değil; daima referans taksonlarla.
ICTV-bağlı: referanslar NCBI/ICTV tanımlı genus+family üyeleri. Yöntem: mash genom mesafesi NJ.
Çıktı: her contig'in 09_Phylogeny/<contig>_ICTV_agac.png"""
import csv, subprocess, json, re, time, os, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
runs = ROOT/"runs"
RCACHE = ROOT/"databases"/"genus_ref_cache"; RCACHE.mkdir(parents=True, exist_ok=True)
SBU = Path("/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj")
SAMPLES = ["200225319","2858","21857478","21663260","4188mrsa","19576470psa_001"]
ENVB, ENVP = "ali-blast", "ali-pharokka"
STATUS = ROOT/"setup"/"per_query.STATUS"
GENUS_CAP, FAM_CAP = 10, 6

def latest(s):
    c=sorted([p for p in runs.iterdir() if p.is_dir() and s in p.name],key=lambda p:p.stat().st_mtime,reverse=True); return c[0]
def clean(x): return re.sub(r'[^A-Za-z0-9]+','_',x)[:40].strip('_')

def esearch_accs(term, cap):
    r=subprocess.run(["conda","run","-n",ENVB,"bash","-c",
        f"esearch -db nuccore -query '{term}' 2>/dev/null | efetch -format acc 2>/dev/null | head -{cap}"],
        capture_output=True,text=True,timeout=180)
    return [a.strip() for a in r.stdout.splitlines() if a.strip()]

def fetch_fasta(acc):
    out=RCACHE/f"{acc}.fasta"
    if out.exists() and out.stat().st_size>0: return out
    r=subprocess.run(["conda","run","-n",ENVB,"efetch","-db","nuccore","-id",acc,"-format","fasta"],
                     capture_output=True,text=True,timeout=200)
    if r.returncode==0 and r.stdout.startswith(">"): out.write_text(r.stdout); return out
    return None

_mem={}
def genus_members(term, cap):
    if term in _mem: return _mem[term]
    accs=esearch_accs(term, cap*2)
    res=[]
    for a in accs:
        f=fetch_fasta(a)
        if f:
            title=open(f).readline()[1:].strip()
            res.append((a, title, f))
        if len(res)>=cap: break
    _mem[term]=res; return res

def hb(done, total, cur):
    STATUS.write_text(f"status=running query_agac={done}/{total} su_an={cur} hb={time.strftime('%H:%M:%S')}\n")


def main():
    # taksonomi: query -> (genus, family, species, sample)
    tax={}
    tf=SBU/"TUM_TAKSONOMI_RANK.tsv"
    rd=list(csv.reader(open(tf),delimiter="\t")); H={h:i for i,h in enumerate(rd[0])}
    for r in rd[1:]:
        tax[(r[H["ornek"]],r[H["contig"]])]=(r[H["genus"]], r[H["family"]], r[H["en_yakin_tur(species)"]])

    # query fastas + genus grupları
    queries=[]  # (sample, contig, genus, family, fasta)
    for s in SAMPLES:
        run=latest(s); yu=run/"18_Final_Report"/"yayina_uygun"
        for d in sorted([x for x in yu.iterdir() if x.is_dir()],key=lambda x:int(x.name.split("_")[1])):
            g,fam,sp=tax.get((s,d.name),("","",""))
            fa=d/"01_Genome"/f"{d.name}.fasta"
            if not fa.exists(): fa=d/f"{d.name}.fasta"
            queries.append((s,d.name,g,fam,fa))

    # her cins/aile için referansları ÖNCEDEN topla
    genera=sorted({q[2] for q in queries if q[2] and q[2] not in ("(cins?)","")})
    fams=sorted({q[3] for q in queries if q[3] and q[3] not in ("(aile?)","")})
    print(f"{len(queries)} query | {len(genera)} cins | {len(fams)} aile referans toplanacak", flush=True)

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
    ax.set_xlabel("mash genom mesafesi (k-mer) NJ | ★=bu query, ▶=bizim ayni cins, digerleri ICTV referans")
    for sp in ["top","right","left"]: ax.spines[sp].set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.tight_layout(); plt.savefig(sys.argv[2],dpi=190,bbox_inches="tight"); print("OK")
    '''
    Path("/tmp/nj_per_query.py").write_text(NJ)

    def build_tree(genomes, focal_label, title, out_png, out_tree, workdir):
        """genomes: list of (label, fasta_path)"""
        workdir.mkdir(parents=True, exist_ok=True)
        fnas=[]
        for lab,fa in genomes:
            p=workdir/f"{lab}.fna"
            seq="".join(l for l in Path(fa).read_text().splitlines() if not l.startswith(">"))
            p.write_text(f">{lab}\n{seq}\n"); fnas.append(p)
        subprocess.run(["conda","run","-n",ENVP,"mash","sketch","-s","20000","-k","17","-o",str(workdir/"all")]+[str(f) for f in fnas],
                       capture_output=True,timeout=600)
        d=subprocess.run(["conda","run","-n",ENVP,"mash","dist",str(workdir/"all.msh"),str(workdir/"all.msh")],
                         capture_output=True,text=True,timeout=600)
        labels=[f.stem for f in fnas]; idx={l:i for i,l in enumerate(labels)}; n=len(labels)
        M=[[0.0]*n for _ in range(n)]
        for ln in d.stdout.splitlines():
            c=ln.split("\t")
            if len(c)>=3:
                a=Path(c[0]).stem; b=Path(c[1]).stem
                if a in idx and b in idx: M[idx[a]][idx[b]]=float(c[2])
        lower=[[M[i][j] for j in range(i+1)] for i in range(n)]
        jd=workdir/"_d.json"; jd.write_text(json.dumps({"labels":labels,"lower":lower,"focal":focal_label,"title":title}))
        subprocess.run(["conda","run","-n",ENVP,"python","/tmp/nj_per_query.py",str(jd),str(out_png)],capture_output=True,timeout=300)
        shutil.rmtree(workdir, ignore_errors=True)

    total=len(queries)
    for i,(s,c,g,fam,fa) in enumerate(queries,1):
        hb(i,total,f"{s}/{c}")
        focal=f"QUERY_{s}_{c}_{clean(g or 'novel')}"
        genomes=[(focal, fa)]
        # aynı cins kardeş query'lerimiz
        for (s2,c2,g2,fam2,fa2) in queries:
            if (s2,c2)!=(s,c) and g and g2==g:
                genomes.append((f"QUERY_{s2}_{c2}_{clean(g)}", fa2))
        # cins referans üyeleri
        refset=[]
        if g and g not in ("(cins?)",""):
            refset=genus_members(f'"{g}"[Organism] AND complete genome AND viruses[filter]', GENUS_CAP)
        # komşu cins (aile) üyeleri
        famrefs=[]
        if fam and fam not in ("(aile?)",""):
            fammem=genus_members(f'"{fam}"[Organism] AND complete genome AND viruses[filter]', FAM_CAP+GENUS_CAP)
            seen={a for a,_,_ in refset}
            for a,t,f in fammem:
                if a not in seen and (not g or g.lower() not in t.lower()):
                    famrefs.append((a,t,f))
                if len(famrefs)>=FAM_CAP: break
        for a,t,f in refset+famrefs:
            genomes.append((f"REF_{clean(t)}", f))
        if len(genomes)<3:
            continue
        d=latest(s)/"18_Final_Report"/"yayina_uygun"/c/"09_Phylogeny"; d.mkdir(parents=True,exist_ok=True)
        title=f"{s}/{c} ICTV yerlesim — cins {g or 'novel'} / aile {fam or '?'}"
        build_tree(genomes, focal, title, d/f"{c}_ICTV_agac.png", None, runs/"_pq_tmp"/f"{s}_{c}")
        print(f"  [{i}/{total}] {s}/{c}: {len(genomes)} genom (cins {g or 'novel'})", flush=True)

    STATUS.write_text(f"status=done {time.strftime('%H:%M:%S')}\n")
    print("\n=== PER-QUERY ICTV AĞAÇLARI BİTTİ ===")


if __name__=="__main__":
    main()