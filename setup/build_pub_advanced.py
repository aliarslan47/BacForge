#!/usr/bin/env python3
"""YAYIN MİMARİSİ — gelişmiş kısımlar: 08_VIRIDIC, 09_Phylogeny, 10_Comparative, 12_Publication.
build_publication.py'dan SONRA çalışır. Kullanım: python3 setup/build_pub_advanced.py <ornek>"""
import csv, os, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
runs = ROOT/"runs"
CACHE = ROOT/"databases"/"ref_genomes_cache"; CACHE.mkdir(parents=True, exist_ok=True)
ENVB = "ali-blast"
SBU = Path("/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj")

def latest(s):
    c = sorted([p for p in runs.iterdir() if p.is_dir() and s in p.name], key=lambda p: p.stat().st_mtime, reverse=True)
    return c[0] if c else None
def w(p, rows):
    with open(p,"w",newline="") as fh: csv.writer(fh,delimiter="\t").writerows(rows)
def seqlen(fa): return sum(len(l.strip()) for l in open(fa) if not l.startswith(">"))
def cr(cmd, env, t=900):
    return subprocess.run(["conda","run","-n",env]+cmd, capture_output=True, text=True, timeout=t)

def fetch(acc, fmt="fasta"):
    ext = "fasta" if fmt=="fasta" else "gb"
    out = CACHE/f"{acc}.{ext}"
    if out.exists() and out.stat().st_size>0: return out
    r = subprocess.run(["conda","run","-n",ENVB,"efetch","-db","nuccore","-id",acc,"-format",fmt],
                       capture_output=True,text=True,timeout=300)
    if r.returncode==0 and r.stdout.strip(): out.write_text(r.stdout); return out
    return None

def identical_bp(q, s):
    r = cr(["blastn","-query",str(q),"-subject",str(s),"-outfmt","6 qstart qend pident length","-evalue","1e-5"], ENVB, 600)
    ivs=[]
    for ln in r.stdout.splitlines():
        c=ln.split("\t")
        if len(c)>=4:
            qs,qe,pid,al=int(c[0]),int(c[1]),float(c[2]),int(c[3]); ivs.append((min(qs,qe),max(qs,qe),pid/100))
    if not ivs: return 0
    ivs.sort(); pe=0; tot=0.0
    for s2,e,fr in ivs:
        s2=max(s2,pe+1)
        if e>=s2: tot+=(e-s2+1)*fr; pe=e
    return tot

def folder08_viridic(c, d, ncbi_acc, cfa):
    o=d/"08_Intergenomic_VIRIDIC"; o.mkdir(exist_ok=True)
    if not ncbi_acc:
        w(o/f"{c}_viridic.tsv",[["not","NCBI tür kimliği bekliyor (jumbo) — VIRIDIC sonra"]]); return None
    ref=fetch(ncbi_acc,"fasta")
    if not ref:
        w(o/f"{c}_viridic.tsv",[["not",f"referans {ncbi_acc} çekilemedi"]]); return None
    lc,lr=seqlen(cfa),seqlen(ref)
    sim=round(100*(identical_bp(cfa,ref)+identical_bp(ref,cfa))/(lc+lr),2)
    verd=("🔵 bilinen tür" if sim>=95 else "🟢 YENİ TÜR (aynı cins)" if sim>=70 else "🟠 YENİ CİNS" if sim>=50 else "🔴 ÇOK IRAKSAK (yeni cins/üzeri)")
    w(o/f"{c}_viridic.tsv",[["karsilastirma","intergenomik_benzerlik%","ICTV_karar"],
        [f"{c} vs {ncbi_acc}", sim, verd],["ICTV_sinir","tür>=95, cins>=70",""]])
    return sim

def folder10_comparative(c, d, ncbi_acc):
    o=d/"10_Comparative"; o.mkdir(exist_ok=True)
    gbk=d/"04_Annotation"/f"{c}.gbk"
    if not gbk.exists() or not ncbi_acc:
        w(o/"NOT.txt",[["clinker icin contig GBK + referans gerekli; jumbo/eksik ise sonra"]]); return
    refgb=fetch(ncbi_acc,"gb")
    if not refgb:
        (o/"NOT.txt").write_text(f"referans genbank {ncbi_acc} cekilemedi\n"); return
    # clinker iki gbk ister; geçici klasör
    tmp=o/"_in"; tmp.mkdir(exist_ok=True)
    import shutil
    shutil.copy2(gbk, tmp/f"{c}.gbk")
    refout=tmp/f"{ncbi_acc}.gbk"; shutil.copy2(refgb, refout)
    r=cr(["clinker", str(tmp/f"{c}.gbk"), str(refout), "-p", str(o/f"{c}_synteny.html"), "-o", str(o/f"{c}_alignments.csv")], "ali-clinker", 600)
    (o/"clinker.log").write_text(r.stdout+"\n"+r.stderr)
    shutil.rmtree(tmp, ignore_errors=True)

def get_terL_seq(pdir):
    """terL.faa varsa onu; yoksa anotasyondan terminaz(large) CDS'ini phanotate.faa'dan çek (tek-marker bütünlüğü)."""
    tl = pdir/"terL.faa"
    if tl.exists() and tl.stat().st_size > 0:
        txt = tl.read_text().strip().split("\n")
        return "".join(l for l in txt[1:] if not l.startswith(">"))
    merged = pdir/"pharokka_cds_final_merged_output.tsv"; faa = pdir/"phanotate.faa"
    if not (merged.exists() and faa.exists()): return None
    large=[]; anyt=[]
    for line in merged.read_text().splitlines():
        low=line.lower()
        if "terminase" in low:
            gid=line.split("\t")[0]
            (large if "large" in low else anyt).append(gid)
    gid = large[0] if large else (anyt[0] if anyt else None)
    if not gid: return None
    seq=[]; grab=False
    for line in faa.read_text().splitlines():
        if line.startswith(">"): grab = (line[1:].split()[0] == gid)
        elif grab: seq.append(line.strip())
    return "".join(seq) or None

def phylogeny(run, contigs, info):
    """Örnek düzeyi: tüm contig'lerin terminaz (TerL) ağacı."""
    pkg = run/"18_Final_Report"/"YAYIN_PAKETI"; pkg.mkdir(exist_ok=True)
    phy = pkg/"09_Phylogeny_TerL"; phy.mkdir(exist_ok=True)
    concat = phy/"terL_all.faa"; seqs=[]
    for c in contigs:
        body = get_terL_seq(run/"18_Final_Report"/"yayina_uygun"/c/"pharokka")
        if body:
            header=f">{c}_{info.get(c,{}).get('genus','?')}_{info.get(c,{}).get('life','?')}"
            seqs.append(header+"\n"+body)
    if len(seqs)<3:
        (phy/"NOT.txt").write_text(f"Yeterli TerL yok ({len(seqs)}) — ağaç atlandı\n"); return phy, len(seqs)
    concat.write_text("\n".join(seqs)+"\n")
    aln=phy/"terL_aln.fasta"
    r=cr(["mafft","--auto",str(concat)], "annotation_prokka", 600)
    aln.write_text(r.stdout)
    # IQ-TREE
    cr(["iqtree","-s",str(aln),"-m","LG","-bb","1000","-nt","2","-redo","-pre",str(phy/"terL_tree")], "ali-iqtree", 900)
    # ascii ağaç (Bio.Phylo)
    tf=phy/"terL_tree.treefile"
    if tf.exists():
        try:
            r2=cr(["python","-c",
                f"from Bio import Phylo;import io;t=Phylo.read('{tf}','newick');Phylo.draw_ascii(t)"], "ali-pharokka", 120)
            (phy/"terL_tree_ascii.txt").write_text(r2.stdout or r2.stderr)
        except Exception: pass
    return phy, len(seqs)

def main(sample):
    run=latest(sample)
    yu=run/"18_Final_Report"/"yayina_uygun"
    contigs=[d.name for d in sorted([x for x in yu.iterdir() if x.is_dir()], key=lambda x:int(x.name.split("_")[1]))]
    # taksonomi/genus + lifestyle bilgisi
    info={}
    taxf=SBU/"TUM_TAKSONOMI_RANK.tsv"
    tx={}
    if taxf.exists():
        rd=list(csv.reader(open(taxf),delimiter="\t")); H={h:i for i,h in enumerate(rd[0])}
        for r in rd[1:]:
            if r[H["ornek"]]==sample: tx[r[H["contig"]]]=(r[H["genus"]] or "?", r[H["en_yakin_tur(species)"]])
    print(f"=== {sample}: 08 VIRIDIC + 10 clinker (per contig) ===", flush=True)
    pub_rows=[["contig","uzunluk","tamlik%","genus","yasam_tarzi","terapotik","intergenomik%","ICTV"]]
    for c in contigs:
        d=yu/c
        # NCBI accession
        acc=""; raw=d/f"{c}_blast_NCBI_nt.tsv"
        if raw.exists() and raw.stat().st_size>0: acc=raw.read_text().splitlines()[0].split("\t")[1]
        cfa=d/"01_Genome"/f"{c}.fasta"
        sim=folder08_viridic(c,d,acc,cfa)
        folder10_comparative(c,d,acc)
        # lifestyle
        life=""
        lf=d/"06_Lifestyle"/f"{c}_lifestyle.tsv"
        if lf.exists():
            kv={r[0]:r[1] for r in csv.reader(open(lf),delimiter="\t") if len(r)>=2}; life=kv.get("yasam_tarzi","")
        sv=""
        svf=d/"07_Safety_AMR_Virulence"/f"{c}_safety_karar.tsv"
        if svf.exists():
            kv={r[0]:r[1] for r in csv.reader(open(svf),delimiter="\t") if len(r)>=2}; sv=kv.get("TERAPOTIK_KARAR","")
        st={r[0]:r[1] for r in csv.reader(open(d/f"{c}_istatistik.tsv"),delimiter="\t") if len(r)>=2}
        ictv = ("🔵tür" if (sim and sim>=95) else "🟢yeni tür" if (sim and sim>=70) else "🟠yeni cins" if (sim and sim>=50) else "🔴çok ıraksak" if sim else "—")
        info[c]={"genus":tx.get(c,("?",""))[0],"life":"L" if "Lytic" in life else "T"}
        pub_rows.append([c,st.get("uzunluk_bp",""),st.get("CheckV_tamlik%",""),tx.get(c,("?",""))[0],life,sv,sim if sim else "—",ictv])
        print(f"  {c}: VIRIDIC={sim}  {ictv}", flush=True)
    # 09 filogeni (örnek düzeyi)
    print("=== 09 Phylogeny (TerL ağacı) ===", flush=True)
    phy,n=phylogeny(run, contigs, info)
    print(f"  TerL ağacı: {n} dizi", flush=True)
    # 12 publication paketi
    pkg=run/"18_Final_Report"/"YAYIN_PAKETI"; pkg.mkdir(exist_ok=True)
    w(pkg/"YAYIN_TABLOSU.tsv", pub_rows)
    # her contig 12_Publication: konsolide + announcement
    for c in contigs:
        d=yu/c; o=d/"12_Publication"; o.mkdir(exist_ok=True)
        row=[r for r in pub_rows[1:] if r[0]==c]
        if row:
            r=row[0]
            w(o/f"{c}_ozet.tsv",[["alan","deger"],["contig",c],["uzunluk_bp",r[1]],["tamlik%",r[2]],
                ["genus",r[3]],["yasam_tarzi",r[4]],["terapotik",r[5]],["intergenomik%",r[6]],["ICTV",r[7]]])
            # GenBank submission dosyası referansı
            gbk=d/"04_Annotation"/f"{c}.gbk"
            ann=(f"Bacteriophage {c} ({r[1]} bp, %{r[2]} tam, {r[3]} cinsi), {sample} klinik örneğinden "
                 f"izole edildi. Genom {('lytic' if 'Lytic' in r[4] else 'ılımlı')} yaşam tarzı gösterir; "
                 f"AMR ve virülans geni saptanmadı. En yakın akrabaya intergenomik benzerlik %{r[6]} "
                 f"({r[7].strip('🔵🟢🟠🔴')}). GenBank: {gbk.name}.")
            (o/f"{c}_announcement_taslak.txt").write_text(ann+"\n")
    # --- son faz: her contig'e 09_Phylogeny + analiz/ kopya + tarih tazele ---
    import shutil
    phy=pkg/"09_Phylogeny_TerL"
    for c in contigs:
        d=yu/c; o=d/"09_Phylogeny"; o.mkdir(exist_ok=True)
        tl=d/"pharokka"/"terL.faa"
        if tl.exists() and tl.stat().st_size>0: shutil.copy2(tl,o/f"{c}_terL.faa")
        for f in ["terL_tree.treefile","terL_tree_ascii.txt","terL_tree.contree"]:
            if (phy/f).exists(): shutil.copy2(phy/f,o/f)
        (o/"NOT.txt").write_text("Bu contig'in terminaz (TerL) dizisi + ornek-duzeyi TerL agaci.\n")
    dest=SBU/sample/"analiz"
    if dest.exists():
        if (dest/"18_Final_Report").exists(): shutil.rmtree(dest/"18_Final_Report")
        shutil.copytree(run/"18_Final_Report", dest/"18_Final_Report")
        now=time.time()
        for p in [dest/"18_Final_Report", *(dest/"18_Final_Report").rglob("*")]:
            try: os.utime(p,(now,now))
            except OSError: pass
        print(f"  analiz/ kopyalandı + tarih tazelendi: {dest}")
    print("\n=== GELİŞMİŞ KISIM BİTTİ ===")
    for r in pub_rows[1:]:
        print(f"  {r[0]:<12}{str(r[1]):>8}bp genus={r[3]:<16} {r[7]:<14} {r[5][:38]}")

if __name__=="__main__":
    main(sys.argv[1] if len(sys.argv)>1 else "4188mrsa")
