#!/usr/bin/env python3
"""Her yayına uygun faj için RANKLI taksonomi (NCBI Taxonomy: realm→...→genus→species)
ekle. Kaynak: en yakın NCBI core_nt türünün taxid'i -> ranklı lineage. Çıktı:
- her dossier'a <c>_taksonomi.tsv
- her örneğe 08_Taxonomy/taksonomi_rank.tsv
- sbu-faj/TUM_TAKSONOMI_RANK.tsv
+ analiz/ kopyala, tarih tazele."""
import csv, subprocess, time, os, shutil
from pathlib import Path

runs = Path("runs")
SBU = Path("/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj")
SAMPLES = ["200225319","2858","21857478","21663260","4188mrsa","19576470psa_001"]
ENV = "ali-blast"
RANKS = ["realm","kingdom","phylum","class","order","family","subfamily","genus"]

def latest(s):
    c=sorted([p for p in runs.iterdir() if p.is_dir() and s in p.name],key=lambda p:p.stat().st_mtime,reverse=True)
    return c[0] if c else None

# novelty verdict (proxy+ICTV) ve VIRIDIC
nov={}
f=SBU/"YENI_TUR_ADAYLARI.tsv"
if f.exists():
    rd=list(csv.reader(open(f),delimiter="\t")); hd=rd[0]; idx={h:i for i,h in enumerate(hd)}
    for r in rd[1:]:
        nov[(r[idx["ornek"]],r[idx["contig"]])]=(r[idx["intergenomik_proxy%"]],r[idx["ICTV_karar"]])
vir={}
f=SBU/"YENI_TUR_VIRIDIC_DOGRULAMA.tsv"
if f.exists():
    for r in list(csv.reader(open(f),delimiter="\t"))[1:]:
        if len(r)>=7 and r[6] not in ("-","efetch yok"): vir[(r[0],r[1])]=r[6]

_cache={}
def ncbi_taxo(acc):
    if acc in _cache: return _cache[acc]
    res={"species":"","ranks":{}}
    try:
        taxid=subprocess.run(["conda","run","-n",ENV,"bash","-c",
            f"efetch -db nuccore -id {acc} -format docsum 2>/dev/null | xtract -pattern DocumentSummary -element TaxId"],
            capture_output=True,text=True,timeout=180).stdout.strip().split("\n")[0].strip()
        if taxid:
            out=subprocess.run(["conda","run","-n",ENV,"bash","-c",
                f'efetch -db taxonomy -id {taxid} -format xml 2>/dev/null | xtract -pattern Taxon -element ScientificName -block "*/Taxon" -sep "|" -element Rank,ScientificName'],
                capture_output=True,text=True,timeout=180).stdout.strip()
            if out:
                parts=out.split("\t")
                res["species"]=parts[0]
                for tok in parts[1:]:
                    if "|" in tok:
                        rk,nm=tok.split("|",1)
                        if rk in RANKS: res["ranks"][rk]=nm
    except Exception as e:
        print(f"  !! taxo hata {acc}: {e}")
    _cache[acc]=res; return res

master=[["ornek","contig","host","realm","kingdom","phylum","class","order","family","subfamily","genus","en_yakin_tur(species)","intergenomik%","taksonomik_durum"]]
for s in SAMPLES:
    run=latest(s)
    if not run: continue
    yu=run/"18_Final_Report"/"yayina_uygun"
    if not yu.exists(): continue
    # geNomad lineage yedeği
    gtax=run/"08_Taxonomy"/"genomad"/"contigs.filtered_annotate"/"contigs.filtered_taxonomy.tsv"
    glin={}
    if gtax.exists():
        for r in csv.DictReader(open(gtax),delimiter="\t"): glin[r["seq_name"]]=r.get("lineage","")
    sample_rows=[master[0]]
    for d in sorted([x for x in yu.iterdir() if x.is_dir()],key=lambda x:int(x.name.split("_")[1])):
        c=d.name
        raw=d/f"{c}_blast_NCBI_nt.tsv"
        acc=sp=""
        if raw.exists() and raw.stat().st_size>0:
            rr=raw.read_text().splitlines()[0].split("\t"); acc=rr[1]; sp=rr[8]
        # INPHARED host
        host=""
        mf=d/"pharokka"/"pharokka_top_hits_mash_inphared.tsv"
        if mf.exists():
            md=list(csv.reader(open(mf),delimiter="\t"))
            if len(md)>1: host=dict(zip(md[0],md[1])).get("Host","")
        ranks={r:"" for r in RANKS}; species=""
        if acc:
            t=ncbi_taxo(acc); species=t["species"]
            for r in RANKS: ranks[r]=t["ranks"].get(r,"")
        # geNomad backbone ile boşları doldur (Realm..Class + Family)
        if glin.get(c):
            g=glin[c].replace("[NCBI:core_nt]","").split(";")
            gmap={"Duplodnaviria":"realm","Heunggongvirae":"kingdom","Uroviricota":"phylum","Caudoviricetes":"class"}
            for tok in g:
                tok=tok.strip()
                if tok in gmap and not ranks[gmap[tok]]: ranks[gmap[tok]]=tok
                if tok.endswith("viridae") and not ranks["family"]: ranks["family"]=tok
        viri=vir.get((s,c),""); prox,ictv=nov.get((s,c),("",""))
        simshow=viri if viri else prox
        # taksonomik durum
        if not acc: durum="tür kimliği bekliyor (jumbo)"
        elif viri and float(viri)<70: durum=f"YENİ CİNS/üzeri — en yakın cins {ranks.get('genus') or '?'} (%{viri})"
        elif (viri and float(viri)<95) or (prox and float(prox)<95): durum=f"YENİ TÜR — {ranks.get('genus') or '?'} cinsinde (en yakın {sp.split(',')[0]})"
        else: durum=f"bilinen türe çok yakın ({sp.split(',')[0]})"
        row=[s,c,host,ranks["realm"],ranks["kingdom"],ranks["phylum"],ranks["class"],
             ranks["order"],ranks["family"],ranks["subfamily"],ranks["genus"],
             species or sp.split(",")[0], simshow, durum]
        master.append(row); sample_rows.append(row)
        # per-contig dossier dosyası
        kv=[["rank","atama"]]
        for label,rk in [("Realm","realm"),("Kingdom","kingdom"),("Phylum","phylum"),("Class","class"),
                         ("Order","order"),("Family","family"),("Subfamily","subfamily"),("Genus","genus")]:
            kv.append([label, ranks[rk] or "(sınıflanmamış)"])
        kv += [["Species (en yakın)", species or sp.split(",")[0] or "—"],
               ["Konak (host)", host or "—"],["Baltimore", "Group I (dsDNA)"],
               ["intergenomik benzerlik%", simshow or "—"],["taksonomik durum", durum]]
        with open(d/f"{c}_taksonomi.tsv","w",newline="") as fh:
            csv.writer(fh,delimiter="\t").writerows(kv)
    # per-örnek taksonomi
    with open(run/"08_Taxonomy"/"taksonomi_rank.tsv","w",newline="") as fh:
        csv.writer(fh,delimiter="\t").writerows(sample_rows)
    # analiz kopya + tarih tazele
    dest=SBU/s/"analiz"
    if dest.exists():
        (dest/"08_Taxonomy").mkdir(exist_ok=True)
        shutil.copy2(run/"08_Taxonomy"/"taksonomi_rank.tsv", dest/"08_Taxonomy"/"taksonomi_rank.tsv")
        # dossier per-contig taksonomi dosyalarını da kopyala (18_Final_Report zaten kopyalanıyor ama tazele)
        for d in (run/"18_Final_Report"/"yayina_uygun").iterdir():
            if d.is_dir():
                src=d/f"{d.name}_taksonomi.tsv"
                if src.exists():
                    dd=dest/"18_Final_Report"/"yayina_uygun"/d.name
                    dd.mkdir(parents=True,exist_ok=True)
                    shutil.copy2(src, dd/src.name)
        now=time.time()
        for p in [dest/"08_Taxonomy", *(dest/"08_Taxonomy").rglob("*"), dest/"18_Final_Report", *(dest/"18_Final_Report").rglob("*")]:
            try: os.utime(p,(now,now))
            except OSError: pass
    print(f"  {s}: {len(sample_rows)-1} contig taksonomisi yazıldı/kopyalandı")

with open(SBU/"TUM_TAKSONOMI_RANK.tsv","w",newline="") as fh:
    csv.writer(fh,delimiter="\t").writerows(master)
print(f"\nTUM_TAKSONOMI_RANK.tsv yazıldı ({len(master)-1} contig).")
# özet yazdır
for r in master[1:]:
    print(f"  {r[0]}/{r[1]:<12} {r[6] or '?':<14} {r[8] or '(aile?)':<18} {r[10] or '(cins?)':<16} {r[13]}")
