#!/usr/bin/env python3
"""YENİ TÜR ADAYI analizi — 33 yayına uygun faj.
Her contig için: NCBI core_nt en yakın tür + intergenomik benzerlik PROXY (pident*qcov/100,
VIRIDIC-uyumlu eş-uzunlukta) + INPHARED mash en yakın (NCBI'dan bağımsız) + pharokka genom
istatistikleri (uzunluk/GC/kodlama yoğunluğu/CDS/bilinmeyen%/tRNA/integrase) + ICTV kararı.
ICTV: tür sınırı %95, cins sınırı %70 (nükleotid intergenomik benzerlik)."""
import csv
from pathlib import Path

runs = Path("runs")
SBU = Path("/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj")
SAMPLES = ["200225319", "2858", "21857478", "21663260", "4188mrsa", "19576470psa_001"]

def latest(s):
    c = sorted([p for p in runs.iterdir() if p.is_dir() and s in p.name],
               key=lambda p: p.stat().st_mtime, reverse=True)
    return c[0] if c else None

def read_kv(f, sep="\t"):
    d = {}
    if f.exists():
        for r in csv.reader(open(f), delimiter=sep):
            if len(r) >= 2: d[r[0]] = r[1]
    return d

def ictv(prox, qcov):
    if prox == "" : return "—"
    if prox >= 95: return "🔵 bilinen tür"
    if prox >= 70: return "🟢 YENİ TÜR adayı (aynı cins)"
    if prox >= 50: return "🟠 YENİ CİNS adayı"
    return "🔴 ÇOK IRAKSAK — yeni cins/alt-aile (en yüksek değer)"

rows = [["ornek","contig","uzunluk","GC%","kodlama%","CDS","bilinmeyen%","tRNA","integrase(temperate?)",
         "NCBI_en_yakin_tur","NCBI_id%","NCBI_qcov%","intergenomik_proxy%","INPHARED_mash_en_yakin","mash_dist","ICTV_karar"]]
for s in SAMPLES:
    run = latest(s)
    if not run: continue
    yu = run/"18_Final_Report"/"yayina_uygun"
    if not yu.exists(): continue
    for d in sorted([x for x in yu.iterdir() if x.is_dir()], key=lambda x: int(x.name.split("_")[1])):
        c = d.name
        # NCBI
        raw = d/f"{c}_blast_NCBI_nt.tsv"
        sp = pid = qc = ""
        if raw.exists() and raw.stat().st_size > 0:
            r = raw.read_text().splitlines()[0].split("\t")
            sp, pid, qc = r[8], float(r[2]), float(r[4])
        prox = round(pid*qc/100, 1) if pid != "" and qc != "" else ""
        # pharokka istatistik
        P = d/"pharokka"
        lg = {}
        f = P/"pharokka_length_gc_cds_density.tsv"
        if f.exists():
            rd = list(csv.reader(open(f), delimiter="\t"))
            if len(rd) > 1: lg = dict(zip(rd[0], rd[1]))
        length = lg.get("length",""); gc = lg.get("gc_perc",""); dens = lg.get("cds_coding_density","")
        try: gc = f"{float(gc)*100:.1f}" if gc and float(gc) < 1 else gc
        except: pass
        cf = {}
        ff = P/"pharokka_cds_functions.tsv"
        if ff.exists():
            for r in csv.reader(open(ff), delimiter="\t"):
                if len(r) >= 2 and r[0] != "Description": cf[r[0]] = r[1]
        cds = cf.get("CDS",""); unk = cf.get("unknown function",""); trna = cf.get("tRNAs","")
        integ = cf.get("integration and excision","")
        unk_pct = ""
        try: unk_pct = f"{int(unk)/int(cds)*100:.0f}" if cds and int(cds) else ""
        except: pass
        temperate = ("evet(integrase var)" if integ not in ("","0") else "hayır(lytic eğilim)") if integ != "" else ""
        # INPHARED mash
        mash_hit = mash_d = ""
        mf = P/"pharokka_top_hits_mash_inphared.tsv"
        if mf.exists():
            rd = list(csv.reader(open(mf), delimiter="\t"))
            if len(rd) > 1:
                h = dict(zip(rd[0], rd[1]))
                desc = h.get("Description","")
                if "no_inphared" in str(desc) or desc == "":
                    mash_hit = "YOK (INPHARED'da akraba yok!)"
                else:
                    mash_hit = f"{desc} [{h.get('Genus','')}/{h.get('Family','')}]"
                    mash_d = h.get("mash_distance","")
        verdict = ictv(prox, qc) if prox != "" else "— (NCBI bekliyor)"
        rows.append([s,c,length,gc,dens,cds,unk_pct,trna,temperate,sp,pid,qc,prox,mash_hit,mash_d,verdict])

# yaz
out = SBU/"YENI_TUR_ADAYLARI.tsv"
with open(out, "w", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(rows)

# sırala: intergenomik proxy artan (en novel önce), NCBI bekleyenler sonda
def keyf(r):
    p = r[12]
    return (p if isinstance(p, (int, float)) else 9999)
body = sorted(rows[1:], key=keyf)
print(f"YENI_TUR_ADAYLARI.tsv yazıldı ({len(rows)-1} contig)\n")
print(f"{'ornek/contig':<26}{'uzunluk':>8}{'proxy%':>8}  {'tRNA/integrase':<22}{'ICTV karar'}")
print("-"*110)
for r in body:
    name = f"{r[0]}/{r[1]}"
    prox = r[12] if r[12] != "" else "—"
    ti = f"t{r[7] or '?'}/{('temp' if 'evet' in str(r[8]) else 'lytic' if 'hayır' in str(r[8]) else '?')}"
    print(f"{name:<26}{str(r[2]):>8}{str(prox):>8}  {ti:<22}{r[15]}   [{r[9][:34]} {r[10]}%/{r[11]}%]")
