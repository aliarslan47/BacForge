"""M14 -- Genomic Context & Closest-N Comparison (HEDEFLI clinker)
Spec §18: clinker tum-genom icin DEGIL, secili homolog bolgeler icin. Bu modul AMR (M08) ve
virulans (M09) genlerini hedef alir; sorgu ve en-yakin referanslarin Bakta-GBK'lerinde o genin
+-W komsulugunu cikarir; lokus basina clinker -> sinteni HTML. (Hizli + biyolojik olarak anlamli.)
KATI: en az bir lokus icin gercek clinker HTML uretilmeden PASS YOK. Hedef/referans yoksa NOT_APPLICABLE.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from .base import Module
from .. import util

WINDOW = 15000     # gen etrafinda +-15 kb
MAX_TARGETS = 6    # en fazla lokus (AMR + virulans karisik)
N_REFS = 5         # query + en yakin 5 tur ile karsilastir

# ali-comparative (Biopython) icinde calisir: hedef genin komsulugunu her genomdan kes.
SLICER = r'''
import sys, json
from Bio import SeqIO
q_gbk, refs_json, targets_json, window, outdir = sys.argv[1:6]
window = int(window); outdir = __import__("pathlib").Path(outdir)
refs = json.load(open(refs_json))          # [{"label":..,"gbk":..}]
targets = json.load(open(targets_json))    # ["blaOXA-66", ...]

def load(gbk):
    return list(SeqIO.parse(gbk, "genbank"))

def find_gene(records, gene):
    g = gene.lower()
    for rec in records:
        for f in rec.features:
            if f.type != "CDS":
                continue
            names = [x.lower() for x in f.qualifiers.get("gene", [])]
            prods = [x.lower() for x in f.qualifiers.get("product", [])]
            if any(g == n or g in n for n in names) or any(g in p for p in prods):
                return rec, int(f.location.start), int(f.location.end)
    return None

def slab(rec, s, e, label):
    a = max(0, s - window); b = min(len(rec.seq), e + window)
    sub = rec[a:b]
    sub.id = label[:16]; sub.name = label[:16]
    sub.annotations["molecule_type"] = "DNA"
    return sub

qrecs = load(q_gbk)
refrecs = {r["label"]: load(r["gbk"]) for r in refs}
summary = []
for gene in targets:
    hit = find_gene(qrecs, gene)
    if not hit:
        summary.append({"gene": gene, "found_in_query": False}); continue
    rec, s, e = hit
    tdir = outdir / gene.replace("/", "_").replace("'", "")
    tdir.mkdir(parents=True, exist_ok=True)
    SeqIO.write(slab(rec, s, e, "QUERY"), tdir / "QUERY.gbk", "genbank")
    have = ["QUERY"]
    for label, rr in refrecs.items():
        rh = find_gene(rr, gene)
        if rh:
            rrec, rs, re_ = rh
            SeqIO.write(slab(rrec, rs, re_, label), tdir / (label + ".gbk"), "genbank")
            have.append(label)
    summary.append({"gene": gene, "found_in_query": True, "genomes": have, "dir": str(tdir)})
json.dump(summary, open(outdir / "targets_summary.json", "w"), indent=2)
print("SLICE_DONE", sum(1 for x in summary if x.get("found_in_query")))
'''


class GenomicContextModule(Module):
    number = "14"
    name = "genomic_context"
    folder = "M14_GENOMIC_CONTEXT"
    enabled_key = "clinker"

    def inputs(self):
        return [self.ctx.run_dir / "M06_GENOME_ANNOTATION" / "annotation.gbk"]

    def outputs(self):
        return [self.out_dir / "M14_summary.json"]

    def _targets(self) -> list[str]:
        genes = []
        amr = self.ctx.run_dir / "M08_AMR" / "amr.json"
        if amr.exists():
            try:
                for g in json.load(open(amr)).get("amr_genes", []):
                    gs = g.get("gene_symbol")
                    if gs and gs not in genes:
                        genes.append(gs)
            except Exception:
                pass
        amr_n = len(genes)
        vir = self.ctx.run_dir / "M09_VIRULENCE" / "virulence.json"
        if vir.exists():
            try:
                for g in json.load(open(vir)).get("virulence_genes", []):
                    gs = g.get("gene_symbol")
                    if gs and gs not in genes:
                        genes.append(gs)
            except Exception:
                pass
        # AMR'den birkac + virulanstan birkac (dengeli), MAX_TARGETS ile sinirli
        amr_part = genes[:amr_n][:MAX_TARGETS - 2] if amr_n else []
        vir_part = genes[amr_n:][:MAX_TARGETS - len(amr_part)]
        return (amr_part + vir_part)[:MAX_TARGETS]

    def run(self):
        self.check_inputs()
        std_dir = self.sub_dir("04_standardized")
        work = self.sub_dir("02_work")
        r = self.ctx.runner
        env_bio = util.ENV.get("comparative", "base")

        query_gbk = self.ctx.run_dir / "M06_GENOME_ANNOTATION" / "annotation.gbk"
        ref_json = self.ctx.run_dir / "M05_SPECIES_REFERENCE_IDENTIFICATION" / "closest_5_strains.json"
        refs = []
        if ref_json.exists():
            try:
                for x in (json.load(open(ref_json)) or [])[:N_REFS]:
                    gbk = Path(x.get("fasta_path", "")).with_suffix(".gbff")
                    if gbk.exists():
                        refs.append({"label": x["assembly_accession"], "gbk": str(gbk)})
            except Exception:
                refs = []

        targets = self._targets()
        if not query_gbk.exists() or not refs or not targets:
            reason = (f"Hedefli clinker icin: sorgu GBK({query_gbk.exists()}) + anotasyonlu referans({len(refs)}) "
                      f"+ hedef gen({len(targets)}) gerekir -> NOT_APPLICABLE.")
            self.write_summary(status="NOT_APPLICABLE", details={"reason": reason})
            return

        # 1) Hedef genlerin komsulugunu her genomdan kes (Biopython, ali-comparative)
        (work / "refs.json").write_text(json.dumps(refs), encoding="utf-8")
        (work / "targets.json").write_text(json.dumps(targets), encoding="utf-8")
        helper = work / "slice_loci.py"
        helper.write_text(SLICER, encoding="utf-8")
        loci_dir = work / "loci"
        r.run("slice_loci", ["python", str(helper), str(query_gbk), str(work / "refs.json"),
                             str(work / "targets.json"), str(WINDOW), str(loci_dir)],
              conda_env=env_bio, check=False)

        tsum = loci_dir / "targets_summary.json"
        if not tsum.exists():
            self.write_summary(status="WARNING", warnings=["Lokus cikarimi (Biopython) basarisiz."])
            return
        target_info = json.load(open(tsum))

        # 2) Lokus basina clinker (>=2 genom olan hedefler)
        produced = []
        vis = self.sub_dir("06_visualization")
        for ti in target_info:
            if not ti.get("found_in_query") or len(ti.get("genomes", [])) < 2:
                continue
            gene = ti["gene"]; tdir = Path(ti["dir"])
            html = vis / f"clinker_{gene.replace('/','_').replace(chr(39),'')}.html"
            # QUERY DAIMA USTTE: once QUERY.gbk, sonra referanslar (clinker girdi sirasini korur)
            q = tdir / "QUERY.gbk"
            refs_gbk = sorted(str(p) for p in tdir.glob("*.gbk") if p.name != "QUERY.gbk")
            gbks = ([str(q)] if q.exists() else []) + refs_gbk
            prov = r.run(f"clinker_{gene[:16]}", ["clinker", *gbks, "-p", str(html)],
                         conda_env=util.ENV.get("clinker", "base"), check=False)
            if html.exists() and html.stat().st_size > 1000:
                produced.append({"gene": gene, "genomes": ti["genomes"], "html": str(html)})

        with open(std_dir / "genomic_context.json", "w", encoding="utf-8") as f:
            json.dump({"targets_examined": targets, "loci_with_synteny": produced}, f, indent=2, ensure_ascii=False)
        with open(std_dir / "gene_neighborhoods.tsv", "w", encoding="utf-8") as f:
            f.write("Target_gene\tGenomes_compared\tClinker_HTML\n")
            for p in produced:
                f.write(f"{p['gene']}\t{len(p['genomes'])}\t{Path(p['html']).name}\n")

        if produced:
            self.write_summary(
                status="PASS",
                statistics={"loci_compared": len(produced),
                            "genes": [p["gene"] for p in produced]},
                details={"doi": "10.1093/bioinformatics/btab007",
                         "clinker_htmls": [p["html"] for p in produced]},
            )
        else:
            self.write_summary(status="WARNING",
                               warnings=["Hedef genler icin >=2 genomda ortak lokus bulunamadi; clinker uretilmedi."])
