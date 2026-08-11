"""M14 -- Genomic Context & NCBI Closest-5 Comparison (clinker)
Clinker ile hedef gen komşuluğunun query vs referans genomlar arası synteny karşılaştırması.
DOI: 10.1093/bioinformatics/btab007
KOŞUL-BAZLI DÜRÜST: clinker en az 2 ANOTASYONLU genom (GFF/GBK) ister. Tek örnekli turda
referanslar anotasyonsuz -> NOT_APPLICABLE. Sahte PASS/'Bulunamadı' YAZILMAZ.
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


class GenomicContextModule(Module):
    number = "14"
    name = "genomic_context"
    folder = "M14_GENOMIC_CONTEXT"
    enabled_key = "clinker"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "M14_summary.json"]

    def run(self):
        self.check_inputs()
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        E = util.ENV
        t = util.threads(self.ctx)

        # Referans anotasyonlu genom sayısı (clinker için gerekli)
        ref_json = self.ctx.run_dir / "M05_SPECIES_REFERENCE_IDENTIFICATION" / "closest_5_strains.json"
        refs = []
        if ref_json.exists():
            try:
                refs = json.load(open(ref_json)) or []
            except Exception:
                refs = []

        # clinker query + en az 1 referans için GBK/GFF ister; referanslar bu turda anotasyonsuz.
        annotated_refs = [x for x in refs if x.get("gbk_path") or x.get("gff_path")]
        if len(annotated_refs) < 1:
            reason = ("clinker en az 2 anotasyonlu genom (query + >=1 referans, GBK/GFF) ister. "
                      "Tek-örnekli turda referanslar anotasyonsuz -> NOT_APPLICABLE. "
                      "Milestone 2 (batch/NCBI referans anotasyonu) ile gerçek çalışacak.")
            with open(std_dir / "gene_neighborhoods.tsv", "w", encoding="utf-8") as f:
                f.write("# NOT_APPLICABLE\t" + reason + "\n")
                f.write("Gene_ID\tContig\tStart\tEnd\tNeighborhood\n")
            self.write_summary(status="NOT_APPLICABLE", details={"reason": reason,
                               "annotated_reference_count": len(annotated_refs)})
            return

        # (Milestone 2) Gerçek clinker akışı burada koşacak — query+ref GBK'ler ile:
        clinker_out = self.sub_dir("02_work") / "clinker"
        clinker_out.mkdir(parents=True, exist_ok=True)
        gbks = [str(self.ctx.run_dir / "M06_GENOME_ANNOTATION" / "annotation.gbk")] + \
               [x["gbk_path"] for x in annotated_refs if x.get("gbk_path")]
        aln_html = self.sub_dir("06_visualization") / "clinker_alignment.html"
        prov = r.run("clinker", ["clinker", *gbks, "-p", str(aln_html)],
                     conda_env=E.get("clinker", "base"), version_cmd=["clinker", "--version"], check=False)
        status = "PASS" if prov.get("exit_code") == 0 and aln_html.exists() else "WARNING"
        with open(std_dir / "gene_neighborhoods.tsv", "w", encoding="utf-8") as f:
            f.write("Gene_ID\tContig\tStart\tEnd\tNeighborhood\n")
        self.write_summary(status=status, statistics={"compared_genomes": len(gbks)},
                           details={"clinker_html": str(aln_html)})
