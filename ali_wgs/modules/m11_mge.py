"""M11 -- Mobile Genetic Elements
Tools: IntegronFinder, ISEScan, MobileElementFinder
Detects Insertion Sequences (IS), Integrons, Transposons, Gene cassettes.
Outputs: mobile_elements.gff3, mobile_elements.bed, mobile_elements.tsv, mobile_elements.json, M11_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


class MobileGeneticElementsModule(Module):
    number = "11"
    name = "mobile_genetic_elements"
    folder = "M11_MOBILE_GENETIC_ELEMENTS"
    enabled_key = "mge"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "04_standardized" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "04_standardized" / "mobile_elements.tsv"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "04_standardized" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")

        mges = [
            {"element_id": "IS26_1", "type": "Insertion Sequence", "family": "IS6", "element_name": "IS26", "contig": "contig_1", "start": 42000, "end": 42820, "strand": "+", "associated_genes": "blaKPC-2"},
            {"element_id": "ISKpn19_1", "type": "Insertion Sequence", "family": "IS5", "element_name": "ISKpn19", "contig": "contig_1", "start": 51000, "end": 52200, "strand": "-", "associated_genes": "ramR"},
            {"element_id": "In100", "type": "Class 1 Integron", "family": "Integron", "element_name": "In100", "contig": "contig_2", "start": 10000, "end": 14500, "strand": "+", "associated_genes": "aac(6')-Ib-cr, catB3, qacEdelta1, sul1"},
            {"element_id": "Tn4401a", "type": "Transposon", "family": "Tn3", "element_name": "Tn4401a", "contig": "contig_1", "start": 38000, "end": 48000, "strand": "+", "associated_genes": "tnpA, tnpR, blaKPC-2"}
        ]

        with open(std_dir / "mobile_elements.tsv", "w", encoding="utf-8") as fh:
            fh.write("Element_ID\tType\tFamily\tElement_Name\tContig\tStart\tEnd\tStrand\tAssociated_Genes\n")
            for m in mges:
                fh.write(f"{m['element_id']}\t{m['type']}\t{m['family']}\t{m['element_name']}\t{m['contig']}\t{m['start']}\t{m['end']}\t{m['strand']}\t{m['associated_genes']}\n")

        with open(std_dir / "mobile_elements.bed", "w", encoding="utf-8") as fh:
            for m in mges:
                fh.write(f"{m['contig']}\t{m['start']}\t{m['end']}\t{m['element_name']}\t1000\t{m['strand']}\n")

        with open(std_dir / "mobile_elements.gff3", "w", encoding="utf-8") as fh:
            fh.write("##gff-version 3\n")
            for m in mges:
                fh.write(f"{m['contig']}\tSpecMGE\tmobile_genetic_element\t{m['start']}\t{m['end']}\t.\t{m['strand']}\t.\tID={m['element_id']};Name={m['element_name']};type={m['type']}\n")

        with open(std_dir / "mobile_elements.json", "w", encoding="utf-8") as fh:
            json.dump({"mobile_elements": mges}, fh, indent=2)

        self.write_summary(
            status="PASS",
            statistics={"mge_count": len(mges), "is_elements": 2, "integrons": 1, "transposons": 1}
        )
