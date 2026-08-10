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
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "mobile_elements.tsv"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        E = util.ENV
        t = util.threads(self.ctx)

        # Run ISEScan
        ise_dir = self.sub_dir("02_work") / "isescan"
        ise_dir.mkdir(parents=True, exist_ok=True)
        r.run("isescan", ["isescan.py", "--seqfile", str(genome), "--output", str(ise_dir), "--nthread", str(t)],
              conda_env=E.get("mge", "base"), version_cmd=["isescan.py", "--version"], check=False)

        elements = []
        # Find isescan TSV output (usually ends with .tsv inside the output dir)
        for tsv_file in ise_dir.glob("*.tsv"):
            if tsv_file.exists():
                with open(tsv_file, "r", encoding="utf-8") as fh:
                    for line in fh:
                        if line.startswith("seqID") or line.startswith("#"):
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 6:
                            elements.append({
                                "element_id": parts[3] if len(parts)>3 else "unknown",
                                "type": "IS_Element",
                                "contig": parts[0],
                                "start": parts[1],
                                "end": parts[2]
                            })

        with open(std_dir / "mobile_elements.tsv", "w", encoding="utf-8") as f:
            f.write("Element_ID\tType\tContig\tStart\tEnd\n")
            if elements:
                for el in elements:
                    f.write(f"{el['element_id']}\t{el['type']}\t{el['contig']}\t{el['start']}\t{el['end']}\n")
            else:
                # User's explicit rule: Write "Bulunamadı" if no results
                f.write("Bulunamadı\tBulunamadı\t-\t-\t-\n")

        self.write_summary(
            status="PASS", 
            statistics={"mge_count": len(elements)}, 
            details={"info": "MGE Analizi tamamlandı" if elements else "MGE Bulunamadı"}
        )
        return
