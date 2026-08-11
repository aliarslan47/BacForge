"""M11 -- Mobile Genetic Elements (ISEScan; IntegronFinder/MEFinder Milestone 2)
KURAL: 'Bulunamadı' satırı + otomatik PASS YOK. Araç çalışıp bulgu bulamazsa -> PASS (count=0).
Araç çalışmazsa -> WARNING. Boş tabloda yalnız başlık yazılır.
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

        ise_dir = self.sub_dir("02_work") / "isescan"
        ise_dir.mkdir(parents=True, exist_ok=True)
        prov = r.run("isescan", ["isescan.py", "--seqfile", str(genome), "--output", str(ise_dir), "--nthread", str(t)],
                     conda_env=E.get("mge", "base"), version_cmd=["isescan.py", "--version"], check=False)
        tool_ran = prov.get("exit_code") == 0

        elements = []
        for tsv_file in ise_dir.rglob("*.tsv"):
            with open(tsv_file, "r", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("seqID") or line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 6:
                        elements.append({
                            "element_id": parts[3] if len(parts) > 3 else "unknown",
                            "type": "IS_Element",
                            "contig": parts[0],
                            "start": parts[1],
                            "end": parts[2],
                        })

        with open(std_dir / "mobile_elements.tsv", "w", encoding="utf-8") as f:
            f.write("Element_ID\tType\tContig\tStart\tEnd\n")
            for el in elements:
                f.write(f"{el['element_id']}\t{el['type']}\t{el['contig']}\t{el['start']}\t{el['end']}\n")
        with open(std_dir / "mobile_elements.json", "w", encoding="utf-8") as f:
            json.dump({"mobile_elements": elements}, f, indent=2, ensure_ascii=False)

        if tool_ran:
            self.write_summary(status="PASS", statistics={"mge_count": len(elements)},
                               details={"tool": "ISEScan", "note": "IntegronFinder/MEFinder: Milestone 2"})
        else:
            self.write_summary(status="WARNING", statistics={"mge_count": len(elements)},
                               warnings=[f"ISEScan başarısız (exit {prov.get('exit_code')}). Log: {prov.get('log')}"])
