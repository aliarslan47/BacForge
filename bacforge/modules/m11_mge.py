"""M11 -- Mobile Genetic Elements (ISEScan; IntegronFinder/MEFinder Milestone 2)
KURAL: 'Bulunamadı' satırı + otomatik PASS YOK. Araç çalışıp bulgu bulamazsa -> PASS (count=0).
Araç çalışmazsa -> WARNING. Boş tabloda yalnız başlık yazılır.
"""
from __future__ import annotations

import json
import shutil
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

        for el in elements:
            el.setdefault("source", "ISEScan")
        is_count = len(elements)

        # IntegronFinder (ali-mge) — integron tespiti (.summary: ID_replicon CALIN complete In0 topology size)
        intf_dir = self.sub_dir("02_work") / "integron_finder"
        if intf_dir.exists():
            shutil.rmtree(intf_dir, ignore_errors=True)  # --outdir mevcutsa integron_finder hata verir
        prov_if = r.run("integron_finder", ["integron_finder", str(genome), "--outdir", str(intf_dir),
                                            "--cpu", str(t), "--mute"],
                        conda_env=E.get("mge", "base"), version_cmd=["integron_finder", "--version"], check=False)
        integron_count = 0
        for summ in intf_dir.rglob("*.summary"):
            for ln in summ.read_text(encoding="utf-8", errors="replace").splitlines():
                if ln.startswith("#") or ln.startswith("ID_replicon") or not ln.strip():
                    continue
                c = ln.split("\t")
                if len(c) >= 4:
                    try:
                        calin, complete, in0 = int(c[1]), int(c[2]), int(c[3])
                    except ValueError:
                        continue
                    tot = calin + complete + in0
                    if tot > 0:
                        integron_count += tot
                        elements.append({"element_id": f"integron_{c[0]}",
                                         "type": f"Integron (complete:{complete},CALIN:{calin},In0:{in0})",
                                         "contig": c[0], "start": "", "end": "", "source": "IntegronFinder"})
        if_ran = prov_if.get("exit_code") == 0

        with open(std_dir / "mobile_elements.tsv", "w", encoding="utf-8") as f:
            f.write("Element_ID\tType\tContig\tStart\tEnd\tSource\n")
            for el in elements:
                f.write(f"{el['element_id']}\t{el['type']}\t{el['contig']}\t{el['start']}\t{el['end']}\t{el.get('source','')}\n")
        with open(std_dir / "mobile_elements.json", "w", encoding="utf-8") as f:
            json.dump({"mobile_elements": elements,
                       "source_counts": {"ISEScan": is_count, "IntegronFinder": integron_count}},
                      f, indent=2, ensure_ascii=False)

        if tool_ran or if_ran:
            ran = [t for t, ok in (("ISEScan", tool_ran), ("IntegronFinder", if_ran)) if ok]
            self.write_summary(status="PASS",
                               statistics={"is_count": is_count, "integron_count": integron_count,
                                           "mge_count": len(elements)},
                               details={"tools": ran})
        else:
            self.write_summary(status="WARNING", statistics={"mge_count": len(elements)},
                               warnings=[f"ISEScan exit {prov.get('exit_code')}, IntegronFinder exit {prov_if.get('exit_code')}"])
