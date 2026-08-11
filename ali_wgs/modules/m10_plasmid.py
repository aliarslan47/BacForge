"""M10 -- Plasmid Analysis
Tools: MOB-suite (mob_recon) & PlasmidFinder
Outputs: plasmids.fasta, plasmids.tsv, replicons.tsv, plasmids.json, M10_summary.json
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from .base import Module
from .. import util


class PlasmidModule(Module):
    number = "10"
    name = "plasmid_analysis"
    folder = "M10_PLASMID"
    enabled_key = "plasmid"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "plasmids.tsv"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        E = util.ENV
        t = util.threads(self.ctx)

        mob_dir = self.sub_dir("02_work") / "mob_recon"
        prov = r.run("mob_recon", ["mob_recon", "--infile", str(genome), "--outdir", str(mob_dir), "--num_threads", str(t)],
                     conda_env=E["mobsuite"], version_cmd=["mob_recon", "--version"], check=False)
        tool_ran = prov.get("exit_code") == 0

        mob_report = mob_dir / "contig_report.txt"
        plasmids = []
        if mob_report.exists():
            with open(mob_report, "r", encoding="utf-8") as fh:
                header = fh.readline().strip().split("\t")
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 6 and parts[3] == "plasmid":
                        plasmids.append({
                            "contig": parts[0],
                            "size": parts[1],
                            "plasmid_id": parts[5] if len(parts) > 5 else "unknown",
                            "rep_type": parts[6] if len(parts) > 6 else ""
                        })

        with open(std_dir / "plasmids.tsv", "w", encoding="utf-8") as f:
            f.write("Plasmid_ID\tContig\tSize\tRep_Type\n")
            for p in plasmids:
                f.write(f"{p['plasmid_id']}\t{p['contig']}\t{p['size']}\t{p['rep_type']}\n")
        with open(std_dir / "plasmids.json", "w", encoding="utf-8") as f:
            json.dump({"plasmids": plasmids}, f, indent=2, ensure_ascii=False)

        # Dürüst durum: araç çalıştıysa (plazmid 0 olsa da) PASS; çalışmadıysa WARNING.
        if tool_ran:
            self.write_summary(status="PASS", statistics={"plasmid_count": len(plasmids)},
                               details={"tool": "MOB-suite (mob_recon)", "note": "PlasmidFinder: Milestone 2"})
        else:
            self.write_summary(status="WARNING", statistics={"plasmid_count": len(plasmids)},
                               warnings=[f"mob_recon başarısız (exit {prov.get('exit_code')}). Log: {prov.get('log')}"])
        return
