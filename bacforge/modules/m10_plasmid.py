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

        for p in plasmids:
            p.setdefault("source", "MOB-suite")

        # PlasmidFinder (abricate --db plasmidfinder) — replikon tespiti (M04 genome üstünde)
        pf_out = self.sub_dir("02_work") / "abricate_plasmidfinder.tsv"
        prov_pf = r.run("abricate_plasmidfinder", ["abricate", "--db", "plasmidfinder", "--nopath", str(genome)],
                        conda_env=E.get("virulence", "base"), version_cmd=["abricate", "--version"],
                        stdout_path=str(pf_out), check=False)
        replicons = []
        if pf_out.exists() and pf_out.stat().st_size > 0:
            lines = pf_out.read_text(encoding="utf-8", errors="replace").splitlines()
            if lines:
                hdr = lines[0].lstrip("#").split("\t"); ix = {h.strip(): i for i, h in enumerate(hdr)}
                gv = lambda p, c: (p[ix[c]] if c in ix and ix[c] < len(p) else "")
                for ln in lines[1:]:
                    if not ln.strip():
                        continue
                    p = ln.split("\t")
                    replicons.append({"plasmid_id": gv(p, "GENE"), "contig": gv(p, "SEQUENCE"),
                                      "size": "", "rep_type": gv(p, "GENE"),
                                      "identity": gv(p, "%IDENTITY"), "source": "PlasmidFinder"})
        pf_ran = prov_pf.get("exit_code") == 0
        plasmids_all = plasmids + replicons

        with open(std_dir / "plasmids.tsv", "w", encoding="utf-8") as f:
            f.write("Plasmid_ID\tContig\tSize\tRep_Type\tSource\n")
            for p in plasmids_all:
                f.write(f"{p['plasmid_id']}\t{p['contig']}\t{p['size']}\t{p['rep_type']}\t{p.get('source','')}\n")
        with open(std_dir / "plasmids.json", "w", encoding="utf-8") as f:
            json.dump({"plasmids": plasmids_all, "replicons": replicons,
                       "source_counts": {"MOB-suite": len(plasmids), "PlasmidFinder": len(replicons)}},
                      f, indent=2, ensure_ascii=False)

        # Dürüst durum: en az bir araç çalıştıysa (0 bulsa da) PASS; hiçbiri değilse WARNING.
        if tool_ran or pf_ran:
            ran = [t for t, ok in (("MOB-suite", tool_ran), ("PlasmidFinder", pf_ran)) if ok]
            self.write_summary(status="PASS",
                               statistics={"plasmid_count": len(plasmids), "replicon_count": len(replicons)},
                               details={"tools": ran})
        else:
            self.write_summary(status="WARNING", statistics={"plasmid_count": 0},
                               warnings=[f"mob_recon exit {prov.get('exit_code')}, PlasmidFinder exit {prov_pf.get('exit_code')}"])
        return
