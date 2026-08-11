"""M09 -- Virulence Analysis
Screening against VFDB (Virulence Factor Database).
Outputs: virulence_genes.tsv, virulence.json, M09_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


class VirulenceModule(Module):
    number = "09"
    name = "virulence"
    folder = "M09_VIRULENCE"
    enabled_key = "vfdb"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "virulence_genes.tsv"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        E = util.ENV

        vf_out = self.sub_dir("03_native_outputs") / "abricate_vfdb.tsv"
        prov = r.run("abricate_vfdb", ["abricate", "--db", "vfdb", str(genome)],
                     conda_env=E.get("virulence", "base"), version_cmd=["abricate", "--version"],
                     stdout_path=str(vf_out), check=False)
        tool_ran = prov.get("exit_code") == 0

        virulence_genes = []
        if vf_out.exists():
            with open(vf_out, "r", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 14:
                        virulence_genes.append({
                            "contig": parts[1],
                            "start": parts[2],
                            "end": parts[3],
                            "gene_symbol": parts[5],
                            "coverage": parts[9],
                            "identity": parts[10],
                            "category": parts[13] if len(parts) > 13 else "unknown"
                        })

        with open(std_dir / "virulence_genes.tsv", "w", encoding="utf-8") as fh:
            fh.write("Gene_Symbol\tCategory\tIdentity\tCoverage\tContig\tStart\tEnd\n")
            for v in virulence_genes:
                fh.write(f"{v['gene_symbol']}\t{v['category']}\t{v['identity']}\t{v['coverage']}\t{v['contig']}\t{v['start']}\t{v['end']}\n")

        with open(std_dir / "virulence.json", "w", encoding="utf-8") as fh:
            json.dump({"virulence_genes": virulence_genes}, fh, indent=2)

        if tool_ran:
            self.write_summary(
                status="PASS",
                statistics={"virulence_gene_count": len(virulence_genes)},
                details={"database": "VFDB (abricate)",
                         "categories": sorted(set(v["category"] for v in virulence_genes))},
            )
        else:
            self.write_summary(
                status="WARNING",
                statistics={"virulence_gene_count": len(virulence_genes)},
                warnings=[f"abricate/VFDB başarısız (exit {prov.get('exit_code')}). Log: {prov.get('log')}"],
            )
