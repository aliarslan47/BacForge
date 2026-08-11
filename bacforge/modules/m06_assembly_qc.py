"""Modül 06 — Assembly QC.
INPUT : 05_Assembly/assembly.fasta (+ okuma, coverage için)
OUTPUT: 06_Assembly_QC/quast/report.tsv, coverage.tsv, circular.tsv
"""
from __future__ import annotations

from pathlib import Path

from .base import Module
from .. import util


class AssemblyQCModule(Module):
    number = "06"
    name = "assembly_qc"
    folder = "06_Assembly_QC"
    enabled_key = "assembly_qc"

    def inputs(self):
        return [util.assembly_fasta(self.ctx)]

    def outputs(self):
        return [self.out_dir / "quast" / "report.tsv"]

    def run(self):
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        asm = util.assembly_fasta(self.ctx)

        # QUAST
        r.run("quast", ["quast", str(asm), "-o", str(self.out_dir / "quast"),
                       "-t", str(t), "--silent"],
              conda_env=E["quast"], version_cmd=["quast", "--version"], check=False)

        # Coverage (okuma varsa)
        if not util.is_fasta_input(self.ctx):
            preset = {"ONT": "map-ont", "PacBio_HiFi": "map-hifi"}.get(
                util.platform(self.ctx), "sr")
            reads = util.reads_for_assembly(self.ctx)
            cov = self.out_dir / "coverage.tsv"
            cmd = (f"minimap2 -ax {preset} -t {t} '{asm}' "
                   + " ".join(f"'{x}'" for x in reads)
                   + f" | samtools sort -@ {t} -o '{self.out_dir}/aln.bam' - "
                   + f"&& samtools depth -a '{self.out_dir}/aln.bam' "
                   + "| awk '{s+=$3; n++} END{if(n>0) printf \"mean_coverage\\t%.1f\\n\", s/n}' "
                   + f"> '{cov}'")
            r.run("coverage", ["bash", "-c", cmd], conda_env=E["core"],
                  version_cmd=["minimap2", "--version"], check=False)

        # Circular (Flye assembly_info.txt'ten)
        info = Path(self.ctx.run_dir) / "05_Assembly" / "assembly_info.txt"
        if info.exists():
            circ = self.out_dir / "circular.tsv"
            with open(info) as fh, open(circ, "w") as out:
                out.write("contig\tlength\tcoverage\tcircular\n")
                for i, line in enumerate(fh):
                    if i == 0:
                        continue
                    c = line.split("\t")
                    if len(c) >= 4:
                        out.write(f"{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\n")
