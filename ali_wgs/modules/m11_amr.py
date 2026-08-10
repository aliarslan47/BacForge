"""Modül 11 — AMR (AMRFinderPlus). Bakteri ve faj contig'lerinde ortak (lizojenik ARG).
INPUT : 07_Contig_Filtering/contigs.filtered.fasta
OUTPUT: 11_AMR/amr.tsv
"""
from __future__ import annotations

from .base import Module
from .. import util


class AMRModule(Module):
    number = "11"
    name = "amr"
    folder = "11_AMR"
    enabled_key = "amr"

    def inputs(self):
        return [util.filtered_contigs(self.ctx)]

    def outputs(self):
        return [self.out_dir / "amr.tsv"]

    def run(self):
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        dbp = self.ctx.config["paths"]["db"]
        contigs = util.filtered_contigs(self.ctx)
        r.run("amrfinder", ["amrfinder", "-n", str(contigs), "--plus",
                           "--threads", str(t), "-o", str(self.outputs()[0]),
                           "-d", f"{dbp}/amrfinderplus/latest"],
              conda_env=E["amrfinder"], version_cmd=["amrfinder", "--version"],
              db_version="amrfinderplus-latest", check=False)
        if not self.outputs()[0].exists():
            self.outputs()[0].write_text("# AMRFinderPlus çalışmadı veya bulgu yok\n")
