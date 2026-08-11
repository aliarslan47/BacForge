"""Modül 03 — QC (platforma göre otomatik).
INPUT : ham okuma (ctx.detection.files) veya FASTA
OUTPUT: 03_QC/qc_stats.tsv (+ NanoPlot/FastQC raporları)
"""
from __future__ import annotations

from pathlib import Path

from .base import Module
from .. import util


class QCModule(Module):
    number = "03"
    name = "qc"
    folder = "03_QC"
    enabled_key = "qc"

    def inputs(self):
        det = util.load_detection(self.ctx)
        return [Path(self.ctx.run_dir) / "01_Input" / "platform.json"] if det else [self.ctx.input_path]

    def outputs(self):
        return [self.out_dir / "qc_stats.tsv"]

    def run(self):
        t = util.threads(self.ctx)
        r = self.ctx.runner
        E = util.ENV

        if util.is_fasta_input(self.ctx):
            fa = util.raw_fasta_files(self.ctx)[0]
            r.run("seqkit_stats", ["seqkit", "stats", "-a", "-T", str(fa)],
                  conda_env=E["ont_qc"], version_cmd=["seqkit", "version"],
                  stdout_path=str(self.outputs()[0]))
            return

        reads = util.raw_read_files(self.ctx)
        plat = util.platform(self.ctx)

        # seqkit istatistik (her platform)
        r.run("seqkit_stats", ["seqkit", "stats", "-a", "-T", *map(str, reads)],
              conda_env=E["ont_qc"], version_cmd=["seqkit", "version"],
              stdout_path=str(self.outputs()[0]))

        if plat in ("ONT", "PacBio_HiFi"):
            r.run("nanoplot", ["NanoPlot", "--fastq", str(reads[0]),
                               "-o", str(self.out_dir / "nanoplot"),
                               "-t", str(t), "--tsv_stats", "--no_static"],
                  conda_env=E["ont_qc"], version_cmd=["NanoPlot", "--version"], check=False)
        elif plat == "Illumina":
            r.run("fastqc", ["fastqc", "-t", str(t), "-o", str(self.out_dir),
                            *map(str, reads)],
                  conda_env=E["ill_qc"], version_cmd=["fastqc", "--version"], check=False)
