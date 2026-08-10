"""Modül 04 — Read Filtering (platforma göre otomatik).
INPUT : ham okuma
OUTPUT: 04_Filtering/filtered.fastq.gz  (Illumina: filtered_R1/_R2.fastq.gz)
Not: genom boyutu assembly öncesi bilinmediğinden coverage-subsample ERTELENDİ;
     Flye yüksek coverage'ı kendi içinde yönetir. Burada kalite/uzunluk filtresi yapılır.
"""
from __future__ import annotations

from .base import Module
from .. import util


class FilteringModule(Module):
    number = "04"
    name = "filtering"
    folder = "04_Filtering"
    enabled_key = "filtering"

    def _is_skip(self):
        return util.is_fasta_input(self.ctx)

    def inputs(self):
        return [] if self._is_skip() else util.raw_read_files(self.ctx)

    def outputs(self):
        if self._is_skip():
            return [self.out_dir / "SKIPPED"]
        if util.is_paired(self.ctx):
            return [self.out_dir / "filtered_R1.fastq.gz", self.out_dir / "filtered_R2.fastq.gz"]
        return [self.out_dir / "filtered.fastq.gz"]

    def run(self):
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        tools = self.ctx.config.get("tools", {})

        if self._is_skip():
            (self.out_dir / "SKIPPED").write_text("FASTA girdi: filtreleme atlandı\n")
            return

        reads = util.raw_read_files(self.ctx)
        plat = util.platform(self.ctx)

        if plat in ("ONT", "PacBio_HiFi"):
            fl = tools.get("filtlong", {})
            ml = fl.get("min_length", 1000)
            kp = fl.get("keep_percent", 95)
            out = self.outputs()[0]
            # filtlong stdout -> gzip; pipe için bash -c (env içinde çalışır)
            cmd = (f"filtlong --min_length {ml} --keep_percent {kp} "
                   f"'{reads[0]}' | gzip > '{out}'")
            r.run("filtlong", ["bash", "-c", cmd], conda_env=E["ont_qc"],
                  version_cmd=["filtlong", "--version"])
        elif plat == "Illumina" and util.is_paired(self.ctx):
            o1, o2 = self.outputs()
            r.run("fastp", ["fastp", "-i", str(reads[0]), "-I", str(reads[1]),
                           "-o", str(o1), "-O", str(o2), "--thread", str(min(t, 16)),
                           "-j", str(self.out_dir / "fastp.json"),
                           "-h", str(self.out_dir / "fastp.html")],
                  conda_env=E["ill_qc"], version_cmd=["fastp", "--version"])
        else:
            # tek-uçlu Illumina veya tanımsız: fastp tek dosya
            o1 = self.outputs()[0]
            r.run("fastp", ["fastp", "-i", str(reads[0]), "-o", str(o1),
                           "--thread", str(min(t, 16))],
                  conda_env=E["ill_qc"], version_cmd=["fastp", "--version"])
