"""M01 -- Read QC & Preprocessing
Short-read QC (fastp / FastQC / MultiQC) -> clean_R1.fastq.gz, clean_R2.fastq.gz
Long-read QC (NanoPlot / Filtlong) -> filtered_long.fastq.gz
Hybrid QC (both short & long branches)
Assembly input -> SKIPPED
"""
from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

from .base import Module
from .. import util


class ReadQCModule(Module):
    number = "01"
    name = "read_qc_preprocessing"
    folder = "M01_READ_QC_PREPROCESSING"
    enabled_key = "qc"

    def inputs(self):
        return [self.ctx.run_dir / "M00_INPUT_AUTO_DETECTION" / "data_type.json"]

    def outputs(self):
        return [self.out_dir / "M01_summary.json"]

    def run(self):
        self.check_inputs()
        data_type = self.ctx.detection.get("data_type", "SHORT_READ")
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        inp = Path(self.ctx.input_path)

        std_dir = self.sub_dir("04_standardized")
        stats = {}

        if data_type == "ASSEMBLY_INPUT":
            self.write_summary(status="SKIPPED", details={"reason": "Input is FASTA assembly; read QC skipped."})
            return

        # 1. SHORT READ / HYBRID Branch
        if data_type in ("SHORT_READ", "HYBRID"):
            r1, r2 = None, None
            if inp.is_dir():
                files = sorted(list(inp.glob("*.fastq*")) + list(inp.glob("*.fq*")))
                for f in files:
                    if "_R1" in f.name or "_1." in f.name:
                        r1 = f
                    elif "_R2" in f.name or "_2." in f.name:
                        r2 = f
            elif inp.is_file():
                r1 = inp

            clean_r1 = self.sub_dir("04_standardized") / "clean_R1.fastq.gz"
            clean_r2 = self.sub_dir("04_standardized") / "clean_R2.fastq.gz"

            if r1 and r2:
                cmd = [
                    "fastp", "-i", str(r1), "-I", str(r2),
                    "-o", str(clean_r1), "-O", str(clean_r2),
                    "-h", str(self.sub_dir("06_visualization") / "fastp_report.html"),
                    "-j", str(self.sub_dir("04_standardized") / "fastp.json"),
                    "-w", str(t)
                ]
                r.run("fastp", cmd, conda_env=E["illumina_qc"], version_cmd=["fastp", "--version"])
                stats["short_read_qc"] = "COMPLETED"
            elif r1:
                cmd = [
                    "fastp", "-i", str(r1), "-o", str(clean_r1),
                    "-h", str(self.sub_dir("06_visualization") / "fastp_report.html"),
                    "-j", str(self.sub_dir("04_standardized") / "fastp.json"),
                    "-w", str(t)
                ]
                r.run("fastp", cmd, conda_env=E["illumina_qc"], version_cmd=["fastp", "--version"])
                stats["short_read_qc"] = "COMPLETED"

        # 2. LONG READ / HYBRID Branch
        if data_type in ("LONG_READ", "HYBRID"):
            long_fq = None
            if inp.is_dir():
                for f in inp.glob("*"):
                    if any(k in f.name.lower() for k in ["long", "ont", "nanopore", "fastq", "fq"]):
                        long_fq = f
                        break
            elif inp.is_file():
                long_fq = inp

            filtered_long = self.sub_dir("04_standardized") / "filtered_long.fastq.gz"

            if long_fq and long_fq.exists():
                # NanoPlot: ham uzun-okuma QC görseli/istatistiği (uzunluk/kalite dağılımı)
                viz = self.sub_dir("06_visualization")
                r.run("nanoplot", ["NanoPlot", "--fastq", str(long_fq), "-o", str(viz),
                                   "-t", str(t), "-p", "nanoplot_"],
                      conda_env=E["ont_qc"], version_cmd=["NanoPlot", "--version"], check=False)

                cfg_filt = self.ctx.config.get("tools", {}).get("filtlong", {})
                min_len = cfg_filt.get("min_length", 1000)
                keep_pct = cfg_filt.get("keep_percent", 95)
                # filtlong DÜZ-METİN FASTQ verir (stdout); önce .fastq'a yaz, sonra gerçekten gzip'le
                # (aksi halde .gz adlı sıkıştırılmamış dosya Flye/Unicycler'da "Not a gzipped file" hatası verir).
                filtered_plain = self.sub_dir("04_standardized") / "filtered_long.fastq"
                cmd = ["filtlong", "--min_length", str(min_len), "--keep_percent", str(keep_pct), str(long_fq)]
                r.run("filtlong", cmd, conda_env=E["ont_qc"], version_cmd=["filtlong", "--version"], stdout_path=str(filtered_plain))
                with open(filtered_plain, "rb") as _fi, gzip.open(filtered_long, "wb") as _fo:
                    shutil.copyfileobj(_fi, _fo)
                filtered_plain.unlink(missing_ok=True)
                stats["long_read_qc"] = "COMPLETED"

        # Write summary TSV
        with open(std_dir / "qc_statistics.tsv", "w", encoding="utf-8") as fh:
            fh.write("Metric\tValue\n")
            for k, v in stats.items():
                fh.write(f"{k}\t{v}\n")

        self.write_summary(status="PASS", statistics=stats)
