"""M03 -- Genome Assembly
SHORT_READ: SPAdes
LONG_READ: Flye
HYBRID: Unicycler
ASSEMBLY_INPUT: Copy input FASTA directly
Outputs: 04_standardized/draft_genome.fasta, M03_summary.json
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from .base import Module
from .. import util


class AssemblyModule(Module):
    number = "03"
    name = "genome_assembly"
    folder = "M03_GENOME_ASSEMBLY"
    enabled_key = "assembly"

    def inputs(self):
        return [self.ctx.run_dir / "M00_INPUT_AUTO_DETECTION" / "data_type.json"]

    def outputs(self):
        return [self.out_dir / "draft_genome.fasta"]

    def run(self):
        self.check_inputs()
        data_type = self.ctx.detection.get("data_type", "SHORT_READ")
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        inp = Path(self.ctx.input_path)

        std_dir = self.sub_dir("04_standardized")
        draft_fasta = std_dir / "draft_genome.fasta"

        # 1. ASSEMBLY_INPUT
        if data_type == "ASSEMBLY_INPUT" or util.is_fasta_input(self.ctx):
            raw_fastas = util.raw_fasta_files(self.ctx)
            if raw_fastas:
                shutil.copy(raw_fastas[0], draft_fasta)
            else:
                shutil.copy(inp, draft_fasta)
            self.write_summary(status="PASS", details={"assembler": "None (Pre-assembled FASTA)"})
            return

        # 2. LONG_READ (Flye)
        if data_type == "LONG_READ":
            flye_dir = self.sub_dir("02_work") / "flye"
            long_fq = None
            if inp.is_dir():
                for f in inp.glob("*"):
                    if any(k in f.name.lower() for k in ["long", "ont", "nanopore", "fastq", "fq"]):
                        long_fq = f
                        break
            elif inp.is_file():
                long_fq = inp

            if not long_fq:
                # Check M01 filtered output
                long_fq = self.ctx.run_dir / "M01_READ_QC_PREPROCESSING" / "filtered_long.fastq.gz"

            r.run("flye", ["flye", "--nano-hq", str(long_fq), "-o", str(flye_dir), "-t", str(t)],
                  conda_env=E["flye"], version_cmd=["flye", "--version"])

            shutil.copy(flye_dir / "assembly.fasta", draft_fasta)
            self.write_summary(status="PASS", details={"assembler": "Flye"})
            return

        # 3. SHORT_READ (SPAdes / SKESA)
        if data_type == "SHORT_READ":
            clean_r1 = self.ctx.run_dir / "M01_READ_QC_PREPROCESSING" / "clean_R1.fastq.gz"
            clean_r2 = self.ctx.run_dir / "M01_READ_QC_PREPROCESSING" / "clean_R2.fastq.gz"

            if not clean_r1.exists():
                clean_r1 = inp

            spades_dir = self.sub_dir("02_work") / "spades"
            if clean_r2.exists():
                r.run("spades", ["spades.py", "--isolate", "-1", str(clean_r1), "-2", str(clean_r2), "-o", str(spades_dir), "-t", str(t)],
                      conda_env=E["asm_sr"], version_cmd=["spades.py", "--version"])
            else:
                r.run("spades", ["spades.py", "--isolate", "-s", str(clean_r1), "-o", str(spades_dir), "-t", str(t)],
                      conda_env=E["asm_sr"], version_cmd=["spades.py", "--version"])

            shutil.copy(spades_dir / "contigs.fasta", draft_fasta)
            self.write_summary(status="PASS", details={"assembler": "SPAdes"})
            return

        # 4. HYBRID (Unicycler / SPAdes hybrid)
        if data_type == "HYBRID":
            clean_r1 = self.ctx.run_dir / "M01_READ_QC_PREPROCESSING" / "clean_R1.fastq.gz"
            clean_r2 = self.ctx.run_dir / "M01_READ_QC_PREPROCESSING" / "clean_R2.fastq.gz"
            filtered_long = self.ctx.run_dir / "M01_READ_QC_PREPROCESSING" / "filtered_long.fastq.gz"

            uni_dir = self.sub_dir("02_work") / "unicycler"
            r.run("unicycler", ["unicycler", "-1", str(clean_r1), "-2", str(clean_r2), "-l", str(filtered_long), "-o", str(uni_dir), "-t", str(t)],
                  conda_env=E["asm_sr"], version_cmd=["unicycler", "--version"])

            shutil.copy(uni_dir / "assembly.fasta", draft_fasta)
            self.write_summary(status="PASS", details={"assembler": "Unicycler"})
            return
