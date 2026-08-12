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
from ..detect import detect_ont_chemistry


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

        # 2. LONG_READ (Flye --nano-hq + Medaka polishing)
        if data_type == "LONG_READ":
            flye_dir = self.sub_dir("02_work") / "flye"

            # Tercih: M01'in QC-filtrelenmiş uzun okuması; yoksa ham ONT
            filtered_long = self.ctx.run_dir / "M01_READ_QC_PREPROCESSING" / "filtered_long.fastq.gz"
            if filtered_long.exists() and filtered_long.stat().st_size > 0:
                long_fq = filtered_long
            else:
                long_fq = util.find_long_reads(inp)

            if not long_fq or not Path(long_fq).exists():
                raise FileNotFoundError(
                    "LONG_READ: uzun okuma bulunamadı (ne M01 filtered_long ne de ham ONT)."
                )

            # ONT kimyasini ham/filtrelenmis okumadan tespit et -> Flye modu + Medaka modeli
            chem = detect_ont_chemistry(long_fq, self.ctx.config)

            r.run("flye", ["flye", chem["flye_mode"], str(long_fq), "-o", str(flye_dir), "-t", str(t)],
                  conda_env=E["flye"], version_cmd=["flye", "--version"])
            flye_asm = flye_dir / "assembly.fasta"

            # Medaka polishing -- ONT indel duzeltme; kimyaya gore model deneme kaskadi
            # (R9 -> r941 modeli; R10 -> once --bacteria, sonra acik SUP modeli).
            medaka_dir = self.sub_dir("02_work") / "medaka"
            warnings = []
            polisher = None
            polished_ok = False
            for attempt in chem["medaka_attempts"]:
                if medaka_dir.exists():
                    shutil.rmtree(medaka_dir, ignore_errors=True)
                mprov = r.run(
                    "medaka",
                    ["medaka_consensus", "-i", str(long_fq), "-d", str(flye_asm),
                     "-o", str(medaka_dir), "-t", str(t)] + attempt,
                    conda_env=E["medaka"], version_cmd=["medaka", "--version"], check=False,
                )
                cons = medaka_dir / "consensus.fasta"
                if mprov.get("exit_code") == 0 and cons.exists() and cons.stat().st_size > 0:
                    shutil.copy(cons, draft_fasta)
                    polisher = f"Medaka ({' '.join(attempt)})"
                    polished_ok = True
                    break
                warnings.append(
                    f"Medaka denemesi basarisiz ({' '.join(attempt)}, exit {mprov.get('exit_code')})."
                )

            if not polished_ok:
                # Dürüstlük: polishing başarısızsa gizleme, cilalanmamış Flye ile devam et + WARNING
                shutil.copy(flye_asm, draft_fasta)
                warnings.append("Tum Medaka denemeleri basarisiz; cilalanmamis Flye assembly kullanildi.")

            self.write_summary(
                status="PASS" if polished_ok else "WARNING",
                details={"assembler": f"Flye ({chem['flye_mode']})", "polisher": polisher,
                         "polishing_performed": polished_ok,
                         "ont_chemistry": chem["chemistry"],
                         "chemistry_basis": chem["basis"],
                         "chemistry_confidence": chem["confidence"]},
                warnings=warnings,
            )
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
