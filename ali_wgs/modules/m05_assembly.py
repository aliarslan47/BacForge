"""Modül 05 — Assembly (platforma göre otomatik) + polishing.
INPUT : filtrelenmiş okuma (04) — yoksa ham
OUTPUT: 05_Assembly/assembly.fasta  (+ assembly_info.txt varsa)
Dallanma: ONT->Flye->Medaka · Illumina->SKESA · HiFi->hifiasm · FASTA->kopyala(atla)
"""
from __future__ import annotations

import shutil
from pathlib import Path

from .base import Module
from .. import util


class AssemblyModule(Module):
    number = "05"
    name = "assembly"
    folder = "05_Assembly"
    enabled_key = "assembly"

    def inputs(self):
        if util.is_fasta_input(self.ctx):
            return util.raw_fasta_files(self.ctx)
        return util.reads_for_assembly(self.ctx)

    def outputs(self):
        return [self.out_dir / "assembly.fasta"]

    def run(self):
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        final = self.outputs()[0]
        plat = util.platform(self.ctx)

        # --- FASTA girdi: assembly atla, kopyala ---
        if util.is_fasta_input(self.ctx):
            shutil.copy(util.raw_fasta_files(self.ctx)[0], final)
            (self.out_dir / "NOTE.txt").write_text("FASTA girdi: assembly atlandı, kopyalandı\n")
            self._split_contigs(final)
            return

        reads = util.reads_for_assembly(self.ctx)

        if plat in ("ONT", "PacBio_HiFi"):
            self._long_read(plat, reads[0], t, r, E, final)
        elif plat == "Illumina":
            self._illumina(reads, t, r, E, final)
        else:
            raise RuntimeError(f"[05_assembly] desteklenmeyen/tanımsız platform: {plat}")

        # Tüm contig'ler ayrı ayrı dosyalanır (ne olur olmaz saklansın)
        self._split_contigs(final)

    def _split_contigs(self, final):
        """assembly.fasta'daki her contig'i 05_Assembly/contigs/<contig>.fasta olarak yaz."""
        if not final.exists():
            return
        seqs = util.read_fasta(final)
        d = self.out_dir / "contigs"
        d.mkdir(exist_ok=True)
        for name, seq in seqs.items():
            safe = name.replace("/", "_").replace("|", "_")
            util.write_fasta({name: seq}, d / f"{safe}.fasta")
        (d / "INDEX.txt").write_text(
            "Her contig ayrı FASTA olarak burada.\n"
            + "\n".join(f"{n}\t{len(s)} bp" for n, s in seqs.items()) + "\n")

    # --- ONT/HiFi: Flye -> Medaka ---
    def _long_read(self, plat, reads, t, r, E, final):
        flye_dir = self.out_dir / "flye"
        mode = "--pacbio-hifi" if plat == "PacBio_HiFi" else "--nano-hq"
        r.run("flye", ["flye", mode, str(reads), "-o", str(flye_dir),
                       "-t", str(t)],
              conda_env=E["flye"], version_cmd=["flye", "--version"])
        draft = flye_dir / "assembly.fasta"
        info = flye_dir / "assembly_info.txt"
        if info.exists():
            shutil.copy(info, self.out_dir / "assembly_info.txt")  # circular bilgisi

        if plat == "PacBio_HiFi":
            shutil.copy(draft, final)  # HiFi: polishing genelde gereksiz
            return

        # ONT: Medaka polishing
        medaka_dir = self.out_dir / "medaka"
        prov = r.run("medaka", ["medaka_consensus", "-i", str(reads), "-d", str(draft),
                               "-o", str(medaka_dir), "-t", str(min(t, 8))],
                     conda_env=E["medaka"], version_cmd=["medaka", "--version"], check=False)
        polished = medaka_dir / "consensus.fasta"
        if prov["exit_code"] == 0 and polished.exists():
            shutil.copy(polished, final)
        else:  # medaka başarısızsa Flye taslağıyla devam et (pipeline durmaz)
            shutil.copy(draft, final)
            (self.out_dir / "WARN.txt").write_text("Medaka başarısız; Flye taslağı kullanıldı\n")

    # --- Illumina: SKESA ---
    def _illumina(self, reads, t, r, E, final):
        reads_arg = ",".join(str(x) for x in reads)
        r.run("skesa", ["skesa", "--reads", reads_arg, "--cores", str(t),
                       "--contigs_out", str(final)],
              conda_env=E["asm_sr"], version_cmd=["skesa", "--version"], check=False)
        if not final.exists() or final.stat().st_size == 0:
            # SKESA başarısız -> SPAdes isolate fallback
            spades_dir = self.out_dir / "spades"
            r.run("spades", ["spades.py", "--isolate", "-1", str(reads[0]),
                            "-2", str(reads[1]), "-o", str(spades_dir),
                            "-t", str(t)],
                  conda_env=E["asm_sr"], version_cmd=["spades.py", "--version"])
            shutil.copy(spades_dir / "contigs.fasta", final)
