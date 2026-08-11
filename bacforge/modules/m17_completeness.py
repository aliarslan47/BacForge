"""Modül 17 — Completeness (içerik-farkında routing).
INPUT : 08_Taxonomy/{bacterial,viral}.fasta
OUTPUT: 17_Completeness/quality.tsv (bakteri: CheckM2 · faj: CheckV)
"""
from __future__ import annotations

import glob
from pathlib import Path

from .base import Module
from .. import util


class CompletenessModule(Module):
    number = "17"
    name = "completeness"
    folder = "17_Completeness"
    enabled_key = "completeness"

    def _cls(self, name):
        return Path(self.ctx.run_dir) / "08_Taxonomy" / name

    def inputs(self):
        return [self._cls("classification.tsv")]

    def outputs(self):
        return [self.out_dir / "quality.tsv"]

    def run(self):
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        dbp = self.ctx.config["paths"]["db"]
        bacterial = self._cls("bacterial.fasta")
        viral = self._cls("viral.fasta")
        summary = []

        # Bakteri/arkea -> CheckM2
        if util.count_fasta_seqs(bacterial) > 0:
            dmnd = glob.glob(f"{dbp}/checkm2/**/*.dmnd", recursive=True)
            out = self.out_dir / "checkm2"
            cmd = ["checkm2", "predict", "--input", str(bacterial),
                   "--output-directory", str(out), "-t", str(t), "--force"]
            if dmnd:
                cmd += ["--database_path", dmnd[0]]
            r.run("checkm2", cmd, conda_env=E["checkm2"],
                  version_cmd=["checkm2", "--version"], db_version="checkm2", check=False)
            rep = out / "quality_report.tsv"
            if rep.exists():
                summary.append(("CheckM2", rep))

        # Faj/virüs -> CheckV
        if util.count_fasta_seqs(viral) > 0:
            cdb = glob.glob(f"{dbp}/checkv/checkv-db-*")
            out = self.out_dir / "checkv"
            cmd = ["checkv", "end_to_end", str(viral), str(out), "-t", str(t)]
            if cdb:
                cmd += ["-d", cdb[0]]
            r.run("checkv", cmd, conda_env=E["checkv"],
                  version_cmd=["checkv", "--version"], db_version="checkv", check=False)
            rep = out / "quality_summary.tsv"
            if rep.exists():
                summary.append(("CheckV", rep))

        # Birleşik özet
        with open(self.outputs()[0], "w") as out:
            if not summary:
                out.write("# completeness çalıştırılamadı (sınıflanmış contig yok)\n")
            for tool, rep in summary:
                out.write(f"# === {tool} ===\n")
                out.write(Path(rep).read_text())
                out.write("\n")
