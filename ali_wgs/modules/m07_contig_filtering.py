"""Modül 07 — Contig Filtering (LİTERATÜR-TEMELLİ: uzunluk + coverage).
INPUT : 05_Assembly/assembly.fasta (+ assembly_info.txt varsa coverage için)
OUTPUT: 07_Contig_Filtering/contigs.filtered.fasta + filter_report.tsv
Eşikler: config.contig_filter (gerekçe: docs/literature/07_contig_filtering.md)
"""
from __future__ import annotations

import statistics
from pathlib import Path

from .base import Module
from .. import util


class ContigFilteringModule(Module):
    number = "07"
    name = "contig_filtering"
    folder = "07_Contig_Filtering"
    enabled_key = None  # her zaman çalışır

    def inputs(self):
        return [util.assembly_fasta(self.ctx)]

    def outputs(self):
        return [util.filtered_contigs(self.ctx)]

    def _coverages(self) -> dict:
        """Flye assembly_info.txt'ten contig -> coverage (varsa)."""
        info = Path(self.ctx.run_dir) / "05_Assembly" / "assembly_info.txt"
        cov = {}
        if info.exists():
            with open(info) as fh:
                for i, line in enumerate(fh):
                    if i == 0:
                        continue
                    c = line.split("\t")
                    if len(c) >= 3:
                        try:
                            cov[c[0]] = float(c[2])
                        except ValueError:
                            pass
        return cov

    def run(self):
        cfg = self.ctx.config.get("contig_filter", {})
        min_len = int(cfg.get("min_length", 1000))
        min_cov_abs = float(cfg.get("min_coverage_abs", 3))
        min_cov_frac = float(cfg.get("min_coverage_frac", 0.10))
        apply_cov = bool(cfg.get("apply_coverage", True))

        seqs = util.read_fasta(util.assembly_fasta(self.ctx))
        cov = self._coverages() if apply_cov else {}
        use_cov = apply_cov and len(cov) > 0
        cov_floor = 0.0
        if use_cov:
            median_cov = statistics.median(cov.values())
            cov_floor = max(min_cov_abs, min_cov_frac * median_cov)

        kept, report = {}, ["contig\tlength\tcoverage\tkept\treason"]
        for name, seq in seqs.items():
            L = len(seq)
            c = cov.get(name, None)
            reason = "ok"
            keep = True
            if L < min_len:
                keep, reason = False, f"length<{min_len}"
            elif use_cov and c is not None and c < cov_floor:
                keep, reason = False, f"coverage<{cov_floor:.1f}"
            if keep:
                kept[name] = seq
            report.append(f"{name}\t{L}\t{'' if c is None else c}\t{keep}\t{reason}")

        util.write_fasta(kept, self.outputs()[0])
        (self.out_dir / "filter_report.tsv").write_text("\n".join(report) + "\n")
        (self.out_dir / "filter_params.txt").write_text(
            f"min_length={min_len}\nmin_coverage_abs={min_cov_abs}\n"
            f"min_coverage_frac={min_cov_frac}\ncoverage_floor_used="
            f"{cov_floor if use_cov else 'N/A (coverage yok)'}\n"
            f"kept={len(kept)}/{len(seqs)}\n")
