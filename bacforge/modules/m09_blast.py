"""Modül 09 — BLAST (viral contig'ler -> NCBI virüs NT) + AKILLI COVERAGE.
INPUT : 08_Taxonomy/viral.fasta
OUTPUT: 09_BLAST/blast.tsv (ham) + best_hits.tsv (her contig: bp + query/subject/reciprocal coverage)
Hesap: HSP koordinatlarından birleştirilmiş aralıkla gerçek kaplanma (overlap düzeltmeli).
Gerekçe: docs/literature/09_blast.md (orijinal prompt #10 Akıllı BLAST Filtreleme)
"""
from __future__ import annotations

from pathlib import Path

from .base import Module
from .. import util

# koordinatlar + qlen/slen ile coverage hesaplanır
OUTFMT = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen stitle"


def _merged_len(intervals):
    """Örtüşen aralıkları birleştirip toplam kaplanan bp döndür."""
    if not intervals:
        return 0
    ivs = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total = 0
    cur_s, cur_e = ivs[0]
    for s, e in ivs[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    total += cur_e - cur_s + 1
    return total


class BlastModule(Module):
    number = "09"
    name = "blast"
    folder = "09_BLAST"
    enabled_key = "blast"

    def _viral(self):
        return Path(self.ctx.run_dir) / "08_Taxonomy" / "viral.fasta"

    def inputs(self):
        return [self._viral()]

    def outputs(self):
        return [self.out_dir / "best_hits.tsv"]

    def run(self):
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        cfg = self.ctx.config.get("tools", {}).get("blastn", {})
        dbp = self.ctx.config["paths"]["db"]
        db = f"{dbp}/{cfg.get('db', 'blast_viral/ref_viruses_rep_genomes')}"
        viral = self._viral()

        if util.count_fasta_seqs(viral) == 0:
            self._write_empty()
            return

        blast_raw = self.out_dir / "blast.tsv"
        r.run("blastn", ["blastn", "-query", str(viral), "-db", db,
                        "-out", str(blast_raw), "-outfmt", OUTFMT,
                        "-evalue", str(cfg.get("evalue", 1e-10)),
                        "-max_target_seqs", str(cfg.get("max_target_seqs", 5)),
                        "-num_threads", str(t)],
              conda_env=E["blast"], version_cmd=["blastn", "-version"],
              db_version="ref_viruses_rep_genomes", check=False)
        self._best_hits(blast_raw)

    def _write_empty(self):
        cols = "contig\tcontig_len\taligned_bp\tquery_cov%\tref_len\tref_cov%\treciprocal%\tidentity%\tbest_hit\tconfident\n"
        (self.out_dir / "best_hits.tsv").write_text(cols + "# viral contig yok\n")
        (self.out_dir / "blast.tsv").write_text("")

    def _best_hits(self, blast_raw: Path):
        # query -> subject -> {hsps, qlen, slen, stitle, bitsum}
        data = {}
        if blast_raw.exists():
            with open(blast_raw) as fh:
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 13:
                        continue
                    q, s, pid, alen, qs, qe, ss, se, ev, bit, qlen, slen = c[:12]
                    stitle = c[12]
                    d = data.setdefault(q, {}).setdefault(s, {
                        "q": [], "s": [], "pid": [], "alen": [], "bit": 0.0,
                        "qlen": int(qlen), "slen": int(slen), "stitle": stitle})
                    d["q"].append((int(qs), int(qe)))
                    d["s"].append((int(ss), int(se)))
                    d["pid"].append(float(pid))
                    d["alen"].append(int(alen))
                    d["bit"] += float(bit)

        header = ("contig\tcontig_len\taligned_bp\tquery_cov%\tref_len\tref_cov%\t"
                  "reciprocal%\tidentity%\tbest_hit\tconfident")
        lines = [header]
        for q, subs in data.items():
            # en iyi subject = toplam bitscore en yüksek
            best_s = max(subs, key=lambda s: subs[s]["bit"])
            d = subs[best_s]
            qcov_bp = _merged_len(d["q"])
            scov_bp = _merged_len(d["s"])
            qlen, slen = d["qlen"], d["slen"]
            qcov = 100 * qcov_bp / qlen if qlen else 0
            scov = 100 * scov_bp / slen if slen else 0
            recip = min(qcov, scov)
            # uzunluk-ağırlıklı identity
            tot = sum(d["alen"]) or 1
            ident = sum(p * a for p, a in zip(d["pid"], d["alen"])) / tot
            confident = "yes" if (ident >= 70 and qcov >= 50) else "low"
            lines.append(f"{q}\t{qlen}\t{qcov_bp}\t{qcov:.1f}\t{slen}\t{scov:.1f}\t"
                         f"{recip:.1f}\t{ident:.1f}\t{d['stitle']}\t{confident}")
        if len(lines) == 1:
            lines.append("# hit bulunamadı (virüs NT'de eşleşme yok)")
        (self.out_dir / "best_hits.tsv").write_text("\n".join(lines) + "\n")
