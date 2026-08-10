"""Modül 10 — Annotation (içerik-farkında otomatik routing).
INPUT : 08_Taxonomy/{bacterial,viral,plasmid}.fasta
OUTPUT: 10_Annotation/bakta/... (bakteri/plazmid) ve/veya pharokka/... (faj)
Kullanıcı seçimi YOK: hangi sınıf varsa o araç çalışır.
"""
from __future__ import annotations

from pathlib import Path

from .base import Module
from .. import util


class AnnotationModule(Module):
    number = "10"
    name = "annotation"
    folder = "10_Annotation"
    enabled_key = "annotation"

    def _cls(self, name):
        return Path(self.ctx.run_dir) / "08_Taxonomy" / name

    def inputs(self):
        return [self._cls("classification.tsv")]

    def outputs(self):
        return [self.out_dir / "annotation.done"]

    def run(self):
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        dbp = self.ctx.config["paths"]["db"]

        bacterial = self._cls("bacterial.fasta")
        plasmid = self._cls("plasmid.fasta")
        viral = self._cls("viral.fasta")
        ran = []

        # Bakteri + plazmid -> Bakta (Diamond -11 segfault'a karşı dayanıklı)
        for tag, fa in (("chromosome", bacterial), ("plasmid", plasmid)):
            if util.count_fasta_seqs(fa) > 0:
                if self._bakta(tag, fa, dbp, t, r, E):
                    ran.append(f"bakta:{tag}")
                else:
                    (self.out_dir / f"bakta_{tag}.FAILED").write_text(
                        "Bakta --skip-sorf ile de başarısız (bkz 19_Logs)\n")

        # Faj -> Pharokka (CDS-id çakışmasını önlemek için sabit-genişlikli ad)
        if util.count_fasta_seqs(viral) > 0:
            if self._pharokka(viral, dbp, t, r, E):
                ran.append("pharokka:virus")
            else:
                (self.out_dir / "pharokka.FAILED").write_text("Pharokka başarısız (bkz 19_Logs)\n")

        if not ran:
            (self.out_dir / "WARN.txt").write_text(
                "Sınıflanmış contig yok; annotation çalıştırılamadı\n")
        (self.out_dir / "annotation.done").write_text("\n".join(ran) + "\n")

    def _pharokka(self, viral, dbp, t, r, E) -> bool:
        """Viral contig'leri sabit-genişlikli güvenli adla (ctgNNNNN) Pharokka'ya ver -> CDS id çakışması olmaz.
        name_map.tsv: safe<TAB>original (dossier'da geri eşlemek için)."""
        seqs = util.read_fasta(viral)
        namemap, renamed = {}, {}
        for i, (name, seq) in enumerate(sorted(seqs.items()), 1):
            safe = f"ctg{i:05d}"
            namemap[safe] = name
            renamed[safe] = seq
        rn = self.out_dir / "viral_renamed.fasta"
        util.write_fasta(renamed, rn)
        (self.out_dir / "pharokka_namemap.tsv").write_text(
            "safe\toriginal\n" + "\n".join(f"{s}\t{o}" for s, o in namemap.items()) + "\n")
        prov = r.run("pharokka", ["pharokka.py", "-i", str(rn),
                                  "-o", str(self.out_dir / "pharokka"),
                                  "-d", f"{dbp}/pharokka", "-t", str(t), "-f"],
                     conda_env=E["pharokka"], version_cmd=["pharokka.py", "--version"],
                     db_version="pharokka", check=False)
        return prov["exit_code"] == 0

    def _bakta(self, tag, fa, dbp, t, r, E) -> bool:
        """Bakta çalıştır; Diamond -11 (sORF segfault, bakta#424) olursa --skip-sorf ile yeniden dene."""
        out = self.out_dir / f"bakta_{tag}"

        def cmd(skip_sorf):
            c = ["bakta", "--db", f"{dbp}/bakta/db-light", "--output", str(out),
                 "--prefix", tag, "--threads", str(t), "--force"]
            if skip_sorf:
                c.append("--skip-sorf")
            c.append(str(fa))
            return c

        prov = r.run(f"bakta_{tag}", cmd(False), conda_env=E["bakta"],
                     version_cmd=["bakta", "--version"], db_version="bakta-light", check=False)
        if prov["exit_code"] == 0:
            return True
        # Bilinen Diamond -11 segfault -> sORF atlayarak yeniden dene (CDS/tRNA/rRNA/CRISPR/AMR korunur)
        prov = r.run(f"bakta_{tag}_skipsorf", cmd(True), conda_env=E["bakta"],
                     db_version="bakta-light", check=False)
        return prov["exit_code"] == 0
