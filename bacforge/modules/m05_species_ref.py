"""M05 -- Species & Reference Identification (SPECIES-AGNOSTIC, BLAST'li)
Organizmayi VERIDEN kesfeder (sabit tur YOK):
  1) En uzun contig'i NCBI nt'ye remote blastn (timeout'lu)
  2) Timeout/basarisiz -> barrnap ile en tam 16S -> 16S remote blastn (fallback)
  3) BLAST top-hit'lerinden organizma (consensus sscinames)
Sonra kesfedilen turun tam genomlarini `datasets` ile ceker -> FastANI -> GERCEK closest-5.
KATI: hicbir sey uydurulmaz. Kimlik bulunamazsa WARNING; referans cekilemezse closest-5 NA.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import signal
import subprocess
import zipfile
from collections import Counter
from pathlib import Path

from .base import Module
from .. import util

BLAST_FMT = "6 sacc pident length evalue bitscore staxids sscinames stitle"


class SpeciesReferenceIdentificationModule(Module):
    number = "05"
    name = "species_reference_identification"
    folder = "M05_SPECIES_REFERENCE_IDENTIFICATION"
    enabled_key = "ani"

    BLAST_TIMEOUT = 3600              # DAIMA BLAST + sonucu bekle (16S fallback YOK). Uzun bekleme guvenlik siniri.
    MAX_REFS_FETCH = 60               # kesfedilen turden cekilecek tam genom ust siniri
    BLAST_MAX_CONTIG_BP = 1_000_000   # remote nt BLAST icin contig ust siniri (>1Mb remote'ta takilir/timeout)
    CLOSEST_N = 10                    # akrabalik haritasi icin en yakin TEKIL genom sayisi (5-10)

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "closest_5_strains.json"]

    # ---------- yardimcilar ----------
    def _conda(self, env, cmd):
        return ["conda", "run", "-n", env] + cmd

    @staticmethod
    def _killpg(proc):
        """Bir surecin TUM surec-grubunu oldur (conda-run + torunlari) -> orphan kalmaz."""
        if proc is None:
            return
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass
        try:
            proc.wait(timeout=10)
        except Exception:
            pass

    @staticmethod
    def _safe_name(name: str) -> str:
        """Contig id'sini dosya-guvenli hale getir (yol ayraclari vs. -> _)."""
        return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "contig"

    def _export_per_contig(self, genome: Path, contigs_dir: Path, std_dir: Path) -> list[tuple[str, int]]:
        """Her contig'i AYRI FASTA olarak `contigs_dir` altina yazar (manuel BLAST icin) +
        contig uzunluk tablosunu (`contig_lengths.tsv`) std_dir'e yazar. Toplu genome.fasta'ya DOKUNULMAZ.
        Doner: (contig_id, uzunluk) listesi, uzunluga gore azalan."""
        seqs = util.read_fasta(genome)
        lengths = sorted(((n, len(s)) for n, s in seqs.items()), key=lambda x: x[1], reverse=True)
        contigs_dir.mkdir(parents=True, exist_ok=True)
        used = {}
        index_rows = []
        for rank, (name, L) in enumerate(lengths, 1):
            base = self._safe_name(name)
            fn = base
            if fn in used:  # cakisma olursa benzersizlestir
                used[fn] += 1
                fn = f"{base}__{used[fn]}"
            else:
                used[fn] = 0
            fa = contigs_dir / f"{fn}.fasta"
            fa.write_text(f">{name}\n{seqs[name]}\n", encoding="utf-8")
            index_rows.append((rank, name, L, fa.name, "yes" if L < self.BLAST_MAX_CONTIG_BP else "no (>1Mb)"))
        # uzunluk tablosu (std_dir'e; rapora tasinabilir)
        with open(std_dir / "contig_lengths.tsv", "w", encoding="utf-8") as fh:
            fh.write("Rank\tContig\tLength_bp\tPer_contig_fasta\tBLAST_eligible_lt1Mb\n")
            for rank, name, L, fn, elig in index_rows:
                fh.write(f"{rank}\t{name}\t{L}\t{fn}\t{elig}\n")
        # klasore kullanim notu
        (contigs_dir / "README.txt").write_text(
            "Her contig ayri FASTA (manuel yukleme / web BLAST icin). Toplu (ana) dizi: "
            f"../{genome.name} (degistirilmez).\n"
            "Ornek manuel BLAST:\n"
            "  conda run -n ali-blast blastn -query <contig>.fasta -db nt -remote \\\n"
            "    -outfmt '6 sacc pident length evalue bitscore staxids sscinames stitle' -max_target_seqs 5\n"
            f"Not: >{self.BLAST_MAX_CONTIG_BP} bp contig'ler remote nt BLAST'ta takilabilir.\n",
            encoding="utf-8")
        return [(name, L) for _, name, L, _, _ in index_rows]

    def _pick_blast_contig(self, genome: Path, out: Path) -> tuple[str | None, int, bool]:
        """Remote BLAST icin contig sec: 1Mb ALTINDAKI en uzun contig (>1Mb remote'ta takilir).
        1Mb alti hic yoksa en uzun contig'in ilk 1Mb'lik parcasini kullan (truncate).
        Doner: (contig_id, kullanilan_uzunluk, truncated?)."""
        seqs = util.read_fasta(genome)
        if not seqs:
            return None, 0, False
        eligible = {n: s for n, s in seqs.items() if len(s) < self.BLAST_MAX_CONTIG_BP}
        if eligible:
            name = max(eligible, key=lambda k: len(eligible[k]))
            out.write_text(f">{name}\n{seqs[name]}\n", encoding="utf-8")
            return name, len(seqs[name]), False
        # hepsi >=1Mb: en uzunu al, ilk 1Mb'i BLAST'la
        name = max(seqs, key=lambda k: len(seqs[k]))
        frag = seqs[name][: self.BLAST_MAX_CONTIG_BP]
        out.write_text(f">{name}_first{self.BLAST_MAX_CONTIG_BP}bp\n{frag}\n", encoding="utf-8")
        return name, len(frag), True

    def _best_16s(self, genome: Path, work: Path, log_prefix: str) -> Path | None:
        """barrnap ile TUM 16S kopyalarini bul, en UZUN (en tam ~1500bp) olani dondur."""
        gff = work / "rrna.gff"
        try:
            with open(gff, "w") as fh:
                p = subprocess.run(self._conda(util.ENV.get("comparative", "base"), ["barrnap", str(genome)]),
                                   stdout=fh, stderr=subprocess.DEVNULL, timeout=600)
            if p.returncode != 0:
                return None
        except Exception:
            return None
        seqs = util.read_fasta(genome)
        best = None  # (len, seq, header)
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        for line in gff.read_text(encoding="utf-8").splitlines():
            if line.startswith("#") or "16S" not in line:
                continue
            pp = line.split("\t")
            if len(pp) < 7:
                continue
            c, s, e, strand = pp[0], int(pp[3]), int(pp[4]), pp[6]
            if c not in seqs:
                continue
            sub = seqs[c][s - 1:e]
            if strand == "-":
                sub = "".join(comp.get(b, "N") for b in reversed(sub))
            if best is None or len(sub) > best[0]:
                best = (len(sub), sub, f"16S_{c}_{s}_{e}")
        if not best:
            return None
        out = work / "best_16s.fa"
        out.write_text(f">{best[2]}\n{best[1]}\n", encoding="utf-8")
        return out

    def _blast(self, query_fa: Path, out_tsv: Path, timeout: int, runner) -> list | None:
        """remote blastn; timeout/hata -> None. Basari -> parse edilmis satirlar.
        KRITIK: kendi surec-grubunda baslatilir (start_new_session); timeout'ta TUM grup oldurulur
        (os.killpg) -> `conda run` torunu blastn ORPHAN kalmaz (eski hata)."""
        cmd = self._conda(util.ENV.get("blast", "base"),
                          ["blastn", "-query", str(query_fa), "-db", "nt", "-remote",
                           "-outfmt", BLAST_FMT, "-max_target_seqs", "10"])
        proc = None
        try:
            with open(out_tsv, "w") as fh:
                proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.DEVNULL,
                                        start_new_session=True)
                try:
                    rc = proc.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    self._killpg(proc)
                    return None
            if rc != 0:
                return None
        except Exception:
            self._killpg(proc)
            return None
        rows = []
        for line in out_tsv.read_text(encoding="utf-8").splitlines():
            p = line.split("\t")
            if len(p) >= 8:
                rows.append({"sacc": p[0], "pident": p[1], "sciname": p[6], "title": p[7]})
        return rows or None

    def _datasets(self, args: list, timeout: int, capture=False):
        """datasets base env'de kurulu; PATH'te varsa dogrudan, yoksa `conda run -n base`.
        Surec-grubu ile calistirilir; timeout'ta tum grup oldurulur (orphan kalmaz)."""
        for c in ([args] if shutil.which(args[0]) else []) + [["conda", "run", "-n", "base"] + args]:
            proc = None
            try:
                proc = subprocess.Popen(c, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                        text=True, start_new_session=True)
                try:
                    out, _ = proc.communicate(timeout=timeout)
                except subprocess.TimeoutExpired:
                    self._killpg(proc)
                    continue
                if proc.returncode == 0:
                    return out if capture else True
            except Exception:
                self._killpg(proc)
                continue
        return None if capture else False

    def _fetch_species_genomes(self, organism: str, cache: Path, runner) -> list[Path]:
        """datasets summary->accession->download ile kesfedilen turun tam genomlarini ceker (cap).
        Cache'de varsa yeniden cekmez. ('--limit' bu surumde yok; accession listesiyle sinirlanir.)"""
        cache.mkdir(parents=True, exist_ok=True)
        existing = [p for p in cache.rglob("*.fna") if not self._is_annotation_copy(p)]
        if existing:
            return sorted(existing)
        # 1) accession listesi (summary --as-json-lines)
        out = self._datasets(["datasets", "summary", "genome", "taxon", organism,
                              "--assembly-level", "complete", "--as-json-lines"],
                             timeout=300, capture=True)
        if not out:
            return []
        import json as _json
        accs = []
        for line in out.splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                a = _json.loads(line).get("accession")
                if a:
                    accs.append(a)
            except Exception:
                pass
            if len(accs) >= self.MAX_REFS_FETCH:
                break
        if not accs:
            return []
        # 2) accession ile indir
        zip_path = cache / "refs.zip"
        ok = self._datasets(["datasets", "download", "genome", "accession", *accs,
                            "--include", "genome", "--filename", str(zip_path)], timeout=1800)
        if not ok or not zip_path.exists():
            return []
        try:
            with zipfile.ZipFile(zip_path) as zf:
                zf.extractall(cache / "unz")
        except Exception:
            return []
        found = sorted((cache / "unz").rglob("*.fna")) + sorted((cache / "unz").rglob("*.fa"))
        return [p for p in found if not self._is_annotation_copy(p)]

    @staticmethod
    def _slug(name: str) -> str:
        return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()

    @staticmethod
    def _is_annotation_copy(p: Path) -> bool:
        """Bakta/prokka anotasyon ciktisi olan kopya .fna (ornegin .../<acc>_bakta/ref.fna).
        Bunlar closest listesinde ANA genomla AYNI olup TEKRAR yaprak yaratir -> dislanir."""
        s = str(p)
        return "_bakta/" in s or "_prokka/" in s or p.name == "ref.fna"

    @staticmethod
    def _accession(path: str) -> str:
        """Yol icinden gercek assembly accession'ini (GCF_/GCA_) cikar; yoksa dosya adi."""
        m = re.search(r"GC[FA]_\d+\.\d+", str(path))
        return m.group(0) if m else Path(path).stem

    def _parse_rank_fastani(self, ani_out: Path, organism: str, top_n: int) -> list[dict]:
        """FastANI ciktisini oku -> anotasyon kopyalarini at -> accession'a gore DEDUP ->
        ANI'ye gore sirala -> ilk top_n TEKIL genom. Ayni genomun tekrari YOK (dogru akrabalik agaci)."""
        best: dict[str, tuple] = {}   # accession -> (ani, cov, fasta_path)
        for line in ani_out.read_text(encoding="utf-8").splitlines():
            p = line.split("\t")
            if len(p) < 5:
                continue
            refpath = p[1]
            if self._is_annotation_copy(Path(refpath)):
                continue
            try:
                ani, m, tot = float(p[2]), int(p[3]), int(p[4])
            except ValueError:
                continue
            acc = self._accession(refpath)
            cov = round(100.0 * m / tot, 2) if tot else 0.0
            if acc not in best or ani > best[acc][0]:
                best[acc] = (ani, cov, refpath)
        ranked = sorted(best.items(), key=lambda kv: kv[1][0], reverse=True)[:top_n]
        return [{"rank": i, "organism": organism, "strain": acc, "assembly_accession": acc,
                 "ani_percent": round(ani, 4), "query_coverage": cov, "fasta_path": rp}
                for i, (acc, (ani, cov, rp)) in enumerate(ranked, 1)]

    # ---------- ana ----------
    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        work = self.sub_dir("02_work")
        r = self.ctx.runner
        t = util.threads(self.ctx)

        # 0) Her contig'i ayri FASTA + uzunluk tablosu -> ANA contig dosyasinin (genome.fasta) yaninda
        #    `contigs/` altina (manuel yukleme/web BLAST icin kolay erisim). Toplu genome.fasta'ya DOKUNULMAZ.
        contigs_dir = genome.parent / "contigs"
        contig_lens = self._export_per_contig(genome, contigs_dir, contigs_dir)
        # rapor plumbing icin uzunluk tablosunun bir kopyasi M05 std_dir'de de dursun
        try:
            (std_dir / "contig_lengths.tsv").write_text(
                (contigs_dir / "contig_lengths.tsv").read_text(encoding="utf-8"), encoding="utf-8")
        except Exception:
            pass

        # 1) BLAST kimlik: 1Mb ALTINDAKI en uzun contig'i NCBI nt'ye blastn + sonucu BEKLE (16S fallback YOK).
        #    (>1Mb contig remote nt BLAST'ta takilir/timeout -> ya 1Mb-alti sec, ya ilk 1Mb'i BLAST'la)
        query_fa = work / "blast_query_contig.fa"
        blast_contig, blast_bp, truncated = self._pick_blast_contig(genome, query_fa)
        method = "contig_blastn_lt1Mb" + ("_truncated" if truncated else "")
        hits = self._blast(query_fa, work / "blast_query.tsv", self.BLAST_TIMEOUT, r) if blast_contig else None

        # kraken2 (M02) capraz-kontrol
        kraken_species = util.resolve_species(self.ctx)  # resume'da M02 dosyasından da çözer

        organism = None
        if hits:
            names = [h["sciname"].strip() for h in hits if h.get("sciname") and h["sciname"] != "N/A"]
            if names:
                organism = Counter(names).most_common(1)[0][0]
        if not organism:
            organism = kraken_species  # BLAST vermezse M02'ye dus (yine de veri-temelli)

        # 2) closest-5: kesfedilen turun genomlarini cek -> FastANI
        strains = []
        closest_note = None
        if organism:
            cache = Path(self.ctx.config["paths"]["db"]) / "references" / self._slug(organism)
            refs = self._fetch_species_genomes(organism, cache, r)
            if refs:
                reflist = work / "ref_list.txt"
                reflist.write_text("\n".join(str(x) for x in refs) + "\n", encoding="utf-8")
                ani_out = work / "fastani.txt"
                prov = r.run("fastani", ["fastANI", "-q", str(genome), "--rl", str(reflist),
                                         "-o", str(ani_out), "-t", str(t)],
                             conda_env=util.ENV.get("species", "base"),
                             version_cmd=["fastANI", "--version"], check=False)
                if prov.get("exit_code") == 0 and ani_out.exists():
                    # DEDUP: accession'a gore tekil, anotasyon kopyalari haric -> ilk 10 tekil genom
                    strains = self._parse_rank_fastani(ani_out, organism, self.CLOSEST_N)
                    if not strains:
                        closest_note = "FastANI eslesme vermedi."
                else:
                    closest_note = f"FastANI hata (exit {prov.get('exit_code')})."
            else:
                closest_note = f"'{organism}' icin referans genom cekilemedi (datasets)."
        else:
            closest_note = "Organizma kimligi belirlenemedi (BLAST+kraken bos)."

        # yaz: closest_5 (dashboard/rapor tablosu) + closest_10 (akrabalik haritasi) — hepsi TEKIL genom
        def _write_closest(items, stem):
            std_dir.joinpath(f"{stem}.json").write_text(
                json.dumps(items, indent=2, ensure_ascii=False), encoding="utf-8")
            with open(std_dir / f"{stem}.tsv", "w", encoding="utf-8") as fh:
                fh.write("Rank\tOrganism\tAccession\tANI_percent\tQuery_coverage\n")
                for s in items:
                    fh.write(f"{s['rank']}\t{s['organism']}\t{s['assembly_accession']}"
                             f"\t{s['ani_percent']}\t{s['query_coverage']}\n")
        _write_closest(strains[:5], "closest_5_strains")
        _write_closest(strains, "closest_10_strains")
        std_dir.joinpath("species_identification.json").write_text(json.dumps({
            "organism": organism, "identification_method": method,
            "blast_query_contig": blast_contig, "blast_query_bp": blast_bp,
            "blast_query_truncated": truncated,
            "blast_top_hits": hits[:5] if hits else [], "kraken2_species": kraken_species,
            "closest_reference": strains[0] if strains else None,
            "contig_count": len(contig_lens),
            "per_contig_fasta_dir": str(contigs_dir),
        }, indent=2, ensure_ascii=False), encoding="utf-8")

        # tur bilgisini downstream'e tasi
        if organism:
            self.ctx.detection["ncbi_species"] = organism

        # KATI durum: kimlik yoksa WARNING; kimlik var closest yoksa WARNING(NA-note); ikisi de varsa PASS
        if organism and strains:
            status = "PASS"
            warns = []
        elif organism:
            status = "WARNING"
            warns = [closest_note or "Closest-5 hesaplanamadi."]
        else:
            status = "WARNING"
            warns = [closest_note or "Kimlik belirlenemedi."]
        details = {}
        if closest_note:
            details["note"] = closest_note
        details["blast_query_contig"] = blast_contig
        details["blast_query_bp"] = blast_bp
        if truncated:
            details["blast_query_truncated"] = f"contig >1Mb; ilk {self.BLAST_MAX_CONTIG_BP} bp BLAST'landi"
        details["per_contig_fasta_dir"] = str(contigs_dir)
        self.write_summary(status=status,
                           statistics={"organism": organism, "identification_method": method,
                                       "closest_count": len(strains),
                                       "contig_count": len(contig_lens),
                                       "blast_query_bp": blast_bp},
                           warnings=warns,
                           details=details)
