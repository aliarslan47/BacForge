"""Modüllerin paylaştığı yardımcılar: kaynak, okuma çözümleme, platform, FASTA işlemleri."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

FASTQ_INNER = {".fastq", ".fq"}
FASTA_INNER = {".fasta", ".fa", ".fna"}

# Tool -> izole conda env (setup/setup_envs.sh ile birebir). ToolRunner conda_env'i kullanır.
ENV = {
    "core": "bacforge",           # seqkit, minimap2, samtools, weasyprint
    "ont_qc": "ali-ont-qc",      # nanoplot, chopper, filtlong, rasusa, seqkit
    "ill_qc": "ali-illumina-qc", # fastqc, fastp
    "illumina_qc": "ali-illumina-qc",
    "flye": "ali-flye",
    "medaka": "ali-medaka",
    "asm_sr": "ali-assembly-sr", # spades, skesa, unicycler
    "hifiasm": "ali-hifiasm",
    "quast": "ali-quast",
    "genomad": "ali-genomad",
    "bakta": "ali-bakta",
    "pharokka": "ali-pharokka",
    "checkm2": "ali-checkm2",
    "checkv": "ali-checkv",
    "bacphlip": "ali-bacphlip",
    "amrfinder": "ali-amrfinder",
    "blast": "ali-blast",
    "kleborate": "ali-kleborate",
    "mobsuite": "ali-mobsuite",
    "clinker": "ali-clinker",
    "iqtree": "ali-phylogeny",
    "qc": "ali-qc",
    "species": "ali-species",
    "typing": "ali-typing",
    "virulence": "ali-virulence",
    "phylogeny": "ali-phylogeny",
    "comparative": "ali-comparative",
    "defense": "ali-defense",
    "mge": "ali-mge",
}


def db_path(ctx, key: str) -> str:
    return str((ctx.config.get("databases", {}).get(key, {}) or {}).get("path", ""))


def count_fasta_seqs(path) -> int:
    if not Path(path).exists():
        return 0
    n = 0
    with fopen(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def threads(ctx) -> int:
    return int(ctx.resources.get("threads", 4))


def memory_gb(ctx) -> int:
    return int(ctx.resources.get("memory_gb", 8))


def load_detection(ctx) -> dict:
    """ctx.detection doluysa onu, değilse 01_Input/platform.json'u oku (resume için)."""
    if getattr(ctx, "detection", None):
        return ctx.detection
    pj = Path(ctx.run_dir) / "01_Input" / "platform.json"
    if pj.exists():
        with open(pj) as fh:
            ctx.detection = json.load(fh)
            return ctx.detection
    return {}


def platform(ctx) -> str:
    return load_detection(ctx).get("platform", "unknown")


def is_paired(ctx) -> bool:
    return bool(load_detection(ctx).get("paired", False))


def is_fasta_input(ctx) -> bool:
    return load_detection(ctx).get("data") == "fasta"


def _inner_ext(p: Path) -> str:
    return Path(p.stem).suffix.lower() if p.suffix.lower() == ".gz" else p.suffix.lower()


def raw_read_files(ctx) -> list[Path]:
    """Tespitten ham okuma dosyaları (FASTQ). Paired ise sıralı [R1, R2]."""
    files = [Path(f) for f in load_detection(ctx).get("files", [])]
    fastqs = sorted(f for f in files if _inner_ext(f) in FASTQ_INNER)
    return fastqs


def raw_fasta_files(ctx) -> list[Path]:
    files = [Path(f) for f in load_detection(ctx).get("files", [])]
    return sorted(f for f in files if _inner_ext(f) in FASTA_INNER)


def reads_for_assembly(ctx) -> list[Path]:
    """Assembly girdisi: filtrelenmiş varsa onu, yoksa ham okumayı döndür."""
    filt = Path(ctx.run_dir) / "04_Filtering"
    single = filt / "filtered.fastq.gz"
    r1, r2 = filt / "filtered_R1.fastq.gz", filt / "filtered_R2.fastq.gz"
    if r1.exists() and r2.exists():
        return [r1, r2]
    if single.exists():
        return [single]
    return raw_read_files(ctx)


def assembly_fasta(ctx) -> Path:
    return Path(ctx.run_dir) / "05_Assembly" / "assembly.fasta"


def filtered_contigs(ctx) -> Path:
    return Path(ctx.run_dir) / "07_Contig_Filtering" / "contigs.filtered.fasta"


def fopen(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def read_fasta(path) -> dict[str, str]:
    """Basit FASTA okuyucu: {header(boşluğa kadar): sequence}."""
    seqs, name, buf = {}, None, []
    with fopen(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(seqs: dict[str, str], path, width: int = 80):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def total_bases_fastq(path, cap=None) -> int:
    total = 0
    with fopen(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                total += len(line.strip())
                if cap and total >= cap:
                    break
    return total
