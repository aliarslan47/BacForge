"""01_Input: otomatik platform/okuma-tipi/pair tespiti. Sadece stdlib (env'siz çalışır)."""
from __future__ import annotations

import gzip
from pathlib import Path

SIGNAL_EXT = {".fast5", ".pod5"}
FASTA_EXT = {".fasta", ".fa", ".fna"}
FASTQ_EXT = {".fastq", ".fq"}


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _inner_ext(path: Path) -> str:
    """'.gz' ise içteki uzantıyı döndür (reads.fastq.gz -> .fastq)."""
    if path.suffix.lower() == ".gz":
        return Path(path.stem).suffix.lower()
    return path.suffix.lower()


def _sample_fastq(path, max_reads=5000):
    lengths, quals = [], []
    with _open(path) as fh:
        for i, line in enumerate(fh):
            m = i % 4
            if m == 1:
                lengths.append(len(line.strip()))
            elif m == 3:
                q = line.strip()
                if q:
                    quals.append(sum(ord(c) - 33 for c in q) / len(q))
            if len(lengths) >= max_reads and m == 3:
                break
    return lengths, quals


def _n50(lengths):
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half, acc = sum(ordered) / 2, 0
    for length in ordered:
        acc += length
        if acc >= half:
            return length
    return ordered[-1]


def detect_platform(input_path, config: dict | None = None) -> dict:
    cfg = (config or {}).get("input_detection", {}) if config else {}
    short_max = cfg.get("short_read_max_len", 350)
    hifi_q = cfg.get("hifi_min_qual", 25)
    pair_tokens = cfg.get("paired_tokens", ["_R1", "_R2", "_1.", "_2."])

    p = Path(input_path)
    files = sorted(f for f in p.iterdir() if f.is_file()) if p.is_dir() else [p]
    result = {"input": str(p), "files": [str(f) for f in files]}
    exts = {_inner_ext(f) for f in files}

    # 1) Ham sinyal -> basecalling gerekli
    if exts & SIGNAL_EXT:
        result.update(platform="ONT", read_type="long", data="signal",
                      note="FAST5/POD5 sinyal: basecalling (dorado) gerekli")
        return result

    fastqs = [f for f in files if _inner_ext(f) in FASTQ_EXT]
    fastas = [f for f in files if _inner_ext(f) in FASTA_EXT]

    # 2) Sadece FASTA -> önceden assemble edilmiş
    if not fastqs and fastas:
        result.update(platform="unknown", read_type="assembly_input", data="fasta",
                      note="FASTA girdi: assembly atlanır, doğrudan QC/annotation")
        return result
    if not fastqs:
        result.update(platform="unknown", note="Tanınan okuma dosyası bulunamadı")
        return result

    # 3) FASTQ -> uzunluk/kalite ile platform çıkarımı
    lengths, quals = _sample_fastq(fastqs[0])
    mean_len = sum(lengths) / len(lengths) if lengths else 0
    n50 = _n50(lengths)
    mean_q = sum(quals) / len(quals) if quals else 0
    paired = len(fastqs) >= 2 and any(t in f.name for f in fastqs for t in pair_tokens)

    if mean_len <= short_max:
        platform, read_type = "Illumina", "short"
    elif mean_q >= hifi_q and n50 >= 5000:
        platform, read_type = "PacBio_HiFi", "long"
    else:
        platform, read_type = "ONT", "long"

    result.update(
        platform=platform, read_type=read_type, data="reads", paired=paired,
        mean_read_length=round(mean_len, 1), n50=n50, mean_quality=round(mean_q, 1),
        n_reads_sampled=len(lengths),
        note="Eşikler config.input_detection'dan; FAZ 3'te literatürle sabitlenecek",
    )
    return result
