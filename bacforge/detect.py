"""01_Input: otomatik platform/okuma-tipi/pair tespiti. Sadece stdlib (env'siz çalışır)."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

SIGNAL_EXT = {".fast5", ".pod5"}
FASTA_EXT = {".fasta", ".fa", ".fna"}
FASTQ_EXT = {".fastq", ".fq"}

# ONT kimya varsayilanlari: kimya -> (Flye modu, Medaka model deneme sirasi)
# Medaka denemeleri sirali; her biri medaka_consensus'a eklenen ek argumanlar.
# R10: once --bacteria (modern R10.4.1 bakteri icin en iyi), basarisizsa acik SUP modeli.
ONT_CHEMISTRY_DEFAULTS = {
    "r9": {"flye_mode": "--nano-raw", "medaka_attempts": [["-m", "r941_min_sup_g507"]]},
    "r10": {"flye_mode": "--nano-hq",
            "medaka_attempts": [["--bacteria"], ["-m", "r1041_e82_400bps_sup_v5.0.0"]]},
}


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def detect_ont_chemistry(reads_path, config: dict | None = None, sample_headers: int = 200) -> dict:
    """ONT kimyasini (R9.4.1 vs R10.4.x) tespit et -> Flye modu + Medaka modeli.

    Oncelik sirasi (yuksekten dusuk guvene):
      1) config.tools.ont.chemistry ('r9'|'r10') veya config.tools.medaka.model (acik override)
      2) read basligindaki basecall_model_version_id (dorado/guppy modern okumalar) -- ALTIN sinyal
      3) flow_cell_id: FLO-* urun kodu ya da seri prefiksi (FA..=R9.4.1 MinION)
      4) start_time yili (<2022 -> R9, degilse R10) -- DUSUK guven
      5) varsayilan R10 (modern) -- DUSUK guven
    Sonuc dict: chemistry, flye_mode, medaka_attempts, basis, confidence.
    """
    cfg_tools = (config or {}).get("tools", {}) if config else {}
    cfg_ont = cfg_tools.get("ont", {}) or {}
    cfg_medaka = cfg_tools.get("medaka", {}) or {}

    def _pack(chem, basis, confidence, medaka_attempts=None):
        base = ONT_CHEMISTRY_DEFAULTS.get(chem, ONT_CHEMISTRY_DEFAULTS["r10"])
        return {
            "chemistry": chem,
            "flye_mode": base["flye_mode"],
            "medaka_attempts": medaka_attempts or base["medaka_attempts"],
            "basis": basis,
            "confidence": confidence,
        }

    # Sentinel degerler ('auto','none','default','') override DEGILDIR -> otomatik tespit
    _SENTINELS = {"", "auto", "none", "default", "otomatik"}

    def _is_real(v):
        return v is not None and str(v).strip().lower() not in _SENTINELS

    # 1) Acik override (yalnizca gercek bir deger verildiyse)
    model_ov = cfg_medaka.get("model")
    if _is_real(model_ov):
        chem = str(cfg_ont.get("chemistry", "r10")).lower()
        chem = "r9" if chem.startswith("r9") else "r10"
        return _pack(chem, f"config.tools.medaka.model={model_ov}", "override",
                     medaka_attempts=[["-m", str(model_ov).strip()]])
    if _is_real(cfg_ont.get("chemistry")):
        c = str(cfg_ont["chemistry"]).strip().lower()
        chem = "r9" if c.startswith("r9") else "r10"
        return _pack(chem, f"config.tools.ont.chemistry={cfg_ont['chemistry']}", "override")

    # Read basliklarindan sinyal topla
    headers = []
    try:
        with _open(reads_path) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    headers.append(line.strip())
                if len(headers) >= sample_headers:
                    break
    except OSError:
        return _pack("r10", "okuma basligi okunamadi; varsayilan R10", "low")

    blob = "\n".join(headers)

    # 2) basecall_model_version_id (dorado/guppy) -- en guvenilir
    m = re.search(r"basecall_model_version_id=(\S+)", blob)
    if m:
        mv = m.group(1).lower()
        if "r10" in mv or "e82" in mv:
            return _pack("r10", f"basecall_model_version_id={m.group(1)}", "high")
        if "r9" in mv or "r941" in mv:
            return _pack("r9", f"basecall_model_version_id={m.group(1)}", "high")

    # 3) flow_cell_id: urun kodu ya da seri prefiksi
    fc = re.search(r"flow_cell_id=(\S+)", blob)
    if fc:
        fcid = fc.group(1).upper()
        if fcid.startswith(("FLO-MIN114", "FLO-MIN112", "FLO-PRO114", "FLO-PRO112", "FLO-FLG114")):
            return _pack("r10", f"flow_cell_id={fc.group(1)} (R10 urun kodu)", "high")
        if fcid.startswith(("FLO-MIN106", "FLO-MIN107", "FLO-PRO002", "FLO-FLG001")):
            return _pack("r9", f"flow_cell_id={fc.group(1)} (R9 urun kodu)", "high")
        # Seri prefiksi: MinION R9.4.1 flow-cell serileri "FA" ile baslar
        if re.match(r"^FA[A-Z]\d", fcid):
            return _pack("r9", f"flow_cell_id={fc.group(1)} (FA.. seri -> R9.4.1)", "medium")

    # 4) start_time yili
    st = re.search(r"start_time=(\d{4})-", blob)
    if st:
        year = int(st.group(1))
        if year < 2022:
            return _pack("r9", f"start_time yili {year} (<2022 -> R9)", "low")
        return _pack("r10", f"start_time yili {year} (>=2022 -> R10)", "low")

    # 5) Varsayilan
    return _pack("r10", "belirleyici sinyal yok; varsayilan modern R10", "low")


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
