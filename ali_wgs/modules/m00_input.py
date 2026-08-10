"""M00 -- Input & Automatic Data Detection
Validates files, calculates read stats (count, mean, median, N50, paired R1/R2),
detects data type (SHORT_READ, LONG_READ, HYBRID, ASSEMBLY_INPUT), generates checksums and metadata.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from ..detect import detect_platform
from .base import Module


def compute_sha256(filepath: Path, chunk_size: int = 1048576) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


class InputDetectionModule(Module):
    number = "00"
    name = "input_detection"
    folder = "M00_INPUT_AUTO_DETECTION"
    enabled_key = "input_detection"

    def inputs(self):
        return [self.ctx.input_path]

    def outputs(self):
        return [self.out_dir / "data_type.json"]

    def run(self):
        self.check_inputs()

        p = Path(self.ctx.input_path)
        files = sorted([f for f in p.iterdir() if f.is_file()]) if p.is_dir() else [p]

        # 1. Validation & Checksums
        checksums = {}
        validations = []
        for f in files:
            sha = compute_sha256(f)
            checksums[f.name] = sha
            validations.append({
                "file_name": f.name,
                "file_path": str(f),
                "size_bytes": f.stat().st_size,
                "sha256": sha,
                "exists": f.exists()
            })

        val_file = self.sub_dir("01_input") / "file_validation.json"
        with open(val_file, "w", encoding="utf-8") as fh:
            json.dump(validations, fh, indent=2, ensure_ascii=False)

        chk_file = self.sub_dir("01_input") / "checksums.sha256"
        with open(chk_file, "w", encoding="utf-8") as fh:
            for fname, sha in checksums.items():
                fh.write(f"{sha}  {fname}\n")

        # 2. Platform & Data Type Detection
        detection = detect_platform(self.ctx.input_path, self.ctx.config)
        self.ctx.detection = detection

        # Map to spec data types: SHORT_READ, LONG_READ, HYBRID, ASSEMBLY_INPUT
        read_type = detection.get("read_type", "short")
        if read_type == "assembly_input":
            data_type = "ASSEMBLY_INPUT"
        elif read_type == "long":
            data_type = "LONG_READ"
        elif detection.get("paired"):
            data_type = "SHORT_READ"
        else:
            data_type = "SHORT_READ" if detection.get("mean_read_length", 0) <= 350 else "LONG_READ"

        # Check hybrid input if both short and long fastq files are present
        if p.is_dir():
            has_short = any("_R1" in f.name or "_R2" in f.name for f in files if ".fastq" in f.name or ".fq" in f.name)
            has_long = any("long" in f.name.lower() or "ont" in f.name.lower() or "nanopore" in f.name.lower() for f in files)
            if has_short and has_long:
                data_type = "HYBRID"

        detection["data_type"] = data_type
        self.ctx.detection["data_type"] = data_type

        # 3. Standardized outputs
        std_dir = self.sub_dir("04_standardized")
        with open(std_dir / "data_type.json", "w", encoding="utf-8") as fh:
            json.dump({"data_type": data_type, "detected_platform": detection.get("platform")}, fh, indent=2)

        with open(std_dir / "platform_detection.json", "w", encoding="utf-8") as fh:
            json.dump(detection, fh, indent=2, ensure_ascii=False)

        tsv_path = std_dir / "read_statistics.tsv"
        with open(tsv_path, "w", encoding="utf-8") as fh:
            fh.write("Metric\tValue\n")
            fh.write(f"data_type\t{data_type}\n")
            fh.write(f"platform\t{detection.get('platform', 'unknown')}\n")
            fh.write(f"mean_read_length\t{detection.get('mean_read_length', 0)}\n")
            fh.write(f"read_n50\t{detection.get('n50', 0)}\n")
            fh.write(f"paired_end\t{detection.get('paired', False)}\n")

        # 4. Update manifest
        manifest_path = Path(self.ctx.run_dir) / "project_manifest.json"
        if manifest_path.exists():
            with open(manifest_path, "r", encoding="utf-8") as fh:
                man = json.load(fh)
            man["detected_data_type"] = data_type
            man["detected_platform"] = detection.get("platform")
            man["module_status"]["M00"] = "PASS"
            with open(manifest_path, "w", encoding="utf-8") as fh:
                json.dump(man, fh, indent=2, ensure_ascii=False)

        # 5. Write M00_summary.json
        self.write_summary(
            status="PASS",
            statistics={
                "data_type": data_type,
                "platform": detection.get("platform"),
                "mean_read_length": detection.get("mean_read_length"),
                "read_n50": detection.get("n50"),
                "paired": detection.get("paired", False)
            },
            details={"file_count": len(files)}
        )

        return detection
