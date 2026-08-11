"""Modül 01 — Input Detection.
INPUT : ham okuma (FASTQ/FAST5/POD5/BAM) — ctx.input_path
OUTPUT: 01_Input/platform.json (platform, read_type, pair, N50, kalite)
"""
from __future__ import annotations

import json

from ..detect import detect_platform
from .base import Module


class InputDetectionModule(Module):
    number = "01"
    name = "input_detection"
    folder = "01_Input"
    enabled_key = "input_detection"

    def inputs(self):
        return [self.ctx.input_path]

    def outputs(self):
        return [self.out_dir / "platform.json"]

    def run(self):
        self.check_inputs()
        detection = detect_platform(self.ctx.input_path, self.ctx.config)
        self.ctx.detection = detection  # sonraki modüller okur
        with open(self.outputs()[0], "w") as fh:
            json.dump(detection, fh, indent=2, ensure_ascii=False)
        return detection
