"""Tüm modüllerin uyduğu tekdüze spesifikasyon sözleşmesi.
Mxx_MODULE/
├── 01_input/
├── 02_work/
├── 03_native_outputs/
├── 04_standardized/
├── 05_statistics/
├── 06_visualization/
├── 07_logs/
├── 08_metadata/
└── Mxx_summary.json
"""
from __future__ import annotations

from abc import ABC, abstractmethod
import json
from pathlib import Path
import time


class Module(ABC):
    number: str = "00"
    name: str = "base"
    folder: str = "M00_INPUT_AUTO_DETECTION"
    enabled_key: str | None = None  # config.modules altındaki anahtar

    SUBDIRS = []

    def __init__(self, ctx):
        self.ctx = ctx  # RunContext (orchestrator.py)

    @property
    def out_dir(self) -> Path:
        d = Path(self.ctx.run_dir) / self.folder
        d.mkdir(parents=True, exist_ok=True)
        return d

    def sub_dir(self, name: str) -> Path:
        # User requested flat structure. Return out_dir directly.
        return self.out_dir

    def init_subdirs(self):
        pass

    def write_summary(
        self,
        status: str,
        statistics: dict | None = None,
        warnings: list[str] | None = None,
        errors: list[str] | None = None,
        details: dict | None = None,
    ) -> Path:
        """Standardized summary output generator (Mxx_summary.json).
        Status must be one of: PASS, WARNING, FAIL, NOT_APPLICABLE, SKIPPED.
        """
        valid_statuses = {"PASS", "WARNING", "FAIL", "NOT_APPLICABLE", "SKIPPED"}
        if status not in valid_statuses:
            status = "PASS"

        summary = {
            "module_number": self.number,
            "module_name": self.name,
            "folder": self.folder,
            "status": status,
            "timestamp": time.ctime(),
            "statistics": statistics or {},
            "warnings": warnings or [],
            "errors": errors or [],
            "details": details or {},
        }

        sum_path = self.out_dir / f"M{self.number}_summary.json"
        with open(sum_path, "w", encoding="utf-8") as fh:
            json.dump(summary, fh, indent=2, ensure_ascii=False)
        return sum_path

    # --- INPUT/OUTPUT sözleşmesi ---
    @abstractmethod
    def inputs(self) -> list:
        """Bu modülün ihtiyaç duyduğu dosyalar (önceki modülün çıktıları)."""

    @abstractmethod
    def outputs(self) -> list:
        """Bu modülün üreteceği dosyalar (resume kontrolünün dayanağı)."""

    @abstractmethod
    def run(self):
        """Aracı ToolRunner ile çağır, çıktıyı out_dir'e yaz."""

    # --- Ortak davranış ---
    def is_done(self) -> bool:
        outs = self.outputs()
        return bool(outs) and all(Path(p).exists() for p in outs)

    def check_inputs(self):
        missing = [str(p) for p in self.inputs() if not Path(p).exists()]
        if missing:
            raise FileNotFoundError(
                f"[{self.number}_{self.name}] eksik girdi: {missing}. "
                f"Önceki modül çalışmamış olabilir."
            )

    def validate(self) -> bool:
        return self.is_done()

