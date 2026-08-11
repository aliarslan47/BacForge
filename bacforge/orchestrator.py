"""Orchestrator: run dizini, modül sırası, resume, kaynak bağlama. Çekirdeğin kalbi."""
from __future__ import annotations

import json
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

from .config_loader import load_config
from .resources import detect_resources
from .tool_runner import ToolRunner
from .modules import REGISTRY

RUN_SUBDIRS = [
    "M00_INPUT_AUTO_DETECTION",
    "M01_READ_QC_PREPROCESSING",
    "M02_TAXONOMIC_QC",
    "M03_GENOME_ASSEMBLY",
    "M04_POLISHING_GENOME_QC",
    "M05_SPECIES_REFERENCE_IDENTIFICATION",
    "M06_GENOME_ANNOTATION",
    "M07_STRAIN_TYPING",
    "M08_AMR",
    "M09_VIRULENCE",
    "M10_PLASMID",
    "M11_MOBILE_GENETIC_ELEMENTS",
    "M12_PHAGE_CRISPR_DEFENSE",
    "M13_VARIANTS_MUTATIONS",
    "M14_GENOMIC_CONTEXT",
    "M15_COMPARATIVE_GENOMICS",
    "M16_PHYLOGENOMICS",
    "M17_STATISTICS_VISUALIZATION",
    "M18_REPORT_EXPORT",
]


@dataclass
class RunContext:
    config: dict
    resources: dict
    run_dir: Path
    runner: ToolRunner
    input_path: str
    detection: dict = field(default_factory=dict)
    manifest: dict = field(default_factory=dict)


class Orchestrator:
    def __init__(self, config_path=None):
        self.config = load_config(config_path)
        self.resources = detect_resources(self.config)

    def _prepare_run_dir(self, input_path: str, config_path) -> Path:
        work = Path(self.config["paths"]["work"])
        inp = Path(input_path)
        label = inp.stem if inp.is_file() else inp.name
        if label in ("raw", "reads", "fastq"):
            label = inp.parent.name
        run_id = time.strftime("%Y%m%d_%H%M%S") + "_" + label
        run_dir = work / run_id
        for sub in RUN_SUBDIRS:
            (run_dir / sub).mkdir(parents=True, exist_ok=True)

        manifest = {
            "project_id": run_id,
            "project_version": "1.0",
            "spec_version": "bacterial_wgs_antigravity_spec_v1.0",
            "started_at": time.ctime(),
            "input_path": str(input_path),
            "detected_data_type": None,
            "detected_platform": None,
            "module_status": {},
            "tool_versions": {},
            "database_versions": {},
            "reference_accessions": [],
            "errors": [],
            "warnings": [],
        }

        with open(run_dir / "project_manifest.json", "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2, ensure_ascii=False)

        readme_text = f"""# Bacterial WGS Analysis Platform Project

Project ID: {run_id}
Created: {time.ctime()}
Input: {input_path}

This directory contains automated Bacterial WGS bioinformatics analysis outputs structured according to Antigravity Technical Specification v1.0.

Modules M00 to M18 contain standardized JSON/TSV data layer, native tool outputs, statistical reports, and visualizations.
"""
        with open(run_dir / "README.txt", "w", encoding="utf-8") as fh:
            fh.write(readme_text)

        return run_dir

    def _enabled(self, module_cls) -> bool:
        key = module_cls.enabled_key
        if key is None:
            return True
        return bool(self.config.get("modules", {}).get(key, True))

    def run(self, input_path: str, config_path=None, resume: bool = True) -> RunContext:
        run_dir = self._prepare_run_dir(input_path, config_path)
        log_dir = run_dir / "M00_INPUT_AUTO_DETECTION"
        runner = ToolRunner(log_dir)
        ctx = RunContext(self.config, self.resources, run_dir, runner, input_path)

        print(f"== Antigravity Bacterial WGS Platform == run: {run_dir.name}")
        print(f"   Kaynak (agresif): {self.resources['threads']} thread / "
              f"{self.resources['memory_gb']} GB (toplam {self.resources['cores_total']} çekirdek)")

        for module_cls in REGISTRY:
            mod = module_cls(ctx)
            # Update tool runner logs dir to point to current module's out_dir
            ctx.runner.logs_dir = mod.out_dir
            tag = f"M{mod.number}_{mod.name}"

            if not self._enabled(module_cls):
                print(f"   [{tag}] config ile kapalı, atlandı")
                mod.write_summary(status="SKIPPED", details={"reason": "Disabled in configuration"})
                continue
            if resume and mod.is_done():
                print(f"   [{tag}] çıktı mevcut, resume ile atlandı")
                continue
            print(f"   [{tag}] çalışıyor...")
            try:
                mod.run()
                if not mod.validate():
                    mod.write_summary(status="FAIL", errors=[f"[{tag}] validation failed"])
                    raise RuntimeError(f"[{tag}] doğrulama başarısız")
                print(f"   [{tag}] ✓ -> {mod.folder}/")
            except Exception as exc:
                print(f"   [{tag}] ✗ HATA: {exc}")
                mod.write_summary(status="FAIL", errors=[str(exc)])
                # Update manifest
                manifest_path = run_dir / "project_manifest.json"
                if manifest_path.exists():
                    with open(manifest_path, "r", encoding="utf-8") as fh:
                        man = json.load(fh)
                    man["module_status"][f"M{mod.number}"] = "FAIL"
                    man["errors"].append(f"M{mod.number}: {exc}")
                    with open(manifest_path, "w", encoding="utf-8") as fh:
                        json.dump(man, fh, indent=2, ensure_ascii=False)
                if mod.number in ("00", "03", "04", "06"):
                    print(f"   [{tag}] kritik modül hatası, işlem durduruluyor.")
                    raise exc
                else:
                    print(f"   [{tag}] kritik olmayan modül hatası, atlanıyor...")

        print(f"== Analiz Başarıyla Tamamlandı. Çıktı: {run_dir}")
        return ctx
