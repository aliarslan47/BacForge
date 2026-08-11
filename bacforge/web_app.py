"""FastAPI Web Application Server for Antigravity Bacterial WGS Platform.
Provides REST API endpoints for runs, sample submission, real-time log streaming, data visualizers, and project exports.
"""
from __future__ import annotations

import asyncio
import json
import os
from pathlib import Path
import shutil

from fastapi import FastAPI, HTTPException, BackgroundTasks, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from .config_loader import load_config
from .detect import detect_platform
from .orchestrator import Orchestrator

app = FastAPI(
    title="Antigravity Bacterial WGS Platform",
    description="Automated Bacterial WGS Bioinformatics Analysis & Visualization Platform",
    version="1.0.0"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

CONFIG = load_config()
WORK_DIR = Path(CONFIG["paths"]["work"])
WORK_DIR.mkdir(parents=True, exist_ok=True)
STATIC_DIR = Path(__file__).parent / "static"
STATIC_DIR.mkdir(parents=True, exist_ok=True)

app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


@app.get("/", response_class=HTMLResponse)
def index():
    index_file = STATIC_DIR / "index.html"
    if index_file.exists():
        return index_file.read_text(encoding="utf-8")
    return "<h1>Antigravity Bacterial WGS Platform</h1>"


@app.get("/api/runs")
def list_runs():
    runs = []
    if WORK_DIR.exists():
        for d in sorted(WORK_DIR.iterdir(), reverse=True):
            if d.is_dir() and (d / "project_manifest.json").exists():
                try:
                    with open(d / "project_manifest.json", "r", encoding="utf-8") as fh:
                        man = json.load(fh)
                    runs.append({
                        "run_id": d.name,
                        "project_id": man.get("project_id", d.name),
                        "started_at": man.get("started_at"),
                        "input_path": man.get("input_path"),
                        "detected_data_type": man.get("detected_data_type"),
                        "detected_platform": man.get("detected_platform"),
                        "status": "COMPLETED" if "M18" in man.get("module_status", {}) else "RUNNING"
                    })
                except Exception:
                    pass
    return {"runs": runs}


@app.get("/api/run/{run_id}")
def get_run_details(run_id: str):
    run_path = WORK_DIR / run_id
    if not run_path.exists():
        raise HTTPException(status_code=404, detail="Run not found")

    manifest_path = run_path / "project_manifest.json"
    manifest = {}
    if manifest_path.exists():
        with open(manifest_path, "r", encoding="utf-8") as fh:
            manifest = json.load(fh)

    dash_file = run_path / "M17_STATISTICS_VISUALIZATION" / "04_standardized" / "dashboard_data.json"
    dashboard_data = {}
    if dash_file.exists():
        with open(dash_file, "r", encoding="utf-8") as fh:
            dashboard_data = json.load(fh)

    return {
        "run_id": run_id,
        "manifest": manifest,
        "dashboard": dashboard_data
    }


@app.post("/api/detect")
def detect_sample(input_path: str = Form(...)):
    p = Path(input_path)
    if not p.exists():
        raise HTTPException(status_code=400, detail=f"Path '{input_path}' does not exist.")
    detection = detect_platform(str(p), CONFIG)
    return detection


def _run_pipeline_background(input_path: str, override_type: str | None = None):
    orch = Orchestrator()
    orch.run(input_path)


@app.post("/api/submit")
def submit_sample(
    background_tasks: BackgroundTasks,
    input_path: str = Form(None),
    file: UploadFile = File(None)
):
    if file:
        upload_dir = WORK_DIR / "uploads"
        upload_dir.mkdir(parents=True, exist_ok=True)
        dest = upload_dir / file.filename
        with open(dest, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        target = str(dest)
    elif input_path:
        target = input_path
    else:
        raise HTTPException(status_code=400, detail="Provide input_path or file upload")

    background_tasks.add_task(_run_pipeline_background, target)
    return {"status": "SUCCESS", "message": f"Run queued for input: {target}"}


@app.get("/api/export/{run_id}/{export_type}")
def export_file(run_id: str, export_type: str):
    run_path = WORK_DIR / run_id
    if not run_path.exists():
        raise HTTPException(status_code=404, detail="Run not found")

    if export_type == "zip":
        zip_file = run_path / "M18_REPORT_EXPORT" / "04_standardized" / "PROJECT_COMPLETE.zip"
        if zip_file.exists():
            return FileResponse(path=zip_file, filename=f"{run_id}_PROJECT_COMPLETE.zip", media_type="application/zip")
    elif export_type == "html":
        html_file = run_path / "M18_REPORT_EXPORT" / "04_standardized" / "report.html"
        if html_file.exists():
            return FileResponse(path=html_file, filename=f"{run_id}_report.html", media_type="text/html")
    elif export_type == "fasta":
        fa_file = run_path / "M04_POLISHING_GENOME_QC" / "04_standardized" / "genome.fasta"
        if fa_file.exists():
            return FileResponse(path=fa_file, filename=f"{run_id}_genome.fasta", media_type="text/plain")

    raise HTTPException(status_code=404, detail="Requested export file not found")
