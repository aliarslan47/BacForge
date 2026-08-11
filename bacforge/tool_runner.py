"""ToolRunner: her harici araç çağrısı + provenance (tool/sürüm/süre/exit/DB) kaydı."""
from __future__ import annotations

import json
import shlex
import subprocess
import time
from pathlib import Path


class ToolRunner:
    def __init__(self, logs_dir: str | Path):
        self.logs_dir = Path(logs_dir)
        self.logs_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _as_list(cmd):
        return shlex.split(cmd) if isinstance(cmd, str) else list(cmd)

    def _wrap_env(self, cmd_list, conda_env):
        return (["conda", "run", "-n", conda_env] + cmd_list) if conda_env else cmd_list

    def _get_version(self, version_cmd, conda_env):
        if not version_cmd:
            return None
        try:
            vc = self._wrap_env(self._as_list(version_cmd), conda_env)
            out = subprocess.run(vc, capture_output=True, text=True, timeout=60)
            return (out.stdout or out.stderr).strip().splitlines()[0] if (out.stdout or out.stderr) else None
        except Exception as exc:  # sürüm alınamazsa pipeline durmaz
            return f"bilinmiyor ({exc})"

    def run(self, name, cmd, conda_env=None, version_cmd=None,
            db_version=None, cwd=None, check=True, stdout_path=None) -> dict:
        """stdout_path verilirse araç stdout'u oraya (TSV vb.), stderr log'a yazılır."""
        cmd_list = self._wrap_env(self._as_list(cmd), conda_env)
        log_path = self.logs_dir / f"{name}.log"
        prov_path = self.logs_dir / f"{name}.provenance.json"

        version = self._get_version(version_cmd, conda_env)
        start = time.time()
        with open(log_path, "w") as log:
            log.write(f"# CMD: {' '.join(cmd_list)}\n")
            log.flush()
            if stdout_path:
                with open(stdout_path, "w") as out:
                    proc = subprocess.run(cmd_list, stdout=out, stderr=log,
                                          cwd=cwd, text=True)
            else:
                proc = subprocess.run(cmd_list, stdout=log, stderr=subprocess.STDOUT,
                                      cwd=cwd, text=True)
        duration = round(time.time() - start, 2)

        prov = {
            "tool": name,
            "command": " ".join(cmd_list),
            "version": version,
            "db_version": db_version,
            "start_epoch": round(start, 2),
            "duration_sec": duration,
            "exit_code": proc.returncode,
            "log": str(log_path),
        }
        with open(prov_path, "w") as fh:
            json.dump(prov, fh, indent=2, ensure_ascii=False)

        if check and proc.returncode != 0:
            raise RuntimeError(f"[{name}] başarısız (exit {proc.returncode}). Log: {log_path}")
        return prov
