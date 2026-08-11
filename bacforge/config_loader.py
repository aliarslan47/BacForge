"""Config yükleme: ${ENV} değişkenlerini çözer. KURAL: kodda mutlak yol YOK."""
from __future__ import annotations

import os
import re
from pathlib import Path

try:
    import yaml
except ImportError:  # pyyaml çekirdek env'de var; yoksa anlamlı hata
    yaml = None

_ENV_PATTERN = re.compile(r"\$\{([A-Z_][A-Z0-9_]*)\}")


def ensure_env_defaults() -> str:
    """BACFORGE_HOME / _DB / _WORK ayarlanmadıysa makul varsayılan ver (taşınabilirlik)."""
    home = os.environ.get("BACFORGE_HOME") or str(Path(__file__).resolve().parents[1])
    os.environ.setdefault("BACFORGE_HOME", home)
    os.environ.setdefault("BACFORGE_DB", str(Path(home) / "databases"))
    os.environ.setdefault("BACFORGE_WORK", str(Path(home) / "runs"))
    return home


def _expand(value):
    if isinstance(value, str):
        prev, cur = None, value
        while cur != prev:  # iç içe ${A}/${B} için tekrar çöz
            prev = cur
            cur = _ENV_PATTERN.sub(lambda m: os.environ.get(m.group(1), m.group(0)), cur)
        return cur
    if isinstance(value, dict):
        return {k: _expand(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_expand(v) for v in value]
    return value


def load_config(path: str | os.PathLike | None = None) -> dict:
    ensure_env_defaults()
    if path is None:
        path = Path(os.environ["BACFORGE_HOME"]) / "config" / "config.yaml"
    path = Path(path)
    if yaml is None:
        raise RuntimeError("pyyaml kurulu değil. 'conda env create -f environment.yml' çalıştır.")
    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}
    return _expand(raw)
