"""Kaynak tespiti: runtime'da otomatik + AGRESİF. Sabit değer YOK -> taşınabilir."""
from __future__ import annotations

import os


def _total_ram_gb() -> float:
    try:
        import psutil
        return psutil.virtual_memory().total / (1024 ** 3)
    except Exception:
        try:
            pages = os.sysconf("SC_PHYS_PAGES")
            page = os.sysconf("SC_PAGE_SIZE")
            return pages * page / (1024 ** 3)
        except (ValueError, OSError):
            return 8.0  # son çare


def detect_resources(config: dict | None = None) -> dict:
    res = (config or {}).get("resources", {}) if config else {}
    cores = os.cpu_count() or 1
    ram = _total_ram_gb()

    reserve_cores = int(res.get("reserve_cores", 1))
    reserve_mem = float(res.get("reserve_memory_gb", 4))
    thread_frac = float(res.get("thread_fraction", 1.0))
    mem_frac = float(res.get("memory_fraction", 0.90))

    threads = max(1, int(round((cores - reserve_cores) * thread_frac)))
    memory_gb = max(1, int(min(ram * mem_frac, ram - reserve_mem)))

    return {
        "cores_total": cores,
        "ram_gb_total": round(ram, 1),
        "threads": threads,         # araçlara verilecek -j/-t
        "memory_gb": memory_gb,     # araçlara verilecek bellek tavanı
        "profile": res.get("profile", "aggressive"),
    }
