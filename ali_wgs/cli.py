"""CLI: `ali-wgs info | detect | run`. Web katmanı (sonra) bu çekirdeği dışarıdan çağırır."""
from __future__ import annotations

import argparse
import json
import sys

from .config_loader import load_config, ensure_env_defaults
from .resources import detect_resources
from .detect import detect_platform
from .orchestrator import Orchestrator


def _cmd_info(args):
    ensure_env_defaults()
    cfg = {}
    try:
        cfg = load_config(args.config)
    except Exception as exc:
        print(f"(config okunamadı, varsayılan kaynak politikası: {exc})", file=sys.stderr)
    res = detect_resources(cfg)
    print(json.dumps({"resources": res, "paths": cfg.get("paths", {})},
                     indent=2, ensure_ascii=False))


def _cmd_detect(args):
    cfg = {}
    try:
        cfg = load_config(args.config)
    except Exception:
        pass
    print(json.dumps(detect_platform(args.input, cfg), indent=2, ensure_ascii=False))


def _cmd_run(args):
    orch = Orchestrator(args.config)
    orch.run(args.input, config_path=args.config, resume=not args.no_resume)


def _cmd_server(args):
    import uvicorn
    print(f"Starting Antigravity Bacterial WGS Web Server on http://{args.host}:{args.port}")
    uvicorn.run("ali_wgs.web_app:app", host=args.host, port=args.port, reload=args.reload)


def build_parser():
    p = argparse.ArgumentParser(prog="ali-wgs", description="Antigravity Bacterial WGS Bioinformatics Platform")
    p.add_argument("--config", default=None, help="config.yaml yolu (vars: ALI_WGS_HOME/config)")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("info", help="Tespit edilen kaynakları/yolları göster").set_defaults(func=_cmd_info)

    d = sub.add_parser("detect", help="Girdiden platform/okuma-tipi tespit et (env'siz)")
    d.add_argument("--input", required=True, help="okuma dosyası veya dizini")
    d.set_defaults(func=_cmd_detect)

    r = sub.add_parser("run", help="Pipeline'ı uçtan uca çalıştır")
    r.add_argument("--input", required=True, help="okuma dosyası veya dizini")
    r.add_argument("--no-resume", action="store_true", help="resume'u kapat, baştan çalıştır")
    r.set_defaults(func=_cmd_run)

    s = sub.add_parser("server", help="Web Dashboard sunucusunu başlat (FastAPI)")
    s.add_argument("--host", default="0.0.0.0", help="Host IP")
    s.add_argument("--port", type=int, default=8000, help="Port (vars: 8000)")
    s.add_argument("--reload", action="store_true", help="Auto reload")
    s.set_defaults(func=_cmd_server)

    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
