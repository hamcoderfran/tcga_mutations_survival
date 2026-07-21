#!/usr/bin/env python3
"""CLI wrapper: integrate priority 1–12 datasources into exhalepath_atlas."""
from __future__ import annotations

import argparse
from pathlib import Path

from exhalepath.datasources import integrate_all_datasources


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--priorities", default=None, help="e.g. 1,2,3")
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    pri = None
    if args.priorities:
        pri = [int(x) for x in args.priorities.split(",") if x.strip()]
    man = integrate_all_datasources(root=root, offline=args.offline, priorities=pri)
    print(f"Integrated {man['n_datasources']} datasources")
    print(f"Diseases: {man.get('capability', {}).get('n_diseases')}")
    print(f"Manifest: {root / 'data' / 'datasources' / 'INTEGRATION_MANIFEST.json'}")


if __name__ == "__main__":
    main()
