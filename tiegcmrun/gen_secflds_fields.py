#!/usr/bin/env python
"""Regenerate secflds_fields.json, the fallback SECFLDS field list, from the model source.

usage: python gen_secflds_fields.py [TIEGCM_SRC_DIR]     (default: $TIEGCMHOME/src, else ../src)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from misc import scan_secflds_fields, SECFLDS_FIELDS_FILE  # noqa: E402


def main(argv):
    if len(argv) > 1:
        src = argv[1]
    elif os.environ.get("TIEGCMHOME"):
        src = os.path.join(os.environ["TIEGCMHOME"], "src")
    else:
        src = os.path.join(HERE, "..", "src")
    names = scan_secflds_fields(src)
    if not names:
        raise SystemExit(f"No Fortran source in {src}")
    with open(SECFLDS_FIELDS_FILE, "w", encoding="utf-8") as f:
        json.dump({"source": "tiegcm/src (fields.F, diags.F, addfld/mkdiag calls, input.F "
                             "secflds_mandatory); regenerate with gen_secflds_fields.py",
                   "fields": sorted(names, key=lambda n: (n.upper(), n))}, f, indent=1)
        f.write("\n")
    print(f"{len(names)} fields -> {SECFLDS_FIELDS_FILE}")


if __name__ == "__main__":
    main(sys.argv)
