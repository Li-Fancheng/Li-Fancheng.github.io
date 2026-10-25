"""Package the actual source, data and checks, using stable ZIP metadata."""

from pathlib import Path
import hashlib
import json
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parent
FILES = [
    "README.md",
    "requirements.txt",
    "unified_vertex.py",
    "plot_unified_case.py",
    "verify_unified_case.py",
    "package_case.py",
    "unified_vertex_data.json",
    "unified_vertex_matrices.npz",
    "verification.json",
]
manifest = {
    name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in FILES
}
with ZipFile(ROOT / "reproducibility.zip", "w") as out:
    for name in FILES + ["manifest.json"]:
        info = ZipInfo(name, date_time=(2026, 10, 25, 0, 0, 0))
        info.compress_type = ZIP_DEFLATED
        info.external_attr = 0o644 << 16
        raw = (
            (json.dumps(manifest, indent=2) + "\n").encode()
            if name == "manifest.json"
            else (ROOT / name).read_bytes()
        )
        out.writestr(info, raw)
print("Packaged source, numerical outputs and verification.")
