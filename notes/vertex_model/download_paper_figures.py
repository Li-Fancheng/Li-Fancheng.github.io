"""Check original paper figures and fetch missing, checksum-matching files."""

import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    sources = json.loads((ROOT / "paper_figure_sources.json").read_text())
    for source in sources:
        path = ROOT / source["file"]
        if path.exists():
            if digest(path.read_bytes()) != source["sha256"]:
                raise ValueError(f"Local original differs from manifest: {path}")
            continue
        request = Request(source["image_url"], headers={"User-Agent": "VertexModelLectureFigures/1.0"})
        with urlopen(request, timeout=60) as response:
            data = response.read()
        if digest(data) != source["sha256"]:
            raise ValueError(f"Remote image changed; inspect before updating: {source['image_url']}")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    print(f"Verified {len(sources)} unmodified original paper figures.")


if __name__ == "__main__":
    main()
