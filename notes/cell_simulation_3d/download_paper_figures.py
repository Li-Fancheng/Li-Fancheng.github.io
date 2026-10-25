"""Verify local paper figures; fetch missing checksum-matching originals.

--check-only uses only the standard library and never accesses the network.
Restoring the Okuda figure from its published PDF also requires pypdf 6.10.0.
"""

import argparse
import hashlib
import io
import json
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    sources = json.loads((ROOT / 'paper_figure_sources.json').read_text())
    for source in sources:
        path = ROOT / source['file']
        if path.exists():
            data = path.read_bytes()
        else:
            if args.check_only:
                raise FileNotFoundError(path)
            request = Request(source['image_url'], headers={'User-Agent': 'CellSimulation3DFigures/1.0'})
            with urlopen(request, timeout=60) as response:
                data = response.read()
            if source.get('source_format') == 'pdf':
                if hashlib.sha256(data).hexdigest() != source['source_sha256']:
                    raise ValueError('Published source PDF changed; inspect before updating.')
                from pypdf import PdfReader
                page = PdfReader(io.BytesIO(data)).pages[source['pdf_page'] - 1]
                data = next(image.data for image in page.images if image.name == source['pdf_image_name'])
        if hashlib.sha256(data).hexdigest() != source['sha256']:
            raise ValueError(f"Figure differs from recorded source; inspect before updating: {source['file']}")
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
    print(f'Verified {len(sources)} original paper figures.')


if __name__ == '__main__':
    main()
