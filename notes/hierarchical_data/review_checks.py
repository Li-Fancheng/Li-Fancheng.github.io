#!/usr/bin/env python3
"""Validate this review's sources and a built Hugo site, using only Python stdlib.

Run from the checkout:
  python3 notes/hierarchical_data/review_checks.py --site /path/to/hugo/output
The baseline records the RHM pages before this expansion. No files are modified.
"""
from __future__ import annotations
import argparse
from collections import Counter
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote, urljoin, urlsplit

ROOT = Path(__file__).resolve().parents[2]
NOTES = ROOT / 'notes/hierarchical_data'
WRAPPERS = ROOT / 'content/zh/notes/learning_and_generalization/hierarchical_data'
OLD = ['generative_learning', 'language_statistics', 'inference_generation',
       'parsing_latents', 'shared_representations']
REQUIRED = {'saxe14', 'saxe19', 'abbe21', 'abbe22', 'abbe23', 'lofi', 'backlofi',
            'delrocq', 'bfc_full', 'dandi25', 'race', 'nichani', 'wang',
            'tabanelli', 'chen', 'daniely'}

class Page(HTMLParser):
    def __init__(self, path):
        super().__init__()
        self.ids, self.links, self.images, self.tags = [], [], [], []
        self.feed(path.read_text())
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        self.tags.append(tag)
        if 'id' in a: self.ids.append(a['id'])
        if tag == 'a' and a.get('href'): self.links.append(a['href'])
        if tag == 'img' and a.get('src'): self.images.append(a['src'])

def run(site, baseline):
    errors, result = [], {'baseline': baseline, 'pages': []}
    series = json.loads((NOTES/'series.json').read_text())
    slugs = [''] + [row['slug'] for row in series]
    cache = {}
    def page_for(path):
        file = site/path.lstrip('/')
        if file.is_dir(): file /= 'index.html'
        if file not in cache:
            cache[file] = Page(file) if file.is_file() and file.suffix == '.html' else None
        return file, cache[file]
    for slug in slugs:
        path = NOTES/slug/'body.md'
        text = path.read_text()
        name = slug or 'directory'
        wrapper = WRAPPERS/slug/('index.md' if slug else '_index.md')
        canonical = '/notes/hierarchical_data/' + (slug+'/' if slug else '')
        if f'url: "{canonical}"' not in wrapper.read_text(): errors.append(f'{name}: canonical URL')
        headings = re.findall(r'^(#{2,3}) (.+)$', text, re.M)
        for level, title in headings:
            if not re.match(r'\d+\.(?:\d+ )?', title): errors.append(f'{name}: unnumbered {title}')
        if len(re.findall(r'^\$\$$', text, re.M)) % 2: errors.append(f'{name}: display delimiters')
        for display, inline in re.findall(r'\$\$([\s\S]*?)\$\$|(?<!\\)\$([^$\n]+?)(?<!\\)\$',text):
            if re.search(r'<[A-Za-z]',display or inline): errors.append(f'{name}: unsafe literal less-than in math')
        if len(re.findall(r'(?<!\\)\$', text)) % 2: errors.append(f'{name}: inline delimiters')
        source_ids = re.findall(r'\{#([^}]+)\}|<span id="([^"]+)"', text)
        source_ids = [a or b for a,b in source_ids]
        if len(source_ids) != len(set(source_ids)): errors.append(f'{name}: duplicate source IDs')
        if slug in OLD or not slug:
            previous = f'notes/hierarchical_data/{slug+"/" if slug else ""}body.md'
            old = subprocess.check_output(['git','show',f'{baseline}:{previous}'],cwd=ROOT,text=True)
            missing = set(re.findall(r'\{#([^}]+)\}',old))-set(source_ids)
            if missing: errors.append(f'{name}: removed anchors {sorted(missing)}')
            if slug and len(re.findall('[\u4e00-\u9fff]',text)) < .99*len(re.findall('[\u4e00-\u9fff]',old)):
                errors.append(f'{name}: substantial RHM text loss')
        file, html = page_for(canonical)
        if html is None:
            errors.append(f'{name}: no HTML at {file}');continue
        duplicate = [k for k,v in Counter(html.ids).items() if v>1 and k != 'TableOfContents']
        # The site's template intentionally renders TOC twice (desktop/mobile).
        duplicate = [x for x in duplicate if not x.startswith('TableOfContents')]
        if duplicate: errors.append(f'{name}: duplicate rendered IDs {duplicate}')
        if any(re.search(r'[{}\\]',tag) for tag in html.tags): errors.append(f'{name}: math parsed as an HTML tag')
        if not set(source_ids) <= set(html.ids): errors.append(f'{name}: missing rendered source IDs')
        for href in html.links + html.images:
            url = urlsplit(urljoin('https://li-fancheng.github.io'+canonical, href))
            if url.netloc != 'li-fancheng.github.io': continue
            dest, target = page_for(unquote(url.path))
            if not dest.exists(): errors.append(f'{name}: missing {href}');continue
            if url.fragment and target and unquote(url.fragment) not in target.ids:
                errors.append(f'{name}: bad fragment {href}')
        result['pages'].append(dict(slug=name, headings=len(headings),
            display_equations=len(re.findall(r'^\$\$$',text,re.M))//2,
            chinese_characters=len(re.findall('[\u4e00-\u9fff]',text)),
            images=len(html.images), links=len(html.links)))
    refs=json.loads((NOTES/'sources.json').read_text())['references']
    required={r.get('id'):r for r in refs if r.get('id') in REQUIRED}
    if set(required)!=REQUIRED: errors.append('Missing required papers')
    for key, ref in required.items():
        for coverage in ref['coverage']:
            slug,_,fragment=coverage.partition('/#')
            _,html=page_for('/notes/hierarchical_data/'+slug+'/')
            if not html or fragment not in html.ids: errors.append(f'{key}: bad coverage {coverage}')
    routes=json.loads((ROOT/'data/hierarchical_data_routes.json').read_text())
    for anchor,path in routes.items():
        _,html=page_for(path)
        if not html or anchor not in html.ids: errors.append(f'legacy {anchor}: missing target')
    # Independent numerical sanity checks for the explanatory analytic examples.
    residual=0.
    for s in (.3,1.,2.):
        m0=s*1e-3
        def m(t): return s/(1+(s/m0-1)*math.exp(-2*s*t))
        for t in (.1,1.,3.):
            dt=1e-5
            residual=max(residual,abs((m(t+dt)-m(t-dt))/(2*dt)-2*m(t)*(s-m(t))))
    if residual>1e-7: errors.append('logistic ODE residual')
    c0,c1,c2=1.,.7,.2
    Q=[[c0,c1,c2,c2],[c1,c0,c2,c2],[c2,c2,c0,c1],[c2,c2,c1,c0]]
    for v,lam in [([1,1,1,1],c0+c1+2*c2),([1,1,-1,-1],c0+c1-2*c2),([1,-1,0,0],c0-c1),([0,0,1,-1],c0-c1)]:
        if max(abs(sum(Q[i][j]*v[j] for j in range(4))-lam*v[i]) for i in range(4))>1e-12:
            errors.append('block hierarchy eigenvector')
    result.update(required_papers=len(required),legacy_routes=len(routes),
                  logistic_max_residual=residual,errors=sorted(set(errors)))
    print(json.dumps(result,ensure_ascii=False,indent=2))
    return bool(errors)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site',type=Path,required=True)
    parser.add_argument('--baseline',default='5fe1208')
    args=parser.parse_args()
    raise SystemExit(run(args.site.resolve(),args.baseline))
