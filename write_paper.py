"""Build a six-figure research Notebook with executable analysis in an appendix."""
from pathlib import Path
import base64
import json
import re
import nbformat as nb
ROOT=Path(__file__).resolve().parent
CODE={
'setup':'''from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from IPython.display import display, HTML
from run_study import manifest, folder
from analyze_study import read_job, reference
jobs = manifest()
checked, block_rows = [], []
for job in jobs:
    row, blocks, _ = read_job(job)
    checked.append(row)
    block_rows.extend(blocks)
    meta = json.loads((folder(job)/"metadata.json").read_text(encoding="utf8"))
    input_text = (folder(job)/"input.in").read_text(encoding="ascii")
    assert hashlib.sha256(input_text.encode("ascii")).hexdigest() == meta["input_sha256"]
assert len(checked) == 24
figures = json.loads(Path("docs/figure_manifest.json").read_text(encoding="utf8"))
for item in figures:
    actual = hashlib.sha256(Path(item["image"]).read_bytes()).hexdigest()
    assert actual == item["image_sha256"]
display(HTML("<style>" + Path("docs/paper_style.css").read_text(encoding="utf8") + "</style>"))
print("Verified 24 completed trajectories and six selected figure assets.")''',
'statistics':'''runs = pd.DataFrame(checked)
production = runs[runs.group == "production"]
a = production.pivot(index="temperature", columns="replicate", values="a_A")
strain = a.div(a.loc[300], axis=1) - 1
coefficients = strain.loc[800] / 500
estimate = coefficients.mean()
tcrit = 4.3026527299
halfwidth = tcrit * coefficients.std(ddof=1) / np.sqrt(3)
_, refalpha = reference()
bias = 100 * (estimate / refalpha - 1)
blocks = pd.DataFrame(block_rows).pivot(index="job", columns="block_ps", values="a_sem_A")
ratio = blocks[50] / blocks[20]
assert np.isclose(estimate*1e6, 17.8882558243, rtol=1e-8)
print(f"MD coefficient: {estimate*1e6:.3f} ± {halfwidth*1e6:.3f} ×10⁻⁶ K⁻¹ (95% statistical CI)")
print(f"NIST coefficient: {refalpha*1e6:.4f} ×10⁻⁶ K⁻¹; relative deviation: {bias:.2f}%")
print(f"SEM(50 ps)/SEM(20 ps): median {ratio.median():.3f}; observed range {ratio.min():.3f}–{ratio.max():.3f}")
display(pd.DataFrame({"Repeat": coefficients.index, "Coefficient / 10⁻⁶ K⁻¹": coefficients.values*1e6}).style.hide(axis="index"))''',
'protocol':'''# Display recorded inputs; no long simulations are launched.
keys = ("lattice ", "region ", "boundary ", "pair_style ", "pair_coeff ", "timestep ", "fix integrator ", "run ")
for filename in ["outputs/model/build.in", "outputs/production/prod_T300_r1/input.in"]:
    print("\\n" + filename)
    print("\\n".join(line for line in Path(filename).read_text(encoding="ascii").splitlines() if line.startswith(keys)))''',
'rebuild':'''import ast, subprocess, sys
REBUILD_FIGURES = False
for item in figures:
    script = item.get("script")
    if script:
        ast.parse(Path(script).read_text(encoding="utf8"))
        print(f"Figure {item['number']}: {script}")
        if REBUILD_FIGURES:
            subprocess.run([sys.executable, script], check=True)
print("Figure source scripts verified. Figure rebuilding is", "enabled." if REBUILD_FIGURES else "disabled.")'''
}


def prose_cells(text):
    cells=[]
    for section in re.split(r'(?=^#{1,3} )',text,flags=re.M):
        if not section.strip():
            continue
        section=section.strip()
        # Native Markdown remains legible in editors that suppress custom CSS.
        cells.append(nb.v4.new_markdown_cell(section))
    return cells


def main():
    source=(ROOT/'docs/manuscript.md').read_text(encoding='utf8')
    figures=json.loads((ROOT/'docs/figure_manifest.json').read_text(encoding='utf8'))
    by_number={str(item['number']):item for item in figures}
    parts=re.split(r'<!-- (figure:\d+|code:\w+) -->',source)
    cells=[]
    for i,part in enumerate(parts):
        if i%2:
            kind,key=part.split(':')
            if kind=='figure':
                item=by_number[key]
                attachment=f'figure{key}.png'
                data=base64.b64encode((ROOT/item['notebook_image']).read_bytes()).decode('ascii')
                cell=nb.v4.new_markdown_cell(
                    f"![Figure {key}. {item['title']}](attachment:{attachment})\n\n{item['caption']}",
                    metadata={'paper_figure':int(key),'figure_source':item['image']})
                cell.attachments={attachment:{'image/png':data}}
                cells.append(cell)
            else:
                compile(CODE[key],'<notebook:'+key+'>','exec')
                cells.append(nb.v4.new_code_cell(CODE[key],metadata={'tags':['analysis-appendix']}))
        else:
            cells.extend(prose_cells(part))
    book=nb.v4.new_notebook(cells=cells,metadata={
        'kernelspec':{'display_name':'Copper MD (Python)','language':'python','name':'copper-md'},
        'language_info':{'name':'python'},'language':'en',
        'manuscript_revision':'2026-10-20 six-figure journal-style integration',
        'figure_count':6})
    nb.validate(book)
    nb.write(book,ROOT/'铜热膨胀研究.ipynb')
    print('Notebook built:',len(cells),'cells, six embedded figures,',len(CODE),'appendix code cells.')


if __name__=='__main__':
    main()
