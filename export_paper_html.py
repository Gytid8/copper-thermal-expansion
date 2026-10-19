"""Export the executed Notebook with readable typography and optional code."""
from pathlib import Path
import re
import nbformat
from nbconvert import HTMLExporter
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parent

def main():
    notebook = nbformat.read(ROOT/'铜热膨胀研究.ipynb', as_version=4)
    exporter = HTMLExporter(template_name='lab')
    exporter.exclude_input_prompt = True
    exporter.exclude_output_prompt = True
    html, _ = exporter.from_notebook_node(notebook)
    soup = BeautifulSoup(html, 'html.parser')
    soup.body['class'] = soup.body.get('class', []) + ['paper-reader']
    style = soup.new_tag('style')
    style.string = (ROOT/'docs/paper_style.css').read_text(encoding='utf8')
    soup.head.append(style)
    if soup.title:
        soup.title.string = 'Thermal Expansion of Single-Crystal Copper | Reproducible Research'
    for i, wrapper in enumerate(soup.select('.jp-CodeCell .jp-Cell-inputWrapper'), 1):
        details = soup.new_tag('details', attrs={'class': 'paper-code'})
        wrapper.wrap(details)
        summary = soup.new_tag('summary')
        summary.string = f'Analysis code {i:02d} · click to expand'
        details.insert(0, summary)
    for p in soup.select('.jp-RenderedHTMLCommon p'):
        if re.match(r'^(Figure\s+[S\d]|Table\s+\d)', p.get_text().strip()):
            p['class'] = p.get('class', []) + ['paper-caption']
        elif re.match(r'^\[\d+\]', p.get_text().strip()):
            p['class'] = p.get('class', []) + ['paper-reference']
    out = ROOT/'results/copper_thermal_expansion_paper.html'
    out.write_text(str(soup), encoding='utf8')
    print('Reader exported:', out)

if __name__ == '__main__':
    main()
