# Research-article reading edition

`main.tex` is generated from `docs/manuscript.md` and `docs/figure_manifest.json` by `export_journal_pdf.py`. Edit the source Markdown and figure manifest for persistent changes. The PDF uses a generic materials-research article layout: A4, two columns, Times New Roman, numbered equations and references, full-width figures with captions below, and three-line tables with captions above. No journal logo, publication metadata, or acceptance status is implied.

The companion Notebook is the executable document. Its four analysis cells are collected in Appendix A; the PDF omits executable cells. Figure 1 is a disclosed AI-assisted schematic, whereas Figures 2–6 are quantitative plots generated from the recorded dataset. Names and affiliation remain explicit placeholders pending completion by the authors.

Run `export_journal_pdf.py` using the project environment with XeLaTeX installed. The current script uses the locally available TeX Live executable; update that path for another machine. The saved PDF and HTML can be read without rebuilding them.
