# Six-figure manuscript integration and layout review

Date: 20 October 2026.

## Delivered editions

- `铜热膨胀研究.ipynb`: 37 cells, six self-contained figure attachments, two article tables, and four executed appendix code cells. All code cells completed without errors. The selected plot assets are validated against their SHA256 hashes.
- `results/copper_thermal_expansion_paper.html`: single-column research-article reader with six figure/caption objects, three-line tables, numerical references, figure links, and expandable appendix code.
- `results/copper_thermal_expansion_journal.pdf`: 10-page A4 two-column reading edition, with Times New Roman, 9.5 pt body text, 8.5 pt captions, numbered equations, page numbers, full-width figures, and table captions above three-line tables.
- `paper_layout/main.tex`: generated, editable typesetting source; `export_journal_pdf.py` reproduces the PDF using the locally installed XeLaTeX executable.

This is a generic materials-research article layout. SCI is an indexing system, not a single typesetting template. No journal logo, acceptance statement, DOI, publication date, or other fabricated publication metadata has been added.

## Figure integration

1. Crystal model and workflow: inserted in Methods, explicitly described as schematic. The illustration was generated with AI assistance; its atom-count annotation and protocol settings were checked, and its use is disclosed in Data and Code Availability. It is not a trajectory snapshot or a source of numerical measurements.
2. Expansion, fluctuation distributions, and repeat-resolved volumes: replaces the original two-panel mean plot. The text now explains within-trajectory centering, correlated distributions, and repeat-specific volume normalization.
3. NIST comparison, residuals, and coefficients: replaces the original comparison graphic, retaining the 300 K normalization and MD-only statistical interval definition.
4. Settings matrix and coefficient changes: replaces the old bar chart, using r1 consistently as the baseline.
5. Three-dimensional trajectories and endpoint diagnostics: replaces the former supplementary time-series figure. The common volume denominator and display-only rolling means are identified.
6. Block-duration heatmap, endpoint curves, and paired standard errors: replaces the former supplementary block plot. The text distinguishes observed ratios from confidence intervals and avoids claiming complete independence.

Figures are numbered 1–6 in order of first mention; former S1/S2 references have been removed. The selected assets and authoritative captions are recorded in `docs/figure_manifest.json`. Old figures remain available as historical files but are not embedded in the current article.

## Text and references

The existing results were preserved. Additional text explains the newly displayed quantities and their interpretive limits. References 6 and 7 were reordered to match their first appearance (NIST potential entry, then LAMMPS article); the sources themselves are unchanged. Crystallographic directions [100], [010], and [001] remain directions and are not converted into citations.

Coefficient ratios, lattice parameters, and the reference volume in captions use explicit mathematical notation. The MD coefficient remains 17.888 ± 0.015 × 10⁻⁶ K⁻¹ (95% statistical confidence-interval half-width), with a −3.85% deviation from the consistently defined NIST reference. No simulation data or reported physical result was changed during layout work.

## Verification

- The Notebook executed four appendix cells without errors and verified all 24 completed trajectories.
- All six selected native PNG hashes remain unchanged.
- HTML checks confirmed six paired figures/captions, two article tables plus one appendix output table, four expandable code sections, and valid figure anchors.
- The PDF compiled twice with no overfull-box, missing-character, or unresolved-reference warnings.
- The final 10 PDF pages were rendered. A contact sheet and full-size page inspections covered the title/abstract, both text/equation pages, all six figures, both tables, and the bibliography. No observed clipped or overlapping labels, split figure captions, or broken mathematical notation remained.
- The PDF uses vector sources for Figures 2–6; Figure 1 retains the approved raster schematic. The Notebook uses embedded reading-size previews while preserving native PNG and vector PDF files in the repository.
- Local HTML browser screenshot review was not performed because the browser tool blocks local file URLs; structural checks and the visually reviewed PDF provide the recorded validation scope. Notebook CSS support depends on the editor and trust settings; native heading levels remain readable.
- The machine-readable summary is `results/paper_layout_check.json`.

## Remaining author information

Author names and affiliation remain explicit placeholders. The group must supply factual names, affiliation, and contributions before submission. This does not prevent reading or reproducing the current results.

The pre-integration Notebook, prose, and builders are backed up locally in `.tmp/before_six_figure_layout_20261020/`.
