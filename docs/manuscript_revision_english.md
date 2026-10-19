# English manuscript revision — 19 October 2026

The paper was expanded and rewritten in English. All manuscript prose, references, figure labels, legends, captions, table headings, code comments, and displayed analysis text are in English. The existing Chinese Notebook filename is retained to preserve project links.

## Content and scope

- Approximately 4,015 English words before the References section, including the abstract and captions, by a word-token count.
- Added explanations of the NPT choice, the reference-temperature conversion, coefficient definitions, the physical meaning of volume changes, correlation-aware statistics, and the interpretation limits of endpoint checks.
- No new trajectories, invented data, or new scientific citations were introduced. All 24 completed trajectories retain their existing sampling definitions.
- The key result remains 17.888 ± 0.015 × 10⁻⁶ K⁻¹ (95% statistical interval half-width), with a −3.85% deviation from the consistently defined NIST reference.

## Presentation

- Main title: 32 px; section headings: 25 px; subsection headings: 20 px; body: 17 px with 1.72 line spacing.
- Serif body text; sans-serif tables, captions, and figure labels; selective bold for central findings and definitions.
- HTML reader: constrained column width, distinct caption and reference sizing, and nine expandable analysis-code sections. The Notebook retains all executable code.
- Native Markdown heading levels remain usable when an editor suppresses custom CSS. Exact Notebook font rendering depends on the editor and its trust settings.

## Verification

- Nine executable code cells completed without errors; five PNG figures are embedded.
- All five regenerated English figure images were visually inspected; labels were not clipped.
- No CJK characters were found in the manuscript or displayed code outputs.
- The HTML export was structurally checked, including its headings, figures, tables, and code-fold elements. Automatic browser viewing of the local HTML was blocked because the browser tool disallows file URLs; no workaround was used. Consequently the full HTML page has not received a browser screenshot review.
- The current bibliography metadata audit is `reference_verification_20261019.json`; earlier primary-source and full-text checks remain in `reference_verification.json` and `reference_audit.md`. Failed metadata requests are explicitly recorded rather than reported as successful.
- Original prose, builder, and Notebook backups are in `.tmp/manuscript_before_english_20261019/` (local, ignored by Git).

## Rebuilding the paper

Run `write_paper.py`, then `execute_notebook.py`, then `export_paper_html.py` using the project Python environment. These commands rebuild the report from completed data and do not launch new MD runs. Run `write_readme.py` to refresh the project entry page.
