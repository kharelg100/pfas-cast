# PFAS Censoring-Aware Screening Tool (PFAS-CAST)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23577354.svg)](https://doi.org/10.5281/zenodo.23577354)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Paper](https://img.shields.io/badge/Paper-Cell%20Reports%20Sustainability-8a2be2.svg)](https://doi.org/10.1016/j.crsus.2026.100893)

**Live tool:** https://kharelg100.github.io/pfas-cast/

PFAS-CAST is a browser-based and Python companion tool for exploratory PFAS benchmark screening. It applies the record-level censored-likelihood multiple-imputation (MI) engine from the paper below to user-supplied PFAS monitoring data and configurable benchmark thresholds, so that results below reporting limits do not bias the screening picture.

> Kharel, G., Kirilenko, A. P., Harvey, O. R. and Zou, C. B. (2026). Censoring-aware models quantify systematic underestimation of PFAS occurrence in US drinking-water monitoring. *Cell Reports Sustainability*. https://doi.org/10.1016/j.crsus.2026.100893

Please cite the paper when using the tool.

## Scope

PFAS-CAST is a screening and prioritization companion, not a regulatory-compliance tool. It does not reproduce the paper's national population-weighted occurrence estimand, does not estimate physical PFAS mass, delivered load, intake dose or toxicity-weighted risk, and does not include the paper's optional Gaussian-copula diagnostic.

Default benchmark values are configurable and reflect the April 2024 PFAS NPDWR framework. EPA announced in May 2025 that it would keep the PFOA and PFOS standards while reconsidering PFHxS, PFNA, HFPO-DA and the Hazard Index mixture. Verify current regulatory status, state requirements and program-specific decision rules before using any benchmark for decision support.

## Two interfaces

| Interface | File | Best for |
|---|---|---|
| Browser tool | `index.html` (identical to `PFAS_CAST_v1_1.html`) | Quick screening of state or system subsets; all processing stays in the browser, nothing is uploaded |
| Python companion | `pfas_cast_screening_tool_v1_1.py` | Large files (for example the full UCMR 5 national file), Excel workbook export, publication figures |

## Quick start

### Browser tool

1. Open https://kharelg100.github.io/pfas-cast/ (or open `index.html` locally).
2. Load a comma- or tab-delimited file. EPA UCMR 5 Data Finder exports are recognized without renaming; `tests/pfas_cast_demo.csv` and `tests/ucmr5_style.txt` are included as examples.
3. Review or adjust the benchmark thresholds, run the screening, and export the results workbook.

### Python companion

```bash
pip install -r pfas_cast_requirements_v1_1.txt

python pfas_cast_screening_tool_v1_1.py \
  --input monitoring_data.csv \
  --unit ug/L \
  --output pfas_cast_results_v1_1.xlsx \
  --m 500 \
  --seed 42
```

Demo run with the bundled data:

```bash
python pfas_cast_screening_tool_v1_1.py --demo --output demo_results.xlsx --m 100 --seed 42
```

## Input data

Column names are auto-detected. The conceptual fields are `system_id`, `location_id` or sampling point (optional), `analyte`, `result_value`, `detect_flag` and `reporting_limit`. Units may be µg/L or ng/L (ppt); all internal calculations use µg/L.

Results with `AnalyticalResultsSign` equal to `<` are treated as below the reporting limit; results with `=` are treated as detections. For the full UCMR 5 national file (about 1.9 million rows) use the Python companion and expect run times of minutes rather than seconds.

## Repository structure

```
index.html                              Browser tool served by GitHub Pages
PFAS_CAST_v1_1.html                     Browser tool, archived release copy
pfas_cast_screening_tool_v1_1.py        Python companion (CLI, workbook and figure export)
pfas_cast_screening_engine_v1_1.py      Shared censored-likelihood MI screening engine
pfas_cast_requirements_v1_1.txt         Python dependencies
pfas_cast_method_references_v1_1.json   Method and regulatory references
pfas_cast_manifest_v1_1.json            SHA-256 checksums of the release files
PFAS_CAST_v1_1_README.md                Release notes as archived on Zenodo
PFAS_CAST_v1_1_QA.md                    Quality-assurance note
tests/                                  Demo input files
```

## Citation

```bibtex
@software{kharel_pfas_cast_2026,
  author    = {Kharel, Gehendra},
  title     = {{PFAS Censoring-Aware Screening Tool (PFAS-CAST)}},
  version   = {1.1.1},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23577354},
  url       = {https://doi.org/10.5281/zenodo.23577354}
}
```

## License

MIT License. Copyright (c) 2025-2026 Gehendra Kharel and contributors.

## Author

Gehendra Kharel, Ph.D., STREAM Lab, Department of Environmental and Geological Sciences, Texas Christian University. https://kharelg100.github.io/
