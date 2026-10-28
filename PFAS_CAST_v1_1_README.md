# PFAS Censoring-Aware Screening Tool (PFAS-CAST) v1.1.1

PFAS-CAST is a browser-based and Python companion tool for exploratory PFAS benchmark screening with the primary censored-likelihood multiple-imputation engine from:

> Kharel, G., Kirilenko, A. P., Harvey, O. R. and Zou, C. B. Censoring-aware models quantify systematic underestimation of PFAS occurrence in US drinking-water monitoring. *Cell Reports Sustainability* (2026). https://doi.org/10.1016/j.crsus.2026.100893

Please cite the paper when using the tool.

## Scope

PFAS-CAST is a **screening and prioritization** companion, not a regulatory-compliance tool. It applies the manuscript's primary record-level censored-likelihood MI engine to user-supplied PFAS monitoring data and configurable benchmark thresholds. It does **not** reproduce the paper's national monitored-CWS population-weighted occurrence estimand, does **not** estimate physical PFAS mass, delivered load, intake dose or toxicity-weighted risk, and does **not** include the optional Gaussian-copula diagnostic because that diagnostic is secondary and non-load-bearing in the manuscript.

Default benchmark values are configurable and reflect the April 2024 PFAS NPDWR framework. EPA announced in May 2025 that it would keep the PFOA and PFOS standards while announcing intent to reconsider PFHxS, PFNA, HFPO-DA and the Hazard Index mixture. Users should verify current regulatory status, state requirements and program-specific decision rules before using any benchmark for decision support.

## Files

- `PFAS_CAST_v1_1.html`: standalone client-side browser tool. No external JavaScript libraries are loaded; uploaded data are processed locally in the browser session.
- `pfas_cast_screening_tool_v1_1.py`: Python companion that exports an Excel workbook, CSV outputs and reviewer-facing figures.
- `pfas_cast_screening_engine_v1_1.py`: shared manuscript-aligned screening engine used by the Python companion.
- `pfas_cast_requirements_v1_1.txt`: Python dependencies.
- `pfas_cast_method_references_v1_1.json`: method and regulatory references.
- `pfas_cast_manifest_v1_1.json`: package checksums.

## Browser tool highlights

- Flexible column detection for system ID, optional location or sample point ID, analyte, result value, detect flag or result sign, and reporting limit. EPA UCMR 5 Data Finder exports (tab-delimited, with `AnalyticalResultValue(µg/L)`, `AnalyticalResultsSign` and `SamplePointID` columns) are recognized without renaming.
- Analyte aliases for common PFAS names.
- Censored-lognormal likelihood fitting with manuscript-style Tier A/B2/B1/C machinery.
- Record-level imputation for censored observations followed by repeated-record aggregation before endpoint screening.
- Screening classifications: always flagged, MI-only screening flag, never flagged and zero-only unexpected.
- Workbook export with systems, endpoint details, fits, import QA, convergence diagnostics, figure-data sheets and references.

## Python companion example

```bash
pip install -r pfas_cast_requirements_v1_1.txt

python pfas_cast_screening_tool_v1_1.py \
  --input monitoring_data.csv \
  --unit ug/L \
  --output pfas_cast_results_v1_1.xlsx \
  --m 500 \
  --seed 42
```

### Demo run

```bash
python pfas_cast_screening_tool_v1_1.py \
  --demo \
  --output demo_pfas_cast_results_v1_1.xlsx \
  --m 100 \
  --seed 42
```

By default, Python figures are written to `<output_stem>_figures` and include PNG, SVG, figure-data CSV files and a `figure_manifest.json`.

## Minimum input columns

The tool can auto-detect common column names. The conceptual input fields are:

- `system_id`
- `location_id` or sampling point (optional; system ID is used if absent)
- `analyte`
- `result_value`
- `detect_flag`
- `reporting_limit`

Input units can be specified as µg/L or ng/L/ppt. All internal calculations use µg/L. Files may be comma- or tab-delimited.

### Using EPA UCMR 5 data directly

Download a state or system subset from the [UCMR 5 Data Finder](https://www.epa.gov/dwucmr/fifth-unregulated-contaminant-monitoring-rule-data-finder) and load it as is. Results with `AnalyticalResultsSign` equal to `<` are treated as below the reporting limit (`MRL`); results with `=` are treated as detections. For the full national file (about 1.9 million rows) use the Python companion rather than the browser tool, and expect run times of minutes rather than seconds.

## Changelog

- **v1.1.1 (October 2026).** Updated citation to the published paper. Added automatic recognition of EPA UCMR 5 Data Finder column names and tab-delimited files in both the browser and Python tools. Fixed two browser-tool defects: a file with no detect-flag column had every row treated as a non-detect (detect status is now inferred from a `<` prefix on the result value or from the result sign column), and a `SamplingLocationID` column was not recognized, so entry points within a system were pooled; the browser and Python tools now give identical zero-substitution flags on the same input. Added the source publication to the method references. No change to the statistical engine.
- **v1.1.0.** Figure outputs and method references added to the Python companion.

## Citation and use statement

When used with the paper, describe PFAS-CAST as an exploratory screening companion implementing the primary marginal censored-likelihood MI engine. Do not describe PFAS-CAST outputs as regulatory-compliance determinations or as physical PFAS mass, flow-weighted load, intake dose or toxicity-weighted risk.

## License

MIT License. © 2025–2026 Gehendra Kharel and contributors.
