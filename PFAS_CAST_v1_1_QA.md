# PFAS-CAST v1.1.1 QA note

Checks run during package assembly (22 October 2026).

## Python companion and engine

```bash
python pfas_cast_screening_engine_v1_1.py --demo --output test_engine.xlsx --m 50 --seed 42 --b2-boot 20 --b1-grid 10
python pfas_cast_screening_tool_v1_1.py --demo --output test_tool.xlsx --m 50 --seed 42 --b2-boot 20 --b1-grid 10
python pfas_cast_screening_tool_v1_1.py --input tests/ucmr5_style.txt --output test_ucmr.xlsx --m 30 --seed 42 --b2-boot 10 --b1-grid 10 --no-figures
```

All three completed. On the UCMR 5-style file (tab-delimited, EPA Data Finder column names, `AnalyticalResultsSign`
column, no detect-flag column), all 288 rows parsed, 60 detections were recognized from the sign column, and
`SamplePointID` was used as the location. Version 1.1.0 rejected the same file (required columns not found).

## Browser tool

Headless Chromium (Playwright) loaded `PFAS_CAST_v1_1.html`, uploaded a file, ran screening with 30 replicates, and
read back the parse QA and classification counts.

- `tests/ucmr5_style.txt`: 288 of 288 rows parsed, 60 detections, 12 systems, 24 locations, no page errors.
  Version 1.1.0 produced no result on this file.
- `tests/pfas_cast_demo.csv` (comma-delimited, `SamplingLocationID` column): 204 of 204 rows parsed, 91 detections,
  17 locations. Zero-substitution flags (3 systems) and MI flags (4 systems) match the Python engine on the same file.
  Version 1.1.0 pooled entry points within each system (8 locations) and flagged 2 systems under zero substitution.

The short smoke tests use reduced replicate, bootstrap and grid settings for speed only; defaults remain M = 500,
B2 bootstrap = 200 and B1 grid = 40 x 40. The statistical engine is unchanged from v1.1.0.
