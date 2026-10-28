#!/usr/bin/env python3
"""
PFAS Censoring-Aware Screening Tool (PFAS-CAST), Python companion v1.1.1
Adds figure outputs and method references on top of the screening engine.

Companion to:
  Kharel, G., Kirilenko, A. P., Harvey, O. R. and Zou, C. B. Censoring-aware models quantify systematic underestimation of PFAS occurrence in US drinking-water monitoring. Cell Reports Sustainability (2026). https://doi.org/10.1016/j.crsus.2026.100893
"""
from __future__ import annotations
import argparse, importlib.util, json, sys, time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE_PATH = Path(__file__).with_name('pfas_cast_screening_engine_v1_1.py')
if not BASE_PATH.exists():
    # fall back to /mnt/data for notebook-generated packages
    alt = Path('/mnt/data/pfas_cast_screening_engine_v1_1.py')
    if alt.exists(): BASE_PATH = alt
    else: raise FileNotFoundError('Missing pfas_cast_screening_engine_v1_1.py beside this v1.1 companion.')
_spec = importlib.util.spec_from_file_location('pfas_cast_v1_1_base', BASE_PATH)
_base = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
sys.modules[_spec.name] = _base
_spec.loader.exec_module(_base)

TOOL_VERSION = '1.1.1'
METHOD_REFERENCES = [{'key': 'source publication', 'citation': 'Kharel, G., Kirilenko, A. P., Harvey, O. R. and Zou, C. B. Censoring-aware models quantify systematic underestimation of PFAS occurrence in US drinking-water monitoring. Cell Reports Sustainability (2026). https://doi.org/10.1016/j.crsus.2026.100893', 'doi': '10.1016/j.crsus.2026.100893', 'url': 'https://doi.org/10.1016/j.crsus.2026.100893', 'supports': 'The censoring-aware multiple-imputation framework implemented by this tool.'}, {'key': 'left-censored environmental data', 'citation': 'Helsel, D. R. Statistics for Censored Environmental Data Using Minitab and R. 2nd edn. Wiley (2012).', 'doi': '10.1002/9781118162729', 'url': 'https://doi.org/10.1002/9781118162729', 'supports': 'General treatment of nondetects and censored environmental measurements.'}, {'key': 'substitution vs censored-data methods', 'citation': 'Antweiler, R. C. & Taylor, H. E. Evaluation of statistical treatments of left-censored environmental data using coincident uncensored data sets: I. Summary statistics. Environmental Science & Technology 42, 3732–3738 (2008).', 'doi': '10.1021/es071301c', 'url': 'https://doi.org/10.1021/es071301c', 'supports': 'Empirical comparison of substitution and censored-data approaches for summary statistics.'}, {'key': 'detection-limit multiple imputation', 'citation': 'Lubin, J. H. et al. Epidemiologic evaluation of measurement data in the presence of detection limits. Environmental Health Perspectives 112, 1691–1696 (2004).', 'doi': '10.1289/ehp.7199', 'url': 'https://doi.org/10.1289/ehp.7199', 'supports': 'Model-based handling of measurements subject to detection limits.'}, {'key': 'multiple imputation', 'citation': 'Rubin, D. B. Multiple Imputation for Nonresponse in Surveys. Wiley (1987).', 'doi': '10.1002/9780470316696', 'url': 'https://doi.org/10.1002/9780470316696', 'supports': 'Foundational multiple-imputation framework and repeated-imputation inference.'}, {'key': 'MI interval pooling', 'citation': 'Barnard, J. & Rubin, D. B. Small-sample degrees of freedom with multiple imputation. Biometrika 86, 948–955 (1999).', 'doi': '10.1093/biomet/86.4.948', 'url': 'https://doi.org/10.1093/biomet/86.4.948', 'supports': 'Reference for finite-M multiple-imputation uncertainty calculations.'}, {'key': 'Gaussian copula context', 'citation': 'Hoff, P. D. Extending the rank likelihood for semiparametric copula estimation. Annals of Applied Statistics 1, 265–283 (2007).', 'doi': '10.1214/07-AOAS107', 'url': 'https://doi.org/10.1214/07-AOAS107', 'supports': 'Background for the manuscript’s secondary Gaussian-copula dependence diagnostic; the screening tool uses the primary marginal engine.'}, {'key': 'UCMR 5 data context', 'citation': 'U.S. Environmental Protection Agency. Fifth Unregulated Contaminant Monitoring Rule Data Finder.', 'doi': '', 'url': 'https://www.epa.gov/dwucmr/fifth-unregulated-contaminant-monitoring-rule-data-finder', 'supports': 'UCMR 5 data context and EPA treatment of below-MRL values in MCL-comparison averages.'}, {'key': 'default PFAS screening benchmarks', 'citation': 'U.S. Environmental Protection Agency. Per- and polyfluoroalkyl substances (PFAS) drinking water regulation information.', 'doi': '', 'url': 'https://www.epa.gov/sdwa/and-polyfluoroalkyl-substances-pfas', 'supports': 'Default April 2024 PFAS NPDWR-framework benchmarks used for exploratory screening.'}, {'key': 'regulatory status warning', 'citation': 'U.S. Environmental Protection Agency. EPA announces it will keep maximum contaminant levels for PFOA, PFOS.', 'doi': '', 'url': 'https://www.epa.gov/newsreleases/epa-announces-it-will-keep-maximum-contaminant-levels-pfoa-pfos', 'supports': 'Warning that users should verify current status of PFHxS, PFNA, HFPO-DA and Hazard Index benchmarks before decision support.'}]

SOFT_COLORS = {'purple':'#5d3a8e','purple_light':'#b9a6d3','amber':'#c79a57','green':'#7aa187','red':'#c47b82','blue':'#7f9fc9','slate':'#7b8798','grid':'#e8edf4'}
plt.rcParams.update({'font.size':11,'axes.titlesize':14,'axes.labelsize':12,'xtick.labelsize':10.5,'ytick.labelsize':10.5,'legend.fontsize':10})

def _select_fit_labels(fit_df: pd.DataFrame) -> set[str]:
    keep: list[str] = []
    def add(series):
        for v in series:
            if v not in keep:
                keep.append(v)
    if fit_df.empty:
        return set()
    add(fit_df.sort_values('arithmetic_mean_ug_l', ascending=False)['analyte'].head(6))
    add(fit_df.sort_values('detection_rate_pct', ascending=False)['analyte'].head(2))
    add(fit_df.loc[fit_df['tier'].astype(str).eq('A'), 'analyte'].head(4))
    add(fit_df.sort_values(['arithmetic_mean_ug_l','detection_rate_pct'], ascending=[False,False])['analyte'].head(8))
    return set(keep)


def classification_counts_df(systems_df: pd.DataFrame) -> pd.DataFrame:
    order = ['always_flagged', 'mi_only_screening_flag', 'never_flagged', 'zero_only_flag_unexpected']
    vc = systems_df['classification'].value_counts()
    return pd.DataFrame({'classification': order, 'count': [int(vc.get(k, 0)) for k in order]})

def top_censoring_gap_df(systems_df: pd.DataFrame, n: int = 15) -> pd.DataFrame:
    cols = ['system_id', 'zero_max_ratio', 'mi_mean_max_ratio', 'delta_mean_minus_zero', 'mi_exceedance_probability']
    out = systems_df.loc[systems_df['classification'] == 'mi_only_screening_flag', cols].copy()
    return out.sort_values('delta_mean_minus_zero', ascending=False).head(n)

def fit_overview_df(fits: Dict[str, Any]) -> pd.DataFrame:
    df = _base.fits_to_df(fits).copy()
    keep = ['analyte', 'tier', 'n_total', 'n_detect', 'detection_rate_pct', 'geometric_mean_ug_l', 'arithmetic_mean_ug_l']
    return df[keep].sort_values(['tier', 'analyte'])

def make_figures(figure_dir: Path, systems_df: pd.DataFrame, convergence_df: pd.DataFrame, fits: Dict[str, Any]) -> Dict[str, str]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    out: Dict[str, str] = {}
    class_df = classification_counts_df(systems_df)
    top_df = top_censoring_gap_df(systems_df, n=15)
    fit_df = fit_overview_df(fits)
    conv_df = convergence_df.copy()
    class_df.to_csv(figure_dir / 'figure_data_classification_counts.csv', index=False)
    top_df.to_csv(figure_dir / 'figure_data_top_censoring_gap.csv', index=False)
    fit_df.to_csv(figure_dir / 'figure_data_fit_overview.csv', index=False)
    conv_df.to_csv(figure_dir / 'figure_data_convergence.csv', index=False)
    # Classification counts
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    colors = [SOFT_COLORS['red'], SOFT_COLORS['amber'], SOFT_COLORS['green'], SOFT_COLORS['slate']]
    ax.bar(class_df['classification'], class_df['count'], color=colors[:len(class_df)])
    ax.set_xlabel('Classification'); ax.set_ylabel('Systems'); ax.set_title('System screening classification summary')
    ax.grid(axis='y', color=SOFT_COLORS['grid'], linewidth=1)
    ax.set_axisbelow(True)
    for i, val in enumerate(class_df['count']):
        ax.text(i, val + max(class_df['count'].max()*0.015, 0.05), str(val), ha='center', va='bottom', fontsize=11)
    fig.tight_layout(); p=figure_dir/'figure1_classification_counts.png'; s=figure_dir/'figure1_classification_counts.svg'; fig.savefig(p,dpi=300,bbox_inches='tight'); fig.savefig(s,bbox_inches='tight'); plt.close(fig); out['classification_counts']=str(p)
    # Scatter
    fig, ax = plt.subplots(figsize=(7.4, 6.4))
    cls_color = {'always_flagged':SOFT_COLORS['red'],'mi_only_screening_flag':SOFT_COLORS['amber'],'never_flagged':SOFT_COLORS['green'],'zero_only_flag_unexpected':SOFT_COLORS['slate']}
    cls_label = {'always_flagged':'Always flagged','mi_only_screening_flag':'MI-only','never_flagged':'Never flagged','zero_only_flag_unexpected':'Zero-only unexpected'}
    for cls, sub in systems_df.groupby('classification', dropna=False):
        ax.scatter(sub['zero_max_ratio'], sub['mi_mean_max_ratio'], alpha=0.82, s=34, color=cls_color.get(cls, SOFT_COLORS['slate']), label=cls_label.get(cls, str(cls)))
    maxv=float(np.nanmax(np.r_[systems_df['zero_max_ratio'].to_numpy(), systems_df['mi_mean_max_ratio'].to_numpy(), [1.1]])); lim=max(1.1,1.05*maxv)
    ax.plot([0,lim],[0,lim], color='#9aa6b2', linestyle='--', linewidth=1.2)
    ax.axvline(1.0,ls='--', color=SOFT_COLORS['amber'], linewidth=1.4)
    ax.axhline(1.0,ls='--', color=SOFT_COLORS['amber'], linewidth=1.4)
    ax.set_xlim(0,lim); ax.set_ylim(0,lim)
    ax.grid(color=SOFT_COLORS['grid'], linewidth=1)
    ax.set_axisbelow(True)
    ax.set_xlabel('Zero-substitution maximum benchmark ratio'); ax.set_ylabel('MI mean maximum benchmark ratio'); ax.set_title('System-level MI-minus-zero shift')
    ax.legend(frameon=True, loc='upper left')
    fig.tight_layout(); p=figure_dir/'figure2_zero_vs_mi_scatter.png'; s=figure_dir/'figure2_zero_vs_mi_scatter.svg'; fig.savefig(p,dpi=300,bbox_inches='tight'); fig.savefig(s,bbox_inches='tight'); plt.close(fig); out['zero_vs_mi_scatter']=str(p)
    # Top shifts
    fig, ax = plt.subplots(figsize=(9.6, 6.0))
    if len(top_df):
        ax.bar([str(x) for x in top_df['system_id']], top_df['delta_mean_minus_zero'], color=SOFT_COLORS['purple'])
        ax.tick_params(axis='x', rotation=55, labelsize=9)
        ax.grid(axis='y', color=SOFT_COLORS['grid'], linewidth=1)
        ax.set_axisbelow(True)
    else:
        ax.text(0.5,0.5,'No MI-only screening flags under current settings',ha='center',va='center',transform=ax.transAxes, fontsize=12)
        ax.set_xticks([])
    ax.set_xlabel('System'); ax.set_ylabel('Δ mean ratio (MI − zero)'); ax.set_title('Top MI-only screening flags by MI-minus-zero shift')
    fig.tight_layout(); p=figure_dir/'figure3_top_censoring_gap_systems.png'; s=figure_dir/'figure3_top_censoring_gap_systems.svg'; fig.savefig(p,dpi=300,bbox_inches='tight'); fig.savefig(s,bbox_inches='tight'); plt.close(fig); out['top_censoring_gap_systems']=str(p)
    # Fit overview
    fig, ax = plt.subplots(figsize=(9.0, 6.4))
    tier_colors = {'A':SOFT_COLORS['purple'],'B2':SOFT_COLORS['blue'],'B1':SOFT_COLORS['amber'],'C':SOFT_COLORS['slate']}
    for tier, sub in fit_df.groupby('tier', dropna=False):
        ax.scatter(sub['detection_rate_pct'], sub['arithmetic_mean_ug_l'], label=f'Tier {tier}', s=48, alpha=0.9, color=tier_colors.get(str(tier), SOFT_COLORS['slate']))
    if len(fit_df):
        x_max = max(5.0, float(fit_df['detection_rate_pct'].max()) * 1.10)
        y_max = max(0.005, float(fit_df['arithmetic_mean_ug_l'].max()) * 1.10)
        ax.set_xlim(0, x_max)
        ax.set_ylim(0, y_max)
        ax.grid(color=SOFT_COLORS['grid'], linewidth=1)
        ax.set_axisbelow(True)
        labels = _select_fit_labels(fit_df)
        placed = []
        rows = fit_df[fit_df['analyte'].isin(labels)].sort_values(['arithmetic_mean_ug_l','detection_rate_pct'], ascending=[False, True]).itertuples(index=False)
        for i, row in enumerate(rows):
            x = float(row.detection_rate_pct)
            y = float(row.arithmetic_mean_ug_l)
            dx = 8 if x < 0.82 * x_max else -8
            dy = -2 if i % 2 == 0 else 11
            tx = x + dx * x_max / 200.0
            ty = y + dy * y_max / 160.0
            tx = min(max(tx, 0.02 * x_max), 0.98 * x_max)
            ty = min(max(ty, 0.03 * y_max), 0.98 * y_max)
            for _ in range(40):
                moved = False
                for px, py in placed:
                    if abs(tx - px) < 0.11 * x_max and abs(ty - py) < 0.05 * y_max:
                        ty = min(0.98 * y_max, py + 0.055 * y_max)
                        moved = True
                if not moved:
                    break
            placed.append((tx, ty))
            ha = 'left' if tx >= x else 'right'
            ax.annotate(str(row.analyte), xy=(x, y), xytext=(tx, ty), textcoords='data', ha=ha, va='center', fontsize=9,
                        arrowprops=dict(arrowstyle='-', lw=0.8, color='#cfd8e3', shrinkA=0, shrinkB=3))
    ax.set_xlabel('Detection rate (%)'); ax.set_ylabel('Fitted arithmetic mean (µg/L)'); ax.set_title('Analyte-level censored-likelihood fit overview')
    ax.legend(frameon=True, fontsize=10, loc='upper right')
    fig.tight_layout(); p=figure_dir/'figure4_fit_overview.png'; s=figure_dir/'figure4_fit_overview.svg'; fig.savefig(p,dpi=300,bbox_inches='tight'); fig.savefig(s,bbox_inches='tight'); plt.close(fig); out['fit_overview']=str(p)
    # Convergence
    fig, ax = plt.subplots(figsize=(8.0,5.1))
    ycol = 'mi_only_screening_flag' if 'mi_only_screening_flag' in conv_df.columns else 'newly_flagged'
    xcol = 'replicates' if 'replicates' in conv_df.columns else 'M'
    ax.plot(conv_df[xcol], conv_df[ycol], marker='o', color=SOFT_COLORS['purple'], linewidth=2.2, markersize=5)
    ax.grid(color=SOFT_COLORS['grid'], linewidth=1)
    ax.set_axisbelow(True)
    ax.set_xlabel('MI replicates'); ax.set_ylabel('MI-only screening flags'); ax.set_title('MI replicate sensitivity / convergence')
    fig.tight_layout(); p=figure_dir/'figure5_convergence.png'; s=figure_dir/'figure5_convergence.svg'; fig.savefig(p,dpi=300,bbox_inches='tight'); fig.savefig(s,bbox_inches='tight'); plt.close(fig); out['convergence']=str(p)
    pd.DataFrame(METHOD_REFERENCES).to_csv(figure_dir/'method_references.csv', index=False)
    with (figure_dir/'figure_manifest.json').open('w', encoding='utf-8') as f:
        json.dump({'tool_version': TOOL_VERSION, 'generated_utc': datetime.now(timezone.utc).isoformat(), 'figures': out, 'method_references': METHOD_REFERENCES}, f, indent=2)
    return out

def write_outputs_v1(out_xlsx: Path, df: pd.DataFrame, qa_df: pd.DataFrame, fits: Dict[str, Any], benchmarks: Dict[str, Dict[str, float]], results: Dict[str, pd.DataFrame], config: Dict[str, Any]) -> None:
    _base.write_outputs(out_xlsx, df, qa_df, fits, benchmarks, results, config)
    with pd.ExcelWriter(out_xlsx, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        classification_counts_df(results['systems']).to_excel(writer, sheet_name='Figure_ClassCounts', index=False)
        top_censoring_gap_df(results['systems'], n=25).to_excel(writer, sheet_name='Figure_TopCensoringShift', index=False)
        fit_overview_df(fits).to_excel(writer, sheet_name='Figure_FitOverview', index=False)
        results['convergence'].to_excel(writer, sheet_name='Figure_Convergence', index=False)
        pd.DataFrame(METHOD_REFERENCES).to_excel(writer, sheet_name='References', index=False)
        for ws in writer.book.worksheets:
            ws.freeze_panes = 'A2'
            for col in ws.columns:
                max_len = max((len(str(cell.value)) if cell.value is not None else 0) for cell in col[:200])
                ws.column_dimensions[col[0].column_letter].width = min(max(max_len + 2, 10), 48)

def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description='PFAS Censoring-Aware Screening Tool (PFAS-CAST), Python companion v1.1.1 with figures and references')
    p.add_argument('--input', type=Path, help='Input CSV. Omit when using --demo.')
    p.add_argument('--output', type=Path, default=Path('pfas_cast_screening_results_v1_1.xlsx'))
    p.add_argument('--unit', default='ug/L')
    p.add_argument('--m', type=int, default=500)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--benchmarks', type=Path)
    p.add_argument('--flag-rule', choices=['mean','probability'], default='mean')
    p.add_argument('--prob-threshold', type=float, default=0.5)
    p.add_argument('--b2-boot', type=int, default=200)
    p.add_argument('--b1-grid', type=int, default=40)
    p.add_argument('--demo', action='store_true')
    p.add_argument('--demo-csv', type=Path, default=Path('pfas_cast_demo.csv'))
    p.add_argument('--figure-dir', type=Path)
    p.add_argument('--no-figures', action='store_true')
    p.add_argument('--system-col', dest='system_id')
    p.add_argument('--location-col', dest='location_id')
    p.add_argument('--analyte-col', dest='analyte')
    p.add_argument('--value-col', dest='value')
    p.add_argument('--detect-col', dest='detect')
    p.add_argument('--rl-col', dest='reporting_limit')
    args = p.parse_args(argv)
    if args.demo:
        input_path = _base.make_demo_data(args.demo_csv, seed=args.seed); print(f'Generated demo CSV: {input_path}')
    elif args.input:
        input_path = args.input
    else:
        p.error('Provide --input CSV or use --demo.'); return 2
    t0=time.time()
    explicit_cols={k:getattr(args,k) for k in _base.COLUMN_ALIASES}
    df, qa_df, qa = _base.parse_input_csv(input_path, args.unit, explicit_cols)
    benchmarks = _base.load_benchmarks(args.benchmarks)
    print(f'Loaded {len(df):,} valid records from {input_path}')
    print(f"Systems={df['system_id'].nunique():,}; analytes={df['analyte'].nunique():,}; detections={int(df['detect'].sum()):,}")
    print('Fitting censored lognormal marginals...')
    fits = _base.fit_all_analytes(df, seed=args.seed, b2_boot=args.b2_boot, b1_grid=args.b1_grid)
    print('Tier counts:', pd.Series([f.tier for f in fits.values()]).value_counts().to_dict())
    print(f'Running record-level MI screening with M={args.m}...')
    results = _base.run_screening(df, fits, benchmarks, args.m, args.seed, args.flag_rule, args.prob_threshold)
    config={'tool_version': TOOL_VERSION, 'base_engine_version': getattr(_base,'TOOL_VERSION','1.0.0'), 'input_file': str(input_path), 'input_sha256': _base.sha256_file(input_path), 'unit_input': args.unit, 'unit_internal': getattr(_base,'INTERNAL_UNIT','ug/L'), 'm': args.m, 'seed': args.seed, 'flag_rule': args.flag_rule, 'prob_threshold': args.prob_threshold, 'benchmarks': benchmarks, 'b2_boot': args.b2_boot, 'b1_grid': args.b1_grid, 'runtime_note': 'Screening and prioritization tool; not a regulatory-compliance determination.', 'method_references': METHOD_REFERENCES}
    write_outputs_v1(args.output, df, qa_df, fits, benchmarks, results, config)
    if not args.no_figures:
        figure_dir = args.figure_dir or args.output.with_suffix('').with_name(args.output.with_suffix('').name + '_figures')
        figs = make_figures(figure_dir, results['systems'], results['convergence'], fits)
        print(f'Figure directory: {figure_dir}')
        for k,v in figs.items(): print(f'  - {k}: {v}')
    print('\nClassification summary:')
    print(results['systems']['classification'].value_counts().to_string())
    print(f'\nWrote: {args.output}')
    print(f'Runtime: {(time.time()-t0):.1f} seconds')
    print('\nRegulatory note: default benchmarks are configurable screening thresholds based on the April 2024 PFAS NPDWR framework; verify current status before decision support.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
