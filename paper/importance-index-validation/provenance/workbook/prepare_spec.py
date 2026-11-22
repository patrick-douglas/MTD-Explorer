from pathlib import Path
import csv
import json
import math

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'results'
WORK = Path(__file__).resolve().parent
summary = json.loads((RESULTS / 'summary.json').read_text())
protocol = json.loads((RESULTS.parent / 'protocol.json').read_text())

LABELS = {
    'representation': 'Abundance representation', 'scenario_id': 'Scenario', 'family': 'Perturbation family',
    'weight_prevalence': 'Prevalence weight', 'weight_mean': 'Mean weight',
    'weight_maximum': 'Maximum weight', 'weight_total': 'Total weight',
    'eligible_taxa': 'Eligible taxa', 'spearman_all': 'Full-ranking Spearman rho',
    'top10_retained': 'Reference top 10 retained', 'top10_jaccard': 'Top 10 Jaccard',
    'top20_retained': 'Reference top 20 retained', 'top20_jaccard': 'Top 20 Jaccard',
    'top50_retained': 'Reference top 50 retained', 'top50_jaccard': 'Top 50 Jaccard',
    'median_abs_shift_original_top20': 'Median absolute rank shift of reference top 20',
    'detail': 'Scenario definition', 'model': 'Score model', 'taxon': 'Sequence-classified taxon label',
    'reference_rank': 'Full-sample model rank', 'historical_default_rank': 'Historical default rank',
    'eligibility_frequency': 'Bootstrap eligibility frequency',
    'top10_frequency': 'Bootstrap top 10 frequency', 'top20_frequency': 'Bootstrap top 20 frequency',
    'top50_frequency': 'Bootstrap top 50 frequency',
    'rank_p025_given_eligible': 'Eligible rank 2.5th percentile',
    'rank_median_given_eligible': 'Eligible rank median',
    'rank_p975_given_eligible': 'Eligible rank 97.5th percentile',
    'score_median_given_eligible': 'Eligible score median',
    'omitted_sample': 'Omitted library', 'omitted_group': 'Omitted library group',
    'historical_default_top20_retained': 'Historical default top 20 retained',
    'eligible_omissions': 'Omissions with taxon eligible (of 15)',
    'top20_omissions': 'Omissions with taxon in top 20 (of 15)',
    'rank_min_given_eligible': 'Minimum eligible rank', 'rank_max_given_eligible': 'Maximum eligible rank',
    'weight_grid_top20_frequency': 'Weight-grid top 20 scenario proportion',
    'weight_grid_rank_min': 'Minimum weight-grid rank', 'weight_grid_rank_max': 'Maximum weight-grid rank',
}
TEXT = {'representation', 'scenario_id', 'family', 'detail', 'model', 'taxon', 'omitted_sample', 'omitted_group'}

def read_table(filename, order):
    with (RESULTS / filename).open(newline='') as stream:
        reader = csv.DictReader(stream)
        assert set(reader.fieldnames) == set(order), (filename, reader.fieldnames, order)
        rows = list(reader)
    values = []
    for row in rows:
        converted = []
        for key in order:
            raw = row[key]
            if key in TEXT:
                converted.append("'" + raw if raw.startswith('=') else raw)
            elif raw == '' or raw.lower() in {'nan', 'none', 'na'}:
                converted.append('n.a.')
            else:
                x = float(raw)
                assert math.isfinite(x), (filename, key, raw)
                converted.append(int(x) if x.is_integer() else x)
        values.append(converted)
    types, formats, widths, wraps = [], {}, [], []
    for i, key in enumerate(order):
        if key in TEXT:
            types.append('text')
        elif 'frequency' in key:
            types.append('proportion')
            formats[i] = '0.0%'
        elif key in {'weight_grid_rank_min', 'weight_grid_rank_max'}:
            types.append('integer')
        elif any(t in key for t in ('jaccard', 'spearman', 'weight_', 'score_', 'percentile')) or key.startswith('rank_p') or key.startswith('rank_median') or key.startswith('median_abs'):
            types.append('number')
            formats[i] = '0.000' if key not in {'rank_median_given_eligible', 'rank_p025_given_eligible', 'rank_p975_given_eligible', 'median_abs_shift_original_top20'} else '0.0'
        else:
            types.append('integer')
        if key == 'taxon':
            widths.append(68); wraps.append(i)
        elif key == 'detail':
            widths.append(68); wraps.append(i)
        elif key == 'representation': widths.append(22)
        elif key == 'scenario_id': widths.append(26)
        elif key == 'family': widths.append(23); wraps.append(i)
        elif key == 'model': widths.append(17)
        elif key == 'omitted_sample': widths.append(25)
        elif key == 'omitted_group': widths.append(23)
        elif key.startswith('weight_') and key not in {'weight_grid_rank_min','weight_grid_rank_max','weight_grid_top20_frequency'}: widths.append(15)
        elif len(LABELS[key]) > 35: widths.append(22)
        elif len(LABELS[key]) > 24: widths.append(19)
        else: widths.append(17)
    return values, types, widths, wraps, formats

def data_sheet(name, title, context, filename, order):
    rows, types, widths, wraps, formats = read_table(filename, order)
    return {'name': name, 'title': title, 'context': context, 'headers': [LABELS[k] for k in order],
            'rows': rows, 'types': types, 'widths': widths, 'wrapColumns': wraps,
            'formats': formats, 'freezeColumns': 3 if 'taxon' in order and order.index('taxon') == 2 else 2,
            'sourceFile': filename, 'sourceColumns': order, 'headerHeight': 58}

summary_rows = []
measures = {
    'top20_retained': 'Within-model top 20 retained (count)',
    'top20_jaccard': 'Top 20 Jaccard (0–1)',
    'spearman_all': 'Full-ranking Spearman rho',
    'median_abs_shift_original_top20': 'Reference top 20 median absolute rank shift',
    'eligible_taxa': 'Eligible taxa (count)',
}
for key, block in summary['sampling_summaries'].items():
    representation, model = key.split('/')
    for measure, label in measures.items():
        stats = block[measure]
        summary_rows.append([representation, model, label] + [stats[x] for x in ['min','p025','median','p975','max']])
sheets = [{
    'name': 'S9 Summary', 'title': 'Supplementary Table S9  Importance index robustness',
    'context': 'Bootstrap comparisons use each model’s own full-sample ranking. 5,000 stratified draws per model; 425 tracked features and 15 libraries.',
    'headers': ['Abundance representation', 'Score model', 'Bootstrap measure', 'Minimum', '2.5th percentile', 'Median', '97.5th percentile', 'Maximum'],
    'rows': summary_rows, 'types': ['text','text','text'] + ['number']*5,
    'widths': [22,17,48,16,18,16,18,16], 'wrapColumns': [2], 'freezeColumns': 2, 'headerHeight': 45,
}]

sheets.append(data_sheet('S9A Scenarios', 'S9A  Weight and detection sensitivity',
    'All 93 specified scenarios per representation; comparisons use the historical default ranking. Weights sum to one.',
    'scenario_metrics.csv', ['representation','scenario_id','family','weight_prevalence','weight_mean','weight_maximum','weight_total','eligible_taxa','top20_retained','top20_jaccard','median_abs_shift_original_top20','spearman_all','top10_retained','top10_jaccard','top50_retained','top50_jaccard','detail']))
for key in ['top10_retained', 'top20_retained', 'top50_retained']:
    i = sheets[-1]['sourceColumns'].index(key)
    sheets[-1]['headers'][i] = LABELS[key].replace('Reference', 'Historical-default')
    sheets[-1]['widths'][i] = 22
sheets.append(data_sheet('S9B Bootstrap taxa', 'S9B  Taxon stability under stratified resampling',
    '5,000 draws per model; frequencies use all draws. Rank and score summaries are conditional on eligibility.',
    'bootstrap_taxon_summary.csv', ['representation','model','taxon','reference_rank','historical_default_rank','eligibility_frequency','top20_frequency','rank_p025_given_eligible','rank_median_given_eligible','rank_p975_given_eligible','top10_frequency','top50_frequency','score_median_given_eligible']))
sheets.append(data_sheet('S9C LOO metrics', 'S9C  Leave one library out',
    'Primary comparisons use each model’s own full-sample ranking. All 15 omissions; the omitted group has four libraries and other groups retain five.',
    'loo_metrics.csv', ['representation','model','omitted_sample','omitted_group','eligible_taxa','top20_retained','top20_jaccard','historical_default_top20_retained','median_abs_shift_original_top20','spearman_all','top10_retained','top10_jaccard','top50_retained','top50_jaccard']))
for key in ['top10_retained', 'top20_retained', 'top50_retained']:
    i = sheets[-1]['sourceColumns'].index(key)
    sheets[-1]['headers'][i] = LABELS[key].replace('Reference', 'Within-model')
    sheets[-1]['widths'][i] = 22
sheets.append(data_sheet('S9D Weight grid taxa', 'S9D  Taxon ranks across the joint weight grid',
    '79 unique normalized weight vectors from 81 joint multiplier settings. Scenario proportions are not probabilities.',
    'weight_grid_taxon_summary.csv', ['representation','taxon','reference_rank','weight_grid_top20_frequency','weight_grid_rank_min','weight_grid_rank_max']))
sheets.append(data_sheet('S9E LOO taxa', 'S9E  Taxon stability to individual library omissions',
    'Counts have denominator 15 omissions. Rank bounds include only omissions in which the taxon was eligible.',
    'loo_taxon_summary.csv', ['representation','model','taxon','reference_rank','historical_default_rank','eligible_omissions','top20_omissions','rank_min_given_eligible','rank_max_given_eligible']))

definitions = [
    ['Scope', protocol['scope']],
    ['Dataset', protocol['dataset']],
    ['Representations', 'relative_percent: within-library percentages of species-level Bracken estimates; estimated_reads: Bracken estimated read counts. Neither is an absolute organism-load measurement.'],
    ['Tracking universe', '425 sequence-classified taxonomic features. The supplied matrices use names rather than numeric TaxIDs. All-zero features are ineligible within a resample.'],
    ['Default score', protocol['score_rule']],
    ['Default component weights', 'Prevalence 0.45; mean 0.30; maximum 0.15; total 0.10.'],
    ['No-total companion model', 'Remove total and renormalize the remaining default weights to prevalence 0.50, mean 1/3 and maximum 1/6.'],
    ['Component dependence', 'For a fixed sample count, total abundance equals sample count times mean abundance. Log transformation and min–max scaling do not make these summaries independent.'],
    ['Historical reproduction', summary['historical_reproduction'] + '. This verifies implementation fidelity, not biological validity.'],
    ['Joint weight grid', protocol['weight_sensitivity']['design'] + ' Exactly 79 unique normalized vectors were evaluated.'],
    ['Named comparisons', ', '.join(protocol['weight_sensitivity']['named_comparisons'])],
    ['Scenario reference', 'All S9A comparisons and S9D ranks use the historical default ranking as reference. Named baselines quantify ordering differences; they do not establish a superior method.'],
    ['Presence sensitivity', 'Use estimated-read thresholds of at least 10 or 50 to define presence only. Keep abundance summaries and the original normalization universe unchanged. These are sensitivity settings, not validated detection cutoffs.'],
    ['10-read scenario', 'This setting leaves presence unchanged because the smallest positive estimated-read entry is at least 10. A no-op scenario does not demonstrate robustness to stricter detection.'],
    ['Universe sensitivity', 'Retain taxa observed in at least 2 or 3 original libraries. Recompute min–max over retained taxa. Keep original abundance values and relative denominators; do not reclose percentages.'],
    ['Bootstrap design', '5,000 draws of five libraries with replacement within each of infected, infection_failed and uninfected. Identical draws are used for both representations and both score models.'],
    ['Bootstrap random number generator', 'NumPy default_rng with PCG64; seed 20261120.'],
    ['Bootstrap calculations', protocol['sampling_stability']['processing']],
    ['Resampling reference', 'S9 Summary, S9B and the primary S9C metrics compare each model with its own full-sample ranking. Historical-default rank/top20 fields are supplied separately.'],
    ['Bootstrap frequencies', 'Top 10/20/50 and eligibility frequencies use all 5,000 draws. A taxon that is absent or ineligible contributes zero to inclusion frequencies.'],
    ['Conditional ranks', 'Per-taxon rank quantiles and rank bounds use eligible draws/omissions only. Inspect eligibility frequency/count beside the rank summaries.'],
    ['Percentile interpretation', 'Empirical resampling percentiles summarize stability conditional on the observed libraries. They are not population confidence intervals or additional independent biological replication.'],
    ['LOO design', 'All 15 single-library omissions; group sizes become 4/5/5. Do not describe LOO as maintaining equal group counts.'],
    ['Full-ranking Spearman rho', protocol['metrics']['excluded_rank_for_full_comparison'] + ' The resulting statistic covers all 425 tracked features.'],
    ['Top-k retention', 'Number of features shared with the reference top k. Jaccard is intersection divided by union of the selected sets.'],
    ['Rank-shift measure', 'Median absolute displacement of the original reference top 20 positions. Consult the accompanying per-taxon eligibility information.'],
    ['Perturbation families', 'Report the joint grid separately from named ablations, presence rules and universe filters. A grid-specific range does not describe all scenarios. All specified adverse cases are retained in S9A.'],
    ['Interpretive limit', 'Internal sensitivity and stability in this dataset; no optimized weights, external validation, confirmed taxonomic identity, biological importance, or accuracy superiority is established.'],
    ['Software', 'Python ' + summary['software']['python'] + '; NumPy ' + summary['software']['numpy']],
    ['Protocol SHA-256', summary['protocol_sha256']],
    ['Analysis script SHA-256', summary['analysis_script_sha256']],
    ['Historical benchmark commit', protocol['benchmark_commit']],
    ['Versioned ranking source', 'https://github.com/patrick-douglas/MTD-Explorer/blob/' + protocol['benchmark_commit'] + '/aux_scripts/exploratory/MTD_exploratory_detected_microbiome_rank.R'],
    ['Methodological reference', 'Saisana M, Saltelli A, Tarantola S (2005). Uncertainty and sensitivity analysis techniques as tools for the quality assessment of composite indicators. JRSS A168(2):307–323. https://doi.org/10.1111/j.1467-985X.2005.00350.x'],
    ['Complete records', 'The reproducibility bundle retains protocol.json, input matrices, weight_design.csv, all scenario/LOO summaries, bootstrap draw identities and replicate-level records.'],
]
for source in summary['inputs']:
    definitions.append(['Input ' + source['file'], 'SHA-256 ' + source['sha256'] + '; ' + str(source['bytes']) + ' bytes'])
sheets.append({'name': 'S9F Definitions', 'title': 'S9F  Definitions and provenance',
    'context': 'Definitions apply to all S9 parts. Table values are static outputs of the archived Python analysis.',
    'headers': ['Item', 'Definition or source'], 'rows': definitions,
    'types': ['text','text'], 'widths': [42,126], 'wrapColumns': [0,1], 'headerHeight': 30})

spec = {'outputName': 'Supplementary_Table_S9.xlsx', 'sheets': sheets}
(WORK / 'workbook_spec.json').write_text(json.dumps(spec, ensure_ascii=False, indent=2))
print(json.dumps({'sheets': [{'name': s['name'], 'rows': len(s['rows']), 'columns':len(s['headers'])} for s in sheets]},indent=2))
