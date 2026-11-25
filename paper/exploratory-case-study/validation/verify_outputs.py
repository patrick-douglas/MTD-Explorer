#!/usr/bin/env python3
"""Verify definitive tables without R or reference databases."""
import csv, math
from pathlib import Path
O=Path(__file__).resolve().parents[1]/'outputs'
def rows(p):
    with (O/p).open(newline='') as f: return list(csv.reader(f,delimiter='\t'))
def check(ok, msg):
    if not ok: raise ValueError(msg)
r=rows('temp/bracken_raw_results/bracken_species_all'); cols=[i for i,n in enumerate(r[0]) if n.endswith('_num')]; names=[r[0][i].removeprefix('Report_').removesuffix('.species.bracken_num') for i in cols]; x={a[0]:[float(a[i]) for i in cols] for a in r[1:]}
check(len(x)==104 and len(cols)==6,'matrix dimensions')
sums=[sum(v[i] for v in x.values()) for i in range(len(cols))]
for mode in ('absolute','relative'):
    v={k:([100*a/b for a,b in zip(y,sums)] if mode=='relative' else y) for k,y in x.items()}
    stats={k:[sum(a>0 for a in y)/len(y),math.log10(1+sum(y)/len(y)),math.log10(1+max(y)),math.log10(1+sum(y))] for k,y in v.items()}
    lo=[min(a[i] for a in stats.values()) for i in range(4)]; hi=[max(a[i] for a in stats.values()) for i in range(4)]
    calc={k:.45*a[0]+sum(w*((a[i]-lo[i])/(hi[i]-lo[i]) if hi[i]>lo[i] else 0) for i,w in [(1,.30),(2,.15),(3,.10)]) for k,a in stats.items()}
    t=rows(f'exploratory/taxonomy/detected_microbiome/species_{mode}_with_samples_distribution/detected_microbiome_species_ranked_by_TPI.tsv'); h=t[0]; d=[dict(zip(h,a)) for a in t[1:]]; check(set(calc)=={a['taxon'] for a in d},'TPI taxa'); err=max(abs(calc[a['taxon']]-float(a['tpi_score'])) for a in d); check(err<5.1e-7,'TPI error'); check(all(a['tpi_score']==a['importance_score'] for a in d),'legacy score'); print(f'{mode} TPI: {len(d)} taxa; max error {err:.9g}')
for mode in ('Abundance','relAbundance'):
    for term in ('go','kegg'):
        p=f'hmn_genefamily_abundance_files/humann_genefamilies_{mode}_{term}'
        raw=rows(p+'.tsv'); translated=rows(p+'_translated.tsv'); mp=rows(p+'_translated_id_map.tsv'); check(len(raw)==len(translated)==len(mp),'functional row count'); m=[dict(zip(mp[0],a)) for a in mp[1:]]; orig={a[0]:a[1:] for a in raw[1:]}; out={a[0]:a[1:] for a in translated[1:]}; check(len(orig)==len(out)==len(m),'unique IDs'); check({a['original_id'] for a in m}==set(orig),'original IDs'); check({a['Gene_Family'] for a in m}==set(out),'translated IDs')
        for a in m:
            before=orig[a['original_id']]; after=out[a['Gene_Family']]; check(len(before)==len(after)==6,'sample count'); check(all(math.isclose(float(i),float(j),rel_tol=1e-12,abs_tol=1e-12) for i,j in zip(before,after)),'abundance preservation'); check((a['eligible_for_functional_de']=='TRUE')==bool(__import__('re').fullmatch(r'GO:\d{7}' if term=='go' else r'K\d{5}',a['functional_id'])),'eligibility')
        print(f'{mode} {term}: {len(m)} rows, {sum(a["eligible_for_functional_de"]=="TRUE" for a in m)} eligible')
q=rows('exploratory/taxonomy/microbiome_abundance_qc/species_relative/microbiome_species_sample_qc.tsv'); check({a[0] for a in q[1:]}==set(names) and {a[1] for a in q[1:]}=={'Zebrafish'},'QC metadata'); print('QC sample IDs and group: PASS')
print('All selected output checks passed.')
