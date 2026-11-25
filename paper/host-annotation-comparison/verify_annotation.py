#!/usr/bin/env python3
"""Offline identifier compatibility audit, not a biological accuracy benchmark."""
import argparse, csv, hashlib, json, re, sys
from pathlib import Path
from collections import defaultdict

def readcsv(path, sep=','):
    with path.open(newline='') as f: return list(csv.DictReader(f, delimiter=sep))
def gmt(path):
    d={}
    with path.open() as f:
        for r in csv.reader(f, delimiter='\t'):
            if not r: continue
            assert len(r)>=3 and r[0] not in d, (path,r[0])
            d[r[0]]=set(r[2:])
    return d

def union(d): return set().union(*d.values())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(__doc__); ap.add_argument('--original',type=Path,required=True); ap.add_argument('--explorer',type=Path,required=True); ap.add_argument('--run',type=Path,required=True); ap.add_argument('--reference',type=Path,required=True,help='functional_annotation directory'); ap.add_argument('--out',type=Path,required=True); a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    O,E,B,A=a.original,a.explorer,a.run,a.reference
    paths={'original_catalog':O/'HostSpecies.csv','explorer_catalog':E/'HostSpecies.csv','original_gmt':O/'Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt','host_counts':B/'host_counts.txt','cache':B/'Host_DEG/gene_ID_cache.csv','tpm':B/'Host_DEG/host_counts_TPM.csv','gct':B/'ssGSEA/host.gct','master':A/'custom_taxid_6526_eggNOG_GO_master.gmt','filtered':B/'ssGSEA/custom_gmt_from_host_reference/custom_taxid_6526_eggNOG_GO.gmt','eggnog':A/'taxid_6526_gold_orgdb.emapper.annotations','scores':B/'ssGSEA/ssgsea-results-scores.gct','gene_table':A/'custom_taxid_6526_eggNOG_GO_master_gene_table.tsv','methods':B/'methods/mtd_methods_run_parameters.csv','original_driver':O/'MTD.sh','original_annotation':O/'DEG_Anno_Plot.R','original_builder':O/'Customized_host.sh','original_gct_builder':O/'gct_making.R','original_ssgsea_engine':O/'Tools/ssGSEA2.0/src/ssGSEA2.0.R','original_config':O/'Tools/ssGSEA2.0/config.yaml','explorer_driver':E/'MTD_explorer.sh','explorer_builder':E/'Create_custom_host.sh','explorer_orgdb_builder':E/'aux_scripts/orgdb/build_gold_orgdb_from_gtf_eggnog.R','explorer_gmt_builder':E/'aux_scripts/ssGSEA/build_master_gmt_from_eggnog.py','explorer_gmt_filter':E/'aux_scripts/ssGSEA/filter_master_gmt_for_host_gct.py'}
    old=readcsv(paths['original_catalog']); new=readcsv(paths['explorer_catalog']); oldrow=[r for r in old if r['Taxon_ID']=='6526']; newrow=[r for r in new if r['Taxon_ID']=='6526']; assert not oldrow and len(newrow)==1
    counts=[r['Geneid'] for r in readcsv(paths['host_counts'],'\t')]; cache=[r['ensembl_gene_id'] for r in readcsv(paths['cache'])]; assert len(counts)==len(set(counts))==len(cache)==len(set(cache)) and set(counts)==set(cache)
    universe=set(counts)
    with paths['gct'].open() as f:
        assert next(f).strip()=='#1.2'; dims=next(f).split(); reader=csv.DictReader(f,delimiter='\t'); samples=reader.fieldnames[2:]; h=list(reader)
    genes={r['NAME'] for r in h}; assert len(genes)==len(h)==int(dims[0]) and len(samples)==int(dims[1]); assert genes<=universe
    tpm=readcsv(paths['tpm']); assert len(tpm)==len(h)
    assert all(r['NAME']==t['gene_name'] and r['Description']==t['hybrid_name'] and all(float(r[k])==float(t[k]) for k in samples) for r,t in zip(h,tpm))
    master=gmt(paths['master']); historical=gmt(paths['filtered']); legacy=gmt(paths['original_gmt'])
    # Independent regeneration from the archived eggNOG output using the shared
    # featureCounts/cache universe as reference gene-ID proxy (GTF not in upload).
    rebuilt=defaultdict(set); header=None; nrows=0
    with paths['eggnog'].open() as f:
        for line in f:
            if line.startswith('#query\t'): header=line[1:].rstrip('\n').split('\t'); continue
            if line.startswith('#') or not line.strip(): continue
            assert header is not None; rec=dict(zip(header,line.rstrip('\n').split('\t'))); nrows+=1
            if rec['query'] in universe:
                for term in re.split('[,;]',rec['GOs']):
                    if re.fullmatch(r'GO:\d{7}',term.strip()): rebuilt[term.strip()].add(rec['query'])
    assert dict(rebuilt)==master
    pairs=0
    with paths['gene_table'].open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            assert r['gene_id'] in master[r['go_id']]; pairs+=1
    assert pairs==sum(map(len,master.values()))
    overlap={k:v&genes for k,v in master.items()}; initial={k:v for k,v in overlap.items() if 5<=len(v)<=500}; eligible={k:v for k,v in initial.items() if len(v)>=10}; assert initial==historical
    assert not (genes&union(legacy)); assert not ({g.upper() for g in genes}&{g.upper() for g in union(legacy)})
    with paths['scores'].open() as f:
        assert next(f).strip()=='#1.3'; dims_score=next(f).split(); reader=csv.DictReader(f,delimiter='\t'); scored={}
        for r in reader:
            term=r['id']; assert term not in scored and term in eligible; assert int(r['Signature.set.size'])==len(eligible[term]); assert int(r['No.columns.scored'])==len(samples)
            for s in samples:
                assert set(r['Signature.set.overlap.'+s].split('|'))==eligible[term], (term,s)
            scored[term]=eligible[term]
    assert scored==eligible and int(dims_score[0])==len(scored)
    metrics={'original_catalog_rows':len(old),'explorer_catalog_rows':len(new),'original_taxid_rows':len(oldrow),'explorer_taxid_rows':len(newrow),'reference_count_ids':len(universe),'cache_ids':len(cache),'expression_ids':len(genes),'libraries':len(samples),'legacy_gmt_sets':len(legacy),'legacy_gmt_unique_ids':len(union(legacy)),'legacy_direct_overlap':len(genes&union(legacy)),'legacy_direct_eligible_sets':sum(len(v&genes)>=10 for v in legacy.values()),'legacy_native_GO_coverage':'NA (no configured resource)','eggnog_rows':nrows,'master_GO_sets':len(master),'master_GO_genes':len(union(master)),'master_overlap_genes':len(union(overlap)),'master_overlap_percent':100*len(union(overlap))/len(genes),'filtered_GO_sets_5_500':len(initial),'filtered_GO_genes_5_500':len(union(initial)),'scored_GO_sets_10_500':len(scored),'scored_GO_genes_10_500':len(union(scored)),'scored_GO_percent':100*len(union(scored))/len(genes),'gene_GO_pairs':pairs,'master_rebuilt_from_eggnog':'PASS','historical_GMT_membership':'PASS','historical_score_membership_all_samples':'PASS','TPM_GCT_values':'PASS'}
    with (a.out/'metrics.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['metric','value']); w.writerows(metrics.items())
    with (a.out/'input_sha256.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['role','sha256','basename']); w.writerows((k,sha(v),v.name) for k,v in paths.items())
    with (a.out/'GO_set_audit.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['GO','master_size','expression_overlap','passes_5_500','scored_10_500']); w.writerows((k,len(master[k]),len(v),int(k in initial),int(k in scored)) for k,v in sorted(overlap.items()))
    print(json.dumps(metrics,indent=2)); print('PASS: offline identifier compatibility audit')
if __name__=='__main__': main()
