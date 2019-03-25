# main REST api                           https://rest.ensembl.org/?content-type=text/html
# training from ensembl                   https://notebooks.azure.com/ensembl-training/projects/rest-api-r
# GTEX                                  http://www.ensembl.info/2016/08/15/gtex-eqtl-data-now-in-ensembl/
# core api structure and db scheme       https://www.ensembl.org/info/docs/api/core/index.html#api

# setup ----------------------------------------------------------------------------------------------------------------

import pandas as pd
from dfply import *
import pybiomart as pybm
import re

import os
import gzip
import glob
from pprint import pprint
import subprocess
import multiprocessing as mp
import pickle
import json

# import own functions and classes
from Ensembl_REST_client import *
import uniprot_id_mapper

humanchromosomes = [str(x) for x in range(1,23)] + ['X', 'Y', 'MT']

# read input data ------------------------------------------------------------------------------------------------------
mrm_apos = pd.read_csv("../input_data/mrm_apo_peptides.csv")

# map to ensembl IDs ---------------------------------------------------------------------------------------------------
mrm_apos_id_map = uniprot_id_mapper.uniprot_id_mapper(ids=mrm_apos.Protein, origindb='ACC+ID', targetdb ='ENSEMBL_ID')
mrm_apos_id_map = mrm_apos_id_map.reset_index(drop=True)

# save output
mrm_apos_id_map.to_csv('../data/mrm_apos_id_map.csv', sep='\t', index=False)
mrm_apos_id_map = pd.read_csv('../data/mrm_apos_id_map.csv', sep='\t')

# confirm that all uniprot identifiers were mapped to ensembl IDs
set(mrm_apos.Protein) - set(mrm_apos_id_map.From)

# ensembl biomart lookup -----------------------------------------------------------------------------------------------

server = pybm.Server(host="http://www.ensembl.org")
ens_snp = (server.marts['ENSEMBL_MART_SNP']
                 .datasets['hsapiens_snp'])
ens_genes = (server.marts['ENSEMBL_MART_ENSEMBL']
                 .datasets['hsapiens_gene_ensembl'])
mrm_apos_bm_38 = (ens_genes.query(
    attributes=['ensembl_gene_id','chromosome_name', 'start_position', 'end_position', 'strand'],
    filters={'link_ensembl_gene_id': list(mrm_apos_id_map.To)}) >>
                      rename(id='Gene stable ID',
                             chr='Chromosome/scaffold name',
                             start='Gene start (bp)',
                             end='Gene end (bp)',
                             strand='Strand'))
mrm_apos_bm_canonical_38 = mrm_apos_bm_38[mrm_apos_bm_38.chr.isin(humanchromosomes)]

# # canonical ens IDs are identical between GRCh37 and GRCh38, only positions change.
# ens_genes37 = pybm.Dataset(name='hsapiens_gene_ensembl', host='http://grch37.ensembl.org/')
# mrm_apos_bm = (ens_genes37.query(
#     attributes=['ensembl_gene_id','chromosome_name', 'start_position', 'end_position', 'strand'],
#     filters={'link_ensembl_gene_id': list(mrm_apos_id_map.To)}) >>
#                       rename(id='Gene stable ID',
#                              chr='Chromosome/scaffold name',
#                              start='Gene start (bp)',
#                              end='Gene end (bp)',
#                              strand='Strand'))
# mrm_apos_bm_noncanonical = mrm_apos_bm[~mrm_apos_bm.chr.isin(humanchromosomes)]
# mrm_apos_bm_canonical = mrm_apos_bm[mrm_apos_bm.chr.isin(humanchromosomes)]






# mapping relations in
# 1) uniprot mapping,
# 2) ensembl biomart annotation, and

print(mrm_apos.shape)
print(mrm_apos_id_map.shape)
print(mrm_apos_bm_38.shape)
print(mrm_apos_bm_canonical_38.shape)

# save output
mrm_apos_bm_canonical_38.to_csv('../data/mrm_apos_bm_canonical_38.tsv', sep='\t', index=False)
mrm_apos_bm_canonical_38 = pd.read_csv('../data/mrm_apos_bm_canonical_38.tsv', sep='\t')


# Ensembl REST variants --------------------------------------------------------------------------------------
mrm_apos_rs = {}
mrm_apos_snv = {}

for ens_id in mrm_apos_id_map.To:
    rs, snv = run_variants_by_id(ens_id, feature=['variation', 'structural_variation'])
    mrm_apos_rs[ens_id] = rs
    mrm_apos_snv[ens_id] = snv

# population documentation: https://www.ensembl.org/info/genome/variation/species/populations.html


# save output
with open('../data/mrm_apos_rs.json', 'w') as rs_f:
    rs_f.write(json.dumps(mrm_apos_rs))
with open('../data/mrm_apos_snv.json', 'w') as snv_f:
    snv_f.write(json.dumps(mrm_apos_snv))

with open('../data/mrm_apos_rs.json', 'r') as rs_f:
    mrm_apos_rs = json.load(rs_f)
with open('../data/mrm_apos_snv.json', 'r') as snv_f:
    mrm_apos_snv = json.load(snv_f)


# Ensembl VEP local ----------------------------------------------------------------------------------------------------

# some setup
rs_dir = '../data/mrm_apos_rs/'
vep_out_dir = '../data/vep_out'
config_fp = 'mrm-vep-config.txt'

# helper function
def run_vep_ens_id(ens_id):
    print(ens_id)
    rs_fp = '{0}/{1}_mrm_apos_rs.tsv'.format(rs_dir, ens_id)

    with open(rs_fp, "w") as f:
        for rs in mrm_apos_rs.get(ens_id):
            # for rs in mrm_apos_rs.get(ens_id)[:5]:
            f.write("{0}\n".format(rs))

    vep_out = '{0}/{1}'.format(vep_out_dir, ens_id)
    if os.path.isfile(vep_out):
        print(ens_id + ' already retrieved, so skipped.')
    else:
        subprocess.run(
            '/home/simon/ensembl-vep/vep -i {0} --config {1} -o {2} --format "id" --cache  --offline --verbose --protein --tsl --biotype --variant_class --tab --regulatory --show_ref_allele'.format(
                rs_fp, config_fp, vep_out),
            check=True, text=True, shell=True)

# run VEP locally
subprocess.run("mkdir -p {0}".format(rs_dir), shell=True, check=True)
subprocess.run("mkdir -p {0}".format(vep_out_dir), shell=True, check=True)
with mp.Pool(processes=3) as pool:
    pool.map(run_vep_ens_id, mrm_apos_bm_canonical_38.id.values)

# define coding consequence types
coding_variant_types = (['missense_variant', 'splice_region_variant', 'inframe_deletion', 'frameshift_variant',
                        'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
                        'stop_gained', 'stop_lost', 'start_lost', 'transcript_amplification', 'inframe_insertion',
                        'protein_altering_variant', 'incomplete_terminal_codon_variant'])

# parse VEP output, keeping only coding variants
vep_coding = {}
for ens_id in mrm_apos_bm_canonical_38.id.values:

    #debug
    # ens_id = 'ENSG00000110244'
    print(ens_id)

    vep_fp = '{0}/{1}'.format(vep_out_dir, ens_id)
    with open(vep_fp, "r") as vep_f:
        # vep_output = read_vep_tab(vep_f)

        # skip header lines
        while True:
            line = vep_f.readline()
            if not re.match("##", line):
                headerline = [x.lower().replace(' ', '_') for x in line.split(sep='\t')]
                headerline[0] = headerline[0].replace('#', '')
                break
        vep_output = pd.read_csv(vep_f, sep='\t', names=headerline)

        # separate rows of vep_output.consequence
        old_cols = list(vep_output.columns)

        tmp = (vep_output.consequence.str.split(',')
               .apply(pd.Series)
               .merge(vep_output, left_index=True, right_index=True)
               .melt(id_vars=old_cols, value_name='conseqs', var_name='dummy_var')
               .drop(columns=['consequence', 'dummy_var']))
        tmp = tmp.loc[tmp.conseqs.notna()]
        tmp.sort_values(by=['uploaded_variation', 'allele', 'feature'], inplace=True)

        vep_coding[ens_id] = tmp.loc[(tmp.conseqs.isin(coding_variant_types) &
                                      tmp.uploaded_variation.str.contains('rs') &
                                      (tmp['gene'] == ens_id) &
                                      ((tmp['tsl'] == '1') | (tmp['tsl'] == '2')))]

# save output
pickle.dump(vep_coding, open('../data/vep_coding.p', "wb"))
vep_coding = pickle.load( open('../data/vep_coding.p', "rb"))

# Ensembl REST MAFs ----------------------------------------------------------------------------------------------------

# get MAFs for coding variants
# mrm_coding_mafs = {}
mafs_dir = '../data/mafs_REST'
subprocess.run("mkdir -p {0}".format(mafs_dir), shell=True, check=True)
for ens_id in vep_coding.keys():

    #debug
    # ens_id = 'ENSG00000110244'
    print(ens_id)

    mafs_fp = '{0}/{1}'.format(mafs_dir, ens_id)
    if os.path.isfile(mafs_fp):
        print(ens_id + ' already retrieved, so skipped.')
    else:
        mrm_coding_mafs = run_variant_mafs(rsids=vep_coding[ens_id].uploaded_variation.to_list())
        mrm_coding_mafs.to_csv(mafs_fp, sep='\t')




# merge MAFs back to VEP and LD-expand -------------------------------------------------------------


# approach:
#   1. annotate major_allele column in mafs, calculated per rs and population
#   2. conditional join, implemented as
#       2a. merge on rs only
#       2b. filter on (mafs.major_allele == vep.ref_allele and mafs.pop_allele == vep.allele)
#                   OR (mafs.major_allele == vep.allele and mafs.pop_allele == vep.ref_allele)

# FIRST: check assumptions of approach:
#   1. ref allele is among 2 most frequent alleles
#   2. 3rd allele freq is negligible if present -> true for 99.9% of rs-pop combos
#   extra: ref == major (no loss of information at all) -> true for 99% of rs-pop combos
#
ref_diff_major_checks = {}
assumption1_checks = {}
assumption2_checks = {}

for ens_id in vep_coding.keys():

    #debug
    # ens_id = 'ENSG00000000971'

    print(ens_id)
    gene_vep_coding = vep_coding[ens_id]
    mafs_fp = '{0}/{1}'.format(mafs_dir, ens_id)
    # output of MAF REST contains alleles for rs for which at least one allele is coding in VEP
    gene_coding_rs_mafs = pd.read_csv(mafs_fp, sep='\t')

    gene_assumption1_checks = {}
    gene_assumption2_checks = {}
    gene_ref_diff_major_checks = {}

    for rs in gene_vep_coding.uploaded_variation:
        rs_mafs = gene_coding_rs_mafs.loc[gene_coding_rs_mafs.rs == rs]
        rs_vep = gene_vep_coding.loc[(gene_vep_coding.uploaded_variation == rs)]

        rs_assumption1_checks = {}
        rs_assumption2_checks = {}
        rs_ref_diff_major_checks = {}

        for pop in rs_mafs.population:
            rs_pop_mafs = rs_mafs.loc[rs_mafs.population == pop]
            rs_pop_mafs = rs_pop_mafs.sort_values(axis=0, by='pop_allele_freq', ascending=False)
            rs_pop_alleles_sorted = rs_pop_mafs.pop_allele.tolist()
            rs_pop_allele_freqs_sorted = rs_pop_mafs.pop_allele_freq.tolist()

            # assumption 1 check
            rs_assumption1_checks[pop] = rs_vep.ref_allele.unique() in rs_pop_alleles_sorted[:2]

            # assumption 2 check
            if len(rs_pop_alleles_sorted) < 3:
                rs_assumption2_checks[pop] = True
            elif float(rs_pop_allele_freqs_sorted[2]) < 0.01:
                rs_assumption2_checks[pop] = True
            else:
                rs_assumption2_checks[pop] = False

            # check whether ref allele (VEP) is different from major allele (REST)
            rs_ref_diff_major_checks[pop] = rs_vep.ref_allele.unique() in rs_pop_alleles_sorted[:1]

        gene_assumption1_checks[rs] = pd.Series(rs_assumption1_checks)
        gene_assumption2_checks[rs] = pd.Series(rs_assumption2_checks)
        gene_ref_diff_major_checks[rs] = pd.Series(rs_ref_diff_major_checks)

    assumption1_checks[ens_id] = pd.concat(gene_assumption1_checks)
    assumption2_checks[ens_id] = pd.concat(gene_assumption2_checks)
    ref_diff_major_checks[ens_id] = pd.concat(gene_ref_diff_major_checks)

assumption1_checks = pd.concat(assumption1_checks)
assumption2_checks = pd.concat(assumption2_checks)
ref_diff_major_checks = pd.concat(ref_diff_major_checks)

assumption1_checks.index.set_names(['ens_id', 'rs', 'pop'], inplace=True)
assumption2_checks.index.set_names(['ens_id', 'rs', 'pop'], inplace=True)
ref_diff_major_checks.index.set_names(['ens_id', 'rs', 'pop'], inplace=True)

# save outputs, note implicit conversion from Series (multi-index) to dataframe
assumption1_checks.to_csv('../data/assumption1_checks.tsv', sep='\t', index=True, header=True)
assumption2_checks.to_csv('../data/assumption2_checks.tsv', sep='\t', index=True, header=True)
ref_diff_major_checks.to_csv('../data/ref_diff_major_checks.tsv', sep='\t', index=True, header=True)

assumption1_checks = pd.read_csv('../data/assumption1_checks.tsv', sep='\t')
assumption2_checks = pd.read_csv('../data/assumption2_checks.tsv', sep='\t')
ref_diff_major_checks = pd.read_csv('../data/ref_diff_major_checks.tsv', sep='\t')

# investigate checks
fails1 = assumption1_checks.loc[(assumption1_checks.iloc[:,3] == False)]
fails1.query("pop=='1000GENOMES:phase_3:EUR'").head()
fails1.shape

fails2 = assumption2_checks.loc[(assumption2_checks.iloc[:,3] == False)]
fails2.query("pop=='1000GENOMES:phase_3:EUR'").head()
fails2.shape

ref_diff_major_fails = ref_diff_major_checks.loc[(ref_diff_major_checks.iloc[:,3] == False)]
ref_diff_major_fails.query("pop=='1000GENOMES:phase_3:EUR'").head()
ref_diff_major_fails.shape

# conclusion issues 1&2 not encountered for 1000Genomes EUR
# some cases where ref != major allele, but these are dealt with since in
# these cases ref is always the second most common allele in 1000Genomes EUR

# Now perform the back-merge and ld-expand ---------------------------

ld_partners = {}
vep_coding_mafs = {}

for ens_id in vep_coding.keys():

# ens_id_ld_avail = list(list(ld_partners.index.levels)[0])
# for ens_id in (set(list(vep_coding.keys()))-set(ens_id_ld_avail)):

    # debug
    # ens_id = 'ENSG00000100342'
    print(ens_id)

    gene_vep_coding = vep_coding[ens_id]
    mafs_fp = '{0}/{1}'.format(mafs_dir, ens_id)

    # output of MAF REST contains alleles for rs for which at least one allele is coding in VEP
    gene_coding_rs_mafs = pd.read_csv(mafs_fp, sep='\t')

    # merge back in order to retain only those alleles which result in coding change

    #   1. annotate major_allele column in mafs, calculate per rs and population
    #   2. conditional join
    #   2a. merge on rs only
    #   2b. filter on (mafs.major_allele == vep.ref_allele and mafs.pop_allele == vep.allele)
    #                   OR (mafs.major_allele == vep.allele and mafs.pop_allele == vep.ref_allele)

    # 1. no group_by %>% mutate in pandas -> summarise for loop + merge back into gene_coding_rs_mafs

    rs_list, pop_list, major_allele_list = [], [], []
    for rs in gene_coding_rs_mafs.rs.unique():
        rs_mafs = gene_coding_rs_mafs.loc[gene_coding_rs_mafs.rs == rs]
        for pop in rs_mafs.population.unique():
            rs_pop_mafs = rs_mafs.loc[rs_mafs.population == pop]
            rs_pop_mafs = rs_pop_mafs.sort_values(axis=0, by='pop_allele_freq', ascending=False)
            rs_pop_alleles_sorted = rs_pop_mafs.pop_allele.tolist()
            major_allele = rs_pop_alleles_sorted[0]

            rs_list.append(rs)
            pop_list.append(pop)
            major_allele_list.append(major_allele)

    major_allele_df = pd.DataFrame({'rs': rs_list,
                                    'population': pop_list,
                                    'major_allele': major_allele_list})

    # 2. the merge

    gene_coding_rs_mafs = gene_coding_rs_mafs.merge(major_allele_df)


    full_merged = gene_vep_coding_mafs = (gene_coding_rs_mafs
                                          .merge(gene_vep_coding,
                                                 left_on=['rs'],
                                                 right_on=['uploaded_variation']))
    merged_filtered = full_merged.loc[(((full_merged['major_allele']==full_merged['ref_allele']) & (full_merged['pop_allele']== full_merged['allele'])) |
                     ((full_merged['major_allele'] == full_merged['allele']) & (full_merged['pop_allele'] == full_merged['ref_allele'])))]

    gene_vep_coding_mafs = (merged_filtered
                            .drop(columns=['uploaded_variation', 'minor_allele',
                                           'ancestral_allele', 'allele', 'ref_allele'])
                            .rename(index=str, columns={'pop_allele_freq': 'second_allele_freq',
                                                        'pop_allele': 'second_allele',
                                                        'allele_string': 'ref_allele_string'}))
    vep_coding_mafs[ens_id] = gene_vep_coding_mafs

    # biallelic rs ids, for LD expansion
    # must be biallelic because LD REST endpoint does not return allele information for the LD partner,
    #   nor allows to filter on minor allele of the input
    # 'biallelic' operationalised as: over all (allele-population_frequency) pairs,
    #   all pairs with pop_freq > 0.001 contain one of a set of maximum two alleles
    # must calculate based on mafs before back-merger with vep_coding, because also (esp.) non-coding alt alleles
    #   could cause false positives, when the non-coding alt allele has higher MAF than the coding alt allele
    # note: this is a rather stringent criterion, per-population max biallelic could be sufficient
    # even more stringent: no more than two alleles in allele string, irrespective of corresponding allele freqs

    # get number of alleles
    rs_n_alleles = {}
    for rs in gene_coding_rs_mafs['rs'].unique():
        #debug
        # rs = 'rs1285487070'
        # print(rs)
        rs_mafs = gene_coding_rs_mafs.loc[gene_coding_rs_mafs['rs'] == rs]
        # rs_n_alleles[rs] = len(set([rs_mafs['pop_allele'].iloc[i] for i in range(rs_mafs.shape[0]) if rs_mafs['pop_allele_freq'].iloc[i] > 0.0001]))
        rs_n_alleles[rs] = list(rs_mafs.allele_string.unique())[0].count('/') + 1

    # triallelic = {rs: n_alleles for rs, n_alleles in rs_n_alleles.items() if n_alleles > 2}
    rs_maxbi = {rs: n_alleles for rs, n_alleles in rs_n_alleles.items() if n_alleles < 3}


    # filter for max_biallelic, and population==EUR
    gene_vep_coding_mafs_maxbi_eur = gene_vep_coding_mafs.loc[(
            (gene_vep_coding_mafs['rs'].isin(rs_maxbi.keys())) &
            (gene_vep_coding_mafs['population'] == '1000GENOMES:phase_3:EUR'))]

    # LD-expand with ref 1KG EUR, for increased sensitivity in GWAS and eQTL lookups
    # for some genes, there may be no coding EUR rs allele in 1KG (because very low MAF)
    if gene_vep_coding_mafs_maxbi_eur.shape[0] >0:
        print("\n\nstarting ld_expand for {0} at {1}".format(ens_id, datetime.now().isoformat()),flush=True)
        ld_partners[ens_id] = run_ld_expand(rs_ids=gene_vep_coding_mafs_maxbi_eur['rs'].unique(),
                                            r2_cutoff=0.95,
                                            window_size=500)

# concat over ens_ids
ld_partners = pd.concat(ld_partners)
vep_coding_mafs = pd.concat(vep_coding_mafs)

# export for analysis in R
ld_partners.to_csv('../data/ld_partners', sep='\t',index=True)
vep_coding_mafs.to_csv('../data/vep_coding_mafs', sep='\t',index=True)



