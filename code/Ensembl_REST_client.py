#!/usr/bin/env python

import urllib
# import json
import time
import pandas as pd
import csv
import sys
import requests
import numpy as np
from datetime import datetime
import random
import math

# For reference and to implement new methods:
# see API documentation on https://github.com/Ensembl/ensembl-rest/wiki
# see list of endpoints on: http://rest.ensembl.org


class EnsemblRestClient(object):
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, method='GET', hdrs=None, params=None, http_error_action=None, data=None,
                            jsondata=None, retry_count=0, retry_max=20, sleep=True):

        if sleep:
            random_sleep = round(min(10, np.random.exponential(2)), ndigits=5)
            print("random sleep for {} sec..".format(random_sleep), flush=True)
            time.sleep(random_sleep)

        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if http_error_action is None:
            http_error_action = {'429': 'retry'}
        elif isinstance(http_error_action, dict):
            if '429' not in http_error_action.keys():
                http_error_action['429'] = 'retry'
        else:
            sys.exit('HTTPError_action should be of type dict.')

        # resp_data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            response = requests.request(method=method, url=str(self.server) + endpoint, params=params,
                                        headers=hdrs, data=data, json=jsondata)

            # required to catch the exception (in contrast to urllib requests which raise the exception automatically)
            response.raise_for_status()

            if response is not None:
                resp_data = response.json()
            else:
                print("got response None")
                resp_data = self.perform_rest_action(endpoint=endpoint, hdrs=hdrs, params=params, method=method,
                                                     http_error_action=http_error_action, data=data,
                                                     jsondata=jsondata, retry_count=retry_count, retry_max=retry_max)

            self.req_count += 1



        # https://stackoverflow.com/questions/16511337/correct-way-to-try-except-using-python-requests-module#16511493
        except requests.exceptions.HTTPError as e:
            status_code = str(e.response.status_code)
            action = http_error_action.get(status_code)
            if action == 'retry':
                print("got error {0} with headers {1}. Retrying...".format(status_code, response.headers), flush=True)
                retry_count += 1

                if retry_count > retry_max:
                    resp_data = "retry_max reached"

                # check if we are being rate limited by the server
                elif 'Retry-After' in e.response.headers.keys():
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    resp_data = self.perform_rest_action(endpoint=endpoint, hdrs=hdrs, params=params, method=method,
                                                         http_error_action=http_error_action, data=data,
                                                         jsondata=jsondata, retry_count=retry_count,
                                                         retry_max=retry_max)
                else:
                    resp_data = self.perform_rest_action(endpoint=endpoint, hdrs=hdrs, params=params, method=method,
                                                         http_error_action=http_error_action, data=data,
                                                         jsondata=jsondata, retry_count=retry_count,
                                                         retry_max=retry_max)
            elif action == 'pass':
                print('Ignoring error {0} on endpoint {1}'.format(status_code, endpoint))
                resp_data = e
            else:
                sys.stderr.write(
                    'Request failed with unhandled status code for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(
                        endpoint, e))
                resp_data = e

        return resp_data

    def get_gene_eqtls(self, species, ens_gene_id, p_threshold, http_error_action):

        df = None

        resp = self.perform_rest_action(
            endpoint='/eqtl/id/{0}/{1}'.format(species, ens_gene_id),
            http_error_action=http_error_action
        )

        if isinstance(resp, urllib.error.HTTPError):
            return resp
        elif resp is not None:
            df = pd.DataFrame(resp).query('value < {0}'.format(p_threshold))

        return df

    def get_variants(self, species, symbol):
        genes = self.perform_rest_action(
            endpoint='/xrefs/symbol/{0}/{1}'.format(species, symbol),
            params={'object_type': 'gene'}
        )
        if genes:
            stable_id = genes[0]['id']
            variants = self.perform_rest_action(
                '/overlap/id/{0}'.format(stable_id),
                params={'feature': 'variation'}
            )
            return variants
        return None

    def get_variants_by_id(self, ens_id, feature):
        variants = self.perform_rest_action(
            '/overlap/id/{0}'.format(ens_id),
            params={'feature': feature}
        )
        return variants

    def get_mafs_post(self, rsids, params=None):

        if params is None:
            params = {'pops': 1}

        mafs = self.perform_rest_action(
            endpoint='/variation/homo_sapiens',
            jsondata={'ids': rsids},
            method='POST',
            params=params,
            http_error_action={'504': 'retry',
                               '503': 'retry',
                               '500': 'retry',
                               '507': 'retry',
                               '502': 'retry'}
        )

        return mafs

    def get_vep_post(self, species, rsids, gene_id=None, canonical=True, max_tsl=5, params=None):

        canon_chroms = [str(x) for x in range(1, 23)] + ['X', 'Y', 'MT']

        max_rs = 200
        if len(rsids) > max_rs:
            sys.exit('max {} rs IDs per request'.format(max_rs))
        jsondata = {'ids': rsids}

        if params is None:
            params = {'tsl': 1}
        elif 'tsl' not in params.keys():
            params['tsl'] = 1

        veps = self.perform_rest_action(
            endpoint='/vep/{0}/id'.format(species),
            jsondata=jsondata,
            method='POST',
            params=params,
            http_error_action={'504': 'retry',
                               '503': 'retry',
                               '500': 'retry',
                               '507': 'retry',
                               '502': 'retry'})

        if veps == "retry_max reached":
            filtered = veps
        else:
            filtered = []
            for variant_dict in veps:
                if (((variant_dict['seq_region_name'] in canon_chroms) or (canonical is False))
                        and 'transcript_consequences' in variant_dict.keys()):
                    tmp = variant_dict
                    tx_i = [i
                            for i in range(len(variant_dict['transcript_consequences']))
                            if ((variant_dict['transcript_consequences'][i]['gene_id'] == gene_id
                                 or gene_id is None)
                                and
                                (int(variant_dict['transcript_consequences'][i].get('tsl', 9999)) <= max_tsl))
                            ]
                    tmp['transcript_consequences'] = [variant_dict['transcript_consequences'][i] for i in tx_i]
                    filtered.append(tmp)
        return filtered

    def get_ld_expand_window(self, rs_id, window_size=500, r2_cutoff=0, population='1000GENOMES:phase_3:EUR',
                             species='homo_sapiens'):
        # window size in kb

        params = {'window_size': window_size,
                  'r2': r2_cutoff}

        rs_partners = self.perform_rest_action(
            endpoint= '/ld/{0}/{1}/{2}'.format(species,rs_id, population),
            params=params,
            method='GET',
            sleep=False,
            http_error_action={'504': 'retry',
                               '503': 'retry',
                               '500': 'retry',
                               '507': 'retry',
                               '502': 'retry'}
        )

        rs_partners_df = pd.DataFrame(rs_partners)   #from list of dicts

        return rs_partners_df


def run_vep(species, rsids, gene_id, server='http://rest.ensembl.org'):
    # retrieves variant consequences in the context of a given gene

    max_ids = 200
    veps = []
    for i in range(0, len(rsids), max_ids):
        print("\n\n{0}/{1} at {2}".format(int(i / max_ids)+1, round(len(rsids) / max_ids), datetime.now().isoformat()),
              flush=True)
        query_ids = rsids[i:i + max_ids]

        query_success = False
        while query_success is False:
            client = EnsemblRestClient(server)
            veps_returned = client.get_vep_post(species=species, rsids=query_ids,
                                                gene_id=gene_id, max_tsl=2, canonical=True,
                                                params={
                                                    'tsl': 1,
                                                    'GeneSplicer': 1,
                                                    'MaxEntScan': 1,
                                                    'domains': 1,
                                                    'miRNA': 1})
            if veps_returned != "retry_max reached":
                query_success = True
            else:
                print("max retries reached, reinitializing client..", flush=True)

        veps += veps_returned

    # summarise consequences per variant
    rs, tx_id, tsl, conseq, sift, polyphen, codons, aas, variant_allele, rs_alleles = [], [], [], [], [], [], [], [], [], []
    variant_overview = {}
    for vd in veps:  # variant_dict
        rs_id = vd['input']
        allele_string = vd['allele_string']
        for td in vd['transcript_consequences']:  # transcript_dict
            for ct in td['consequence_terms']:  # consequence term
                rs.append(rs_id)
                tx_id.append(td['transcript_id'])
                tsl.append(td['tsl'])
                conseq.append(ct)
                sift.append(td.get('sift_prediction'))
                polyphen.append(td.get('polyphen_prediction'))
                codons.append(td.get('codons'))
                aas.append(td.get('amino_acids'))
                # important to include exact allele, rs-id uniquely tags location and type of variation, NOT individual allele
                variant_allele.append(td.get('variant_allele'))
                rs_alleles.append(allele_string)

    variant_overview = pd.DataFrame(index=rs,
                                    data={'tx_id': tx_id,
                                          'tsl': tsl,
                                          'conseq': conseq,
                                          'sift': sift,
                                          'polyphen': polyphen,
                                          'codons': codons,
                                          'aas': aas,
                                          'variant_allele': variant_allele,
                                          'rs_alleles': rs_alleles})
    return variant_overview


def run_variants(species, symbol, server = 'http://rest.ensembl.org'):
    client = EnsemblRestClient(server)
    variants = client.get_variants(species, symbol)
    if variants:
        for v in variants:
            print('{seq_region_name}:{start}-{end}:{strand} ==> {id} ({consequence_type})'.format(**v))


def run_variants_by_id(ens_id, server = 'http://rest.ensembl.org', feature='variation'):
    client = EnsemblRestClient(server)
    variants = client.get_variants_by_id(ens_id, feature)
    rsids = [var_dict['id'] for var_dict in variants if var_dict['feature_type'] == 'variation']
    svids = [var_dict['id'] for var_dict in variants if var_dict['feature_type'] != 'variation']

    return [rsids, svids]


def run_variant_mafs(rsids, server='http://rest.ensembl.org'):


    max_ids = 200
    rs_mafs = {}
    # for i in range(0, len(rsids), max_ids):
    for i in range(0, len(rsids), max_ids):
        print(
            "\n\n{0}/{1} at {2}".format(int(i / max_ids) + 1, math.ceil(len(rsids) / max_ids), datetime.now().isoformat()),
            flush=True)
        query_ids = rsids[i:i+max_ids]

        query_success = False
        while query_success is False:
            client = EnsemblRestClient(server)
            mafs_returned = client.get_mafs_post(rsids=query_ids)
            if mafs_returned != "retry_max reached":
                query_success = True
            else:
                print("max retries reached, reinitializing client..", flush=True)
        rs_mafs.update(mafs_returned)

    rs_mafs_list = [{'rs': rsid,
                     'var_class': rs_dict.get('var_class'),
                     'allele_string': rs_dict.get('mappings')[0].get('allele_string'),
                     'minor_allele': rs_dict.get('minor_allele'),
                     'ancestral_allele': rs_dict.get('ancestral_allele'),
                     'synonyms': rs_dict.get('synonyms'),
                     'population': pop_dict.get('population'),
                     'pop_allele': pop_dict.get('allele'),
                     'pop_allele_freq': pop_dict.get('frequency')}
                    for rsid, rs_dict in rs_mafs.items()
                    for pop_dict in rs_dict.get('populations')]
    rs_mafs_df = pd.DataFrame(rs_mafs_list)

    return rs_mafs_df


def run_eqtls(species, ens_gene_ids, p_threshold, mydir, server='http://rest.ensembl.org',
              client=None, http_error_action=None):
    import os
    if not os.path.exists(mydir):
        os.makedirs(mydir)
    if client is None:
        client = EnsemblRestClient(server)

    if isinstance(ens_gene_ids, (list, pd.Series)):

        gene_status = {}
        status_fp = mydir+'gene_status.txt'

        for ens_gene_id in ens_gene_ids:
            # if ens_gene_id in no_eqtl:
            #     print(ens_gene_id + ' reported as not present in EQTL database, so skipping.')
            #     continue
            fp = mydir+'/'+ens_gene_id

            if os.path.isfile(fp):
                print(ens_gene_id+' already retrieved, so skipped.')
                status = 'success'
            else:
                print('retrieving eqtls for ' + ens_gene_id)
                resp = client.get_gene_eqtls(species, ens_gene_id, p_threshold, http_error_action)
                if isinstance(resp, urllib.error.HTTPError):
                    status = str(resp)
                elif resp is None:
                    status = 'no content'
                else:
                    resp.to_json(fp, orient='table')
                    status = 'success'

            gene_status[ens_gene_id] = status

        with open(file=status_fp, mode='w') as status_file:
            writer = csv.writer(status_file, delimiter='\t')
            writer.writerows(gene_status.items())

    else:
        ens_gene_id = ens_gene_ids
        fp = mydir + '/' + ens_gene_id
        print('retrieving eqtls for ' + ens_gene_id)
        eqtl_dict = {ens_gene_id: client.get_gene_eqtls(species, ens_gene_id, p_threshold, http_error_action)}
        eqtl_dict.to_json(fp, orient='table')

    return None


def run_ld_expand(rs_ids, r2_cutoff, population='1000GENOMES:phase_3:EUR', window_size=500,
                  server='http://rest.ensembl.org'):

    client = EnsemblRestClient(server)

    ld_partners = {}
    for i in range(len(rs_ids)):
        rs_id = rs_ids[i]
        print("rs {0}/{1}".format(i+1, len(rs_ids)))
        ld_partners[rs_id] = client.get_ld_expand_window(rs_id, window_size=window_size, r2_cutoff=r2_cutoff,
                                                         population=population)

    return pd.concat(ld_partners)


