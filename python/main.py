
import argparse
import logging
import os
import sys
import collections

import pandas
import numpy as np

import config
import flybaseutil



def example(message):

    print 'Example CLI function.'
    print 'message:', message


def flybase_genes(gene_conversion_table, outfile):
    '''
    Return a list of unique flybase gene ids from the gene_conversion table,
    skipping ids that did not convert.

    outfile: a file path where the genes will be written.  If None, genes will
    be written to standard output.


    '''
    genes = flybaseutil.select_flybase_gene_ids(gene_conversion_table)
    with open(outfile, 'w') as fh:
        for gene in genes:
            fh.write('{}\n'.format(gene))


def synaptomedb_v1_all_genes_file():
    return os.path.join(config.datadir, 'synaptomedb', 'v1.06',
                        'synaptomedb-1.06-all_genes.csv')


def microcosm_v5_human_targets_file():
    return os.path.join(config.microcosm_dir, 'v5',
                        'microcosm-v5-homo_sapiens_targets.txt')


def parse_synaptomedb_all_genes(filename):
    genes = []
    df = pandas.read_csv(filename)
    for i, row in df.iterrows():
        gid = row['gid']
        gene_symbol = row['gene_symbol']
        gene_info = row['gene_info']
        ensembl_gene_ids = row['ensembl_gene_id'].split('; ') if row['ensembl_gene_id'] is not np.nan else []
        ensembl_trans_ids = row['ensembl_trans_id'].split('; ') if row['ensembl_trans_id'] is not np.nan else []
        genes.append({'gid': gid, 'gene_symbol': gene_symbol,  'gene_info': gene_info,  'ensembl_gene_ids': ensembl_gene_ids, 'ensembl_trans_ids': ensembl_trans_ids})

    # with open(filename) as fh:
        # reader =csv.reader
        # for line in fh:
            # # skip comments and blank lines
            # lstripped = line.lstrip()
            # if lstripped.startswith('#') or not lstripped:
                # continue
            # fields

    return genes


def parse_microcosm_human_targets(filename):
    targets = []
    with open(filename) as fh:
        for line in fh:
            # skip comments and blank lines
            lstripped = line.lstrip()
            if lstripped.startswith('#') or not lstripped:
                continue

            fields = line.rstrip('\n').split('\t')
            assert len(fields) == 13
            mir = fields[1]
            ensembl_trans_id = fields[11]
            gene_symbol = fields[12]
            targets.append({'mir': mir, 'gene_symbol': gene_symbol, 'ensembl_trans_id': ensembl_trans_id})

    return targets


def map_human_synaptic_genes_to_human_mirs():
    '''
    Print out a report on the mirs that target each synaptic gene.
    Specifically, print each gene and a list of mirs that target it, if any.
    Results in ~/deploy/vanvactor_mirna/data/20130315_human_miRs_that_target_synaptic_genes.txt
    '''
    # Read and parse synaptomedb file, so we can find out what human synapse
    # genes we have, and their associated ensembl transcript ids.
    genes = parse_synaptomedb_all_genes(synaptomedb_v1_all_genes_file())
    # print 'genes', genes

    # Read and parse microcosm human gene targets file, so we can find out
    # what ensembl transcripts are targeted by what miRs.
    targets = parse_microcosm_human_targets(microcosm_v5_human_targets_file())
    # print 'targets', targets

    # quickly look up the mirs that target a transcript
    trans2mirs = collections.defaultdict(set)
    for target in targets:
        trans2mirs[target['ensembl_trans_id']].add(target['mir'])

    gene2mirs = collections.defaultdict(set)
    mir2genes = collections.defaultdict(set)
    allmirs = set()
    for gene in genes:
        print 'Synaptic gene {}, {}'.format(gene['gene_symbol'], gene['gene_info'])
        if not gene['ensembl_trans_ids']:
            print 'This gene is not targeted by any miRs.'
            print 'Reason: the gene is not annotated with any Ensembl transcript ids.'
            print 'End'
            print
            continue
        mirs = set(mir for trans in gene['ensembl_trans_ids'] for mir in trans2mirs[trans] if mir.startswith('hsa'))
        if not mirs:
            print 'This gene is not targeted by any miRs.'
            print 'Reason: the gene transcripts are not targeted by any miRs.'
            print 'End'
            print
            continue
        for mir in mirs:
            mir2genes[mir].add(gene['gene_symbol'])
            gene2mirs[gene['gene_symbol']].add(mir)
            allmirs.add(mir)

        print 'Mirs that target this gene:', ', '.join(sorted(mirs))
        print 'End'
        print
        # print 'allmirs', allmirs
        # print 'len(allmirs)', len(allmirs)

    # print 'len(allmirs', len(set(t['mir'] for t in targets if t['mir'].startswith('hsa')))

    # for item in sorted((len(mir2genes[m]), m) for m in mir2genes):
        # print item


def main():

    parser = argparse.ArgumentParser(description='')
    subparsers = parser.add_subparsers(dest='action', help='')

    # parser name corresponds to function name
    subparser = subparsers.add_parser('example', help='Example help text')
    subparser.add_argument('message', help='print this')

    # parser name corresponds to function name
    subparser = subparsers.add_parser('flybase_genes', help='Extract flybase')
    subparser.add_argument('gene_conversion_table')
    subparser.add_argument('outfile')

    subparser = subparsers.add_parser('map_human_synaptic_genes_to_human_mirs', help='')

    # invoke a function named by action whose keyword parameters correspond to
    # cli arguments
    args = parser.parse_args()
    kws = dict(vars(args))
    del kws['action']
    return globals()[args.action](**kws)


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception('')
        raise


