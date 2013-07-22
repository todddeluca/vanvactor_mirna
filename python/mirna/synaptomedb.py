

import collections
import csv
import os

import rdflib

import config

from mirna.core import (classified_pred, db_pred, ensembl_gene_db,
                        ensembl_gene_iri, homo_iri, organism_pred,
                        synaptic_iri)

#######################
# SYNAPTOMEDB FUNCTIONS


def download_synaptome_v1():
    raise NotImplementedError()


def synaptomedb_v1_all_genes_file():
    return os.path.join(config.datadir, 'synaptomedb', 'v1.06',
                        'synaptomedb-1.06-all_genes.csv')


def synaptomedb_v1_rdf_path():
    return os.path.join(config.datadir, 'synaptomedb', 'v1.06',
                        'synaptomedb-1.06.nt')


def write_synaptomedb_rdf():
    print 'write_synaptomedb_rdf'
    graph = rdflib.Graph()

    # describe the IRIs as human, synaptic, ensembl gene
    for ensembl_gene_id in gen_synaptomedb_ensembl_genes():
        gene = ensembl_gene_iri(ensembl_gene_id)
        graph.add((gene, classified_pred, synaptic_iri))
        graph.add((gene, db_pred, ensembl_gene_db))
        graph.add((gene, organism_pred, homo_iri))

    with open(synaptomedb_v1_rdf_path(), 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))


def gen_synaptomedb_ensembl_genes():

    genes = parse_synaptomedb_all_genes()
    ensembl_gene_ids = set([ensembl_gene_id for gene in genes for
                            ensembl_gene_id in gene['ensembl_gene_ids']])
    for g in sorted(ensembl_gene_ids):
        yield g


def parse_synaptomedb_all_genes(filename=None):
    '''
    Return a list of genes.  Each gene is a dict containing 'gid' (a synaptome
    gene id?), 'gene_symbol', 'gene_info' (a description of the gene),
    'ensembl_gene_ids' (a list of ensembl gene ids), and 'ensembl_trans_ids' (a
    list of ensembl transcript ids).
    '''
    if filename is None:
        filename = synaptomedb_v1_all_genes_file()

    # These are the fields of the csv file
    headers = ["gid", "gene_symbol", "gene_synonyms", "gene_info", "chr",
               "start", "end", "hgnc", "pdb", "ipi", "unigene", "bp", "mf",
               "cc", "ensembl_gene_id", "ensembl_trans_id", "ensembl_pro_id",
               "interpro_id", "pfam_id", "pir_id", "pirsf_id", "ccds_id",
               "nuc_acc", "pro_acc", "nuc_gi", "pro_gi", "homo_id",
               "homo_gene_id", "homo_ens_gene_id", "homo_pro_id", "pubmed",
               "omim", "path_kegg", "path_biocarta", "path_reactome",
               "path_ecnumber", "path_gsea", "path_hprd", "ref_gen_acc",
               "ref_mrna_acc", "ref_pro_acc", "ref_gen_gi", "ref_nuc_gi",
               "ref_pro_gi"]
    Fields = collections.namedtuple('Fields', headers)

    genes = []
    with open(filename) as fh:
        reader = csv.reader(fh)
        for i, row in enumerate(reader):
            # row 0 is a header row
            if i < 1:
                continue
            row = Fields(*row) # look up values using named attributes
            gid = row.gid
            gene_symbol = row.gene_symbol
            gene_info = row.gene_info
            ensembl_gene_ids = row.ensembl_gene_id.split('; ') if row.ensembl_gene_id else []
            ensembl_trans_ids = row.ensembl_trans_id.split('; ') if row.ensembl_trans_id else []
            genes.append({'gid': gid, 'gene_symbol': gene_symbol,  'gene_info':
                          gene_info,  'ensembl_gene_ids': ensembl_gene_ids,
                          'ensembl_trans_ids': ensembl_trans_ids})
    return genes


