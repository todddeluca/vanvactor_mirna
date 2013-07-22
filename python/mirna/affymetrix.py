
import collections
import csv
import os
import re

import rdflib

import config
from mirna.core import (affymetrix_probeset_db, affymetrix_probeset_iri,
                        db_pred, dro_iri, flybase_gene_db, flybase_gene_iri,
                        organism_pred, see_also_pred)



############
# AFFYMETRIX


def affymetrix_fly_annotations_rdf_path():
    return os.path.join(config.datadir, 'affymetrix', 'na33', 'Drosophila_2.na33.annot.nt')


def write_affymetrix_fly_annotations_rdf():
    print 'write_affymetrix_fly_annotations_rdf'
    graph = rdflib.Graph()
    fly_genes = set()
    for probeset_id, gene_id in gen_affymetrix_fly_probe_to_flybase_mapping():
        gene = flybase_gene_iri(gene_id)
        probeset = affymetrix_probeset_iri(probeset_id)
        # Should this be see_also or something indicating that the probeset
        # binds to the gene transcript.
        graph.add((probeset, see_also_pred, gene))
        graph.add((probeset, db_pred, affymetrix_probeset_db))
        graph.add((probeset, organism_pred, dro_iri))
        fly_genes.add(gene)

    # since genes might appear more than once in the mapping data, collect
    # them into a set and only make triples for them once.
    for gene in fly_genes:
        graph.add((gene, db_pred, flybase_gene_db))
        graph.add((gene, organism_pred, dro_iri))

    with open(affymetrix_fly_annotations_rdf_path(), 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))


def download_affymetrix_fly_annotations_file():
    raise NotImplementedError('''
Go to page http://www.affymetrix.com/support/technical/byproduct.affx?product=fly-20
and download:

    Current NetAffx Annotation Files
    Drosophila_2 Annotations, CSV format, Release 33 (8.4 MB, 10/30/12)
    This will download 'Drosophila_2.na33.annot.csv.zip'.
    Unzip this file, then cd to Drosophila_2.na33.annot.csv/.
    Then copy Drosophila_2.na33.annot.csv to the affymetrix/na33/ data dir.
    ''')


def affymetrix_fly_annotations_file():
    return os.path.join(config.datadir, 'affymetrix', 'na33', 'Drosophila_2.na33.annot.csv')


def gen_affymetrix_fly_probe_to_flybase_mapping():
    '''
    Parse the affymetrix fly annotation file and yield
    tuples of (probeset_id, flybase_gene_id) for each probeset that maps to
    one or more flybase genes.
    Example probeset id: '1641386_at'
    Example flybase gene ids: 'FBgn0010905', 'FBgn0260480'

    Some probeset ids map to no flybase ids.  Some map to multiple ids.  Most
    map to one id.
    '''
    with open(affymetrix_fly_annotations_file()) as fh:
        reader = csv.reader(fh)
        for i, row in enumerate(reader):
            # rows 0 to 18 are comments.
            # row 19 contains column headers.
            if i < 20:
                continue

            if len(row) != 41:
                raise Exception('Unexpected row length', i, len(row), row)

            probeset_id = row[0]
            genes_str = row[24]
            if genes_str == '---':
                genes = []
            else:
                genes = genes_str.split(' /// ')
            for gene in genes:
                if not re.search(r'^FBgn\d+$', gene):
                    raise Exception('Bad flybase gene', probeset_id, gene,
                                    genes_str, row)
                yield probeset_id, gene


def map_affymetrix_fly_probe_ids_to_flybase_genes():
    '''
    Parse the affymetrix fly annotation file and map probeset ids (like
    '1641386_at') to lists of 0 or more flybase gene ids (like ['FBgn0010905',
    'FBgn0260480']).

    Return a dict mapping each probeset id to its possibly empty list of
    flybase genes.
    '''

    probe2genes = collections.defaultdict(list)
    for probe, gene in gen_affymetrix_fly_probe_to_flybase_mapping():
        probe2genes[probe].append(gene)

    for p in probe2genes:
        print p, probe2genes[p]

    return probe2genes


