
import os
import csv

import rdflib

import config
from mirna.core import (affymetrix_probeset_db, affymetrix_probeset_iri,
                        db_pred, dro_iri, mirbase_id_db, mirbase_id_iri,
                        organism_pred, regulates_pred)


# Tissue type used in McNeill CSV file
MUSCLE = 'muscle'
CNS = 'CNS'


#################################
# VAN VACTOR / MCNEILL EXPERIMENT
# ELIZABETH MCNEILL SCREEN

def download_mcneill_screen():
    raise NotImplementedError()


def gen_mcneill_screen_mir_targets():
    '''
    Generate a tuple for every row in the CSV file containing the experimental
    results for affymetrix probesets that were differentially expressed in
    drosophila 'muscle' and 'CNS' tissues when 7 miRs were (individually)
    perturbed.
    '''
    fn = os.path.join(config.datadir, 'mcneill', '20130523',
                      '20130523_vanvactor_fly_tissue_7_mir_targets.csv')
    with open(fn) as fh:
        reader = csv.reader(fh)
        for i, row in enumerate(reader):
            # row 0 is a header row
            if i < 1:
                continue

            if len(row) != 3:
                raise Exception('Unexpected row length', i, len(row), row)

            mir_id, tissue, probeset_id = row
            yield mir_id, tissue, probeset_id


def mcneill_muscle_screen_rdf_path():
    return os.path.join(config.datadir, 'mcneill', '20130523', 'mcneill_muscle_screen.nt')


def mcneill_cns_screen_rdf_path():
    return os.path.join(config.datadir, 'mcneill', '20130523', 'mcneill_cns_screen.nt')


def write_mcneill_muscle_screen_rdf():
    '''
    Write an rdf file for muscle tissue results from Elizabeth McNeill's
    expression array screen of 7 fly miRs.
    '''
    print 'write_mcneill_muscle_screen_rdf'
    write_mcneill_tissue_screen_rdf(MUSCLE, mcneill_muscle_screen_rdf_path())


def write_mcneill_cns_screen_rdf():
    '''
    Write an rdf file for CNS tissue results from Elizabeth McNeill's
    expression array screen of 7 fly miRs.
    '''
    print 'write_mcneill_cns_screen_rdf'
    write_mcneill_tissue_screen_rdf(CNS, mcneill_cns_screen_rdf_path())


def write_mcneill_tissue_screen_rdf(tissue, filename):
    '''
    Example rdf triples showing the regulates predicate and annotating each mir and
    probeset with database and taxon:

        <http://purl.mirbase.org/mirna_id/dme-miR-34> <http://purl.example.com/owl/regulates> <http://purl.affymetrix.com/probeset/1626730_s_at> .
        <http://purl.affymetrix.com/probeset/1626730_s_at> <http://purl.uniprot.org/core/organism> <http://purl.uniprot.org/taxonomy/7227> .
        <http://purl.affymetrix.com/probeset/1626730_s_at> <http://purl.uniprot.org/core/database> <http://purl.example.com/database/affymetrix_probeset> .
        <http://purl.mirbase.org/mirna_id/dme-miR-34> <http://purl.uniprot.org/core/database> <http://purl.example.com/database/mirbase_id> .
        <http://purl.mirbase.org/mirna_id/dme-miR-34> <http://purl.uniprot.org/core/organism> <http://purl.uniprot.org/taxonomy/7227> .
    '''
    graph = rdflib.Graph()
    mirs = set()
    probesets = set()

    for mir_id, tissue_type, probeset_id in gen_mcneill_screen_mir_targets():
        if tissue_type == tissue:
            mir = mirbase_id_iri(mir_id)
            probeset = affymetrix_probeset_iri(probeset_id)
            graph.add((mir, regulates_pred, probeset))
            mirs.add(mir)
            probesets.add(probeset)

    for mir in mirs:
        graph.add((mir, db_pred, mirbase_id_db))
        graph.add((mir, organism_pred, dro_iri))

    for probeset in probesets:
        graph.add((probeset, db_pred, affymetrix_probeset_db))
        graph.add((probeset, organism_pred, dro_iri))

    with open(filename, 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))



