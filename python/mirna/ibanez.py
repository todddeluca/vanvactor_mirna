

import csv
import os

import rdflib

import config
from mirna.core import (db_pred, dro_iri, homo_iri, homologous_to_pred,
                        mirbase_id_db, mirbase_id_iri, organism_pred)


# The two methods used to infer homology in the 2008 Ibanez-Ventoso paper.
FIVE_PRIME = 'five_prime'
FULL_SEQ = 'full_sequence'

###############################
# IBANEZ FLY-HUMAN MIR HOMOLOGS


def ibanez_five_prime_mir_homologs_rdf_path():
    return os.path.join(config.datadir, 'ibanez', '2008', 'ibanez-2008-five_prime-mir-homologs.nt')


def ibanez_fullseq_mir_homologs_rdf_path():
    return os.path.join(config.datadir, 'ibanez', '2008', 'ibanez-2008-fullseq-mir-homologs.nt')


def write_ibanez_five_prime_mir_homologs_rdf():
    print 'write_ibanez_five_prime_mir_homologs_rdf'
    write_method_mir_homologs_rdf(FIVE_PRIME, ibanez_five_prime_mir_homologs_rdf_path())


def write_ibanez_fullseq_mir_homologs_rdf():
    print 'write_ibanez_fullseq_mir_homologs_rdf'
    write_method_mir_homologs_rdf(FULL_SEQ, ibanez_fullseq_mir_homologs_rdf_path())


def write_method_mir_homologs_rdf(method, filename):
    '''
    Write the homologs for a specific method to a file as n-triples rdf.
    method: FIVE_PRIME or FULL_SEQ
    '''
    print 'Populating graph for Ibanez 2008 miR fly-human homologs for method {}'.format(method)
    graph = rdflib.Graph()
    fly_mirs = set()
    human_mirs = set()

    for fly_mir_id, human_mir_id, meth in gen_fly_human_homologs():
        if meth != method:
            continue

        fly_mir = mirbase_id_iri(fly_mir_id)
        human_mir = mirbase_id_iri(human_mir_id)
        graph.add((fly_mir, homologous_to_pred, human_mir))

        # Examine whether there is a many-to-many mapping of these mirs
        if fly_mir in fly_mirs:
            print 'fly mir appears twice', fly_mir_id, human_mir_id, method
        if human_mir in human_mirs:
            print 'human mir appears twice', fly_mir_id, human_mir_id, method

        fly_mirs.add(fly_mir)
        human_mirs.add(human_mir)

    for mir in fly_mirs:
        graph.add((mir, db_pred, mirbase_id_db))
        graph.add((mir, organism_pred, dro_iri))

    for mir in human_mirs:
        graph.add((mir, db_pred, mirbase_id_db))
        graph.add((mir, organism_pred, homo_iri))

    # serialize as n-triples
    with open(filename, 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))


def download_fly_human_homologs():
    raise NotImplementedError('''
Since there is no data file on the internet available for download, the Ibanez
data was hand-copied from the paper into a CSV file.  Upload it from my laptop
to the server location $datadir/ibanez/2008/2008_Ibanez-Ventoso_drosophila_melanogaster_homo_sapiens_mir_orthologs.csv
''')


def gen_fly_human_homologs():
    '''
    From PLoS ONE. 2008; 3(7): e2818.
    Published online 2008 July 30. doi:  10.1371/journal.pone.0002818
    PMCID: PMC2486268
    Sequence Relationships among C. elegans, D. melanogaster and Human microRNAs Highlight the Extensive Conservation of microRNAs in Biology
    Carolina Ibanez-Ventoso, Mehul Vora, and Monica Driscoll

    Table 5: fly-human homologs, using 5' and 70% methods.  Note that the CSV
    version of this table is overly complex b/c it replicates the paper version
    which conflates results for both methods instead of using two separate 
    tables.  This results in empty human columns for the different methods and
    empty drosophila columns.  Saves on ink and paper.

    Generate tuples of drosophila mir id, homo mir id, method) for the ">= 70%
    full sequence homology" and "5-prime sequence homology" methods for
    asserting homology between fly and human mirs from the paper.
    '''
    fn = os.path.join(config.datadir, 'ibanez', '2008',
                      '2008_Ibanez-Ventoso_drosophila_melanogaster_homo_sapiens_mir_orthologs.csv')
    with open(fn) as fh:
        reader = csv.reader(fh)
        current_dro_mir = None
        for i, row in enumerate(reader):
            # row 0 is a header row
            if i < 1:
                continue

            if len(row) != 4:
                raise Exception('Unexpected row length', i, len(row), row)

            mir_group, dro_mir, five_prime_homo_mir, full_homo_mir = row
            if dro_mir and dro_mir != current_dro_mir:
                current_dro_mir = dro_mir

            if not five_prime_homo_mir and not full_homo_mir:
                raise Exception('No human homolog defined for row!', i, row)
            if not current_dro_mir:
                raise Exception('No drosophila homolog defined for row!', i, row)
            if five_prime_homo_mir:
                yield current_dro_mir, five_prime_homo_mir, FIVE_PRIME
            if full_homo_mir:
                yield current_dro_mir, full_homo_mir, FULL_SEQ


def print_fly_human_homologs():
    '''
    See what the generated tuples of human-fly mir homologs look like to make
    sure they do not look wrong.
    '''
    for dro_mir, homo_mir, method in gen_fly_human_homologs():
        print method, dro_mir, homo_mir



