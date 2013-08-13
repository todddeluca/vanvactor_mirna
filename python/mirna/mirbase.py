
import os
import re

import rdflib

import config
from mirna.core import (current_id_pred, db_pred, dro_iri, has_id_pred,
                        homo_iri, mature_mirna_cl, mirbase_acc_db,
                        mirbase_acc_iri, mirbase_id_db, mirbase_id_iri,
                        organism_pred, type_pred)
from mirna.util import download, call # contains download, gunzip, call


#########
# MIRBASE


def aliases_url():
    return 'ftp://mirbase.org/pub/mirbase/19/aliases.txt.gz'


def aliases_dest():
    return os.path.join(config.datadir, 'mirbase', '19',
                        os.path.basename(aliases_url()))


def aliases_filename():
    return aliases_dest().rstrip('.gz')


def mirbase_rdf_path():
    return os.path.join(config.datadir, 'mirbase', '19', 'mirbase.nt')


def download_mirbase_aliases():

    url = aliases_url()
    dest = aliases_dest()
    download(url, dest)
    call('gunzip {}'.format(dest))


def write_mirbase_rdf():
    print 'write_rdf'
    graph = rdflib.Graph()
    human_accs = set()
    human_ids = set()
    fly_accs = set()
    fly_ids = set()

    mature_re = re.compile(r'MIMAT\d+$')
    # link mirbase accessions to mirbase ids and collect the accession and ids
    with open(aliases_filename()) as fh:
        for line in fh:
            # Example lines:
            # MI0000001       cel-let-7L;cel-let-7;
            # MI0000003       cel-mir-1;
            mirbase_acc, mirbase_ids_str = line.strip().split('\t')
            mirbase_ids = mirbase_ids_str.rstrip(';').split(';')
            acc = mirbase_acc_iri(mirbase_acc)
            current_mirbase_id = mirbase_ids[-1]
            # link the mirbase acc to all ids
            for mirbase_id in mirbase_ids:
                mid = mirbase_id_iri(mirbase_id)
                # for now only do human and fly data
                if mirbase_id.startswith('dme') or mirbase_id.startswith('hsa'):
                    # link the mirbase acc to the mirbase id
                    graph.add((acc, has_id_pred, mid))
                    if mirbase_id == current_mirbase_id:
                        # link the mirbase acc to current mirbase id
                        graph.add((acc, current_id_pred, mid))

                    if mirbase_id.startswith('dme'):
                        fly_accs.add(acc)
                        fly_ids.add(mid)
                    elif mirbase_id.startswith('hsa'):
                        human_accs.add(acc)
                        human_ids.add(mid)

    # annotate accessions and ids with their database and organism
    for acc in human_accs:
        if mature_re.search(acc):
            graph.add((acc, type_pred, mature_mirna_cl))
        graph.add((acc, organism_pred, homo_iri))
        graph.add((acc, db_pred, mirbase_acc_db))

    for mid in human_ids:
        graph.add((mid, organism_pred, homo_iri))
        graph.add((mid, db_pred, mirbase_id_db))

    for acc in fly_accs:
        if mature_re.search(acc):
            graph.add((acc, type_pred, mature_mirna_cl))
        graph.add((acc, organism_pred, dro_iri))
        graph.add((acc, db_pred, mirbase_acc_db))

    for mid in fly_ids:
        graph.add((mid, organism_pred, dro_iri))
        graph.add((mid, db_pred, mirbase_id_db))

    with open(rdf_path(), 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))


