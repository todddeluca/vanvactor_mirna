

import os

import rdflib

import config
import temps
from mirna.core import (db_pred, dro_iri, ensembl_transcript_db,
                        ensembl_transcript_iri, flybase_annotation_db,
                        flybase_annotation_iri, homo_iri, mirbase_id_db,
                        mirbase_id_iri, organism_pred, targets_pred)
from mirna.util import call

# Species
# Homo sapiens
HUMAN = '9606'
# Drosophila melanogaster
FLY = '7227'


###########
# MICROCOSM

def microcosm_rdf_path():
    return os.path.join(config.datadir, 'microcosm', 'v5', 'microcosm-v5.nt')


def write_microcosm_rdf():
    print 'write_microcosm_rdf'
    graph = rdflib.Graph()
    fly_mirs = set()
    fly_anno_ids = set()
    human_mirs = set()
    human_ensembl_ids = set()

    # link mir to transcript and collect mirs and transcripts (since each mir
    # or transcript can appear multiple times)
    print 'processing microcosm human mir target predictions'
    for mirbase_id, ensembl_transcript in gen_microcosm_human_predicted_targets():
        mir = mirbase_id_iri(mirbase_id)
        ens = ensembl_transcript_iri(ensembl_transcript)
        graph.add((mir, targets_pred, ens))
        human_mirs.add(mir)
        human_ensembl_ids.add(ens)

    # link each mir to its db and organism
    print 'annotating human mirs'
    for mir in human_mirs:
        graph.add((mir, db_pred, mirbase_id_db))
        graph.add((mir, organism_pred, homo_iri))

    # link each transcript to its db and organism
    print 'annotating human targets'
    for ens in human_ensembl_ids:
        graph.add((ens, db_pred, ensembl_transcript_db))
        graph.add((ens, organism_pred, homo_iri))

    print 'processing microcosm fly mir target predictions'
    for mirbase_id, flybase_annotation_id in gen_microcosm_fly_predicted_targets():
        mir = mirbase_id_iri(mirbase_id)
        anno = flybase_annotation_iri(flybase_annotation_id)
        graph.add((mir, targets_pred, anno))
        fly_mirs.add(mir)
        fly_anno_ids.add(anno)

    print 'annotating fly mirs'
    for mir in fly_mirs:
        graph.add((mir, db_pred, mirbase_id_db))
        graph.add((mir, organism_pred, dro_iri))

    print 'annotating fly targets'
    for anno in fly_anno_ids:
        graph.add((anno, db_pred, flybase_annotation_db))
        graph.add((anno, organism_pred, dro_iri))

    with open(microcosm_rdf_path(), 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))


def download_microcosm_targets(species, version='v5'):
    '''
    species: e.g. mus_musculus
    '''
    species = species.lower() # e.g. 'mus_musculus'
    urlbase = 'ftp://ftp.ebi.ac.uk/pub/databases/microcosm'
    url = urlbase + '/{version}/arch.{version}.txt.{species}.zip'.format(
        species=species, version=version)
    dest = microcosm_species_targets_file(species=species, version=version)

    if os.path.exists(dest):
        raise Exception('Targets already downloaded.', species, dest)

    with temps.tmpdir() as td:
        # download path
        dpath = os.path.join(td, os.path.basename(url))
        # download
        call('curl -o {} {}'.format(dpath, url))
        # unzip
        call('unzip -d {} {}'.format(td, dpath))
        unzip_path = os.path.join(td, '{}.txt.{}'.format(version, species))
        # move to download dir
        os.rename(unzip_path, dest)


def microcosm_v5_mouse_targets_file():
    return microcosm_species_targets_file('mus_musculus')


def microcosm_v5_human_targets_file():
    return microcosm_species_targets_file('homo_sapiens')


def microcosm_species_targets_file(species, version='v5'):
    return os.path.join(config.microcosm_dir, version,
                        'microcosm-{version}-{species}_targets.txt'.format(
                            version=version, species=species))


def gen_microcosm_fly_predicted_targets():
    for target in gen_microcosm_targets(FLY):
        mirbase_id = target['mir']
        flybase_transcript_annotation_id = target['transcript_id']
        yield mirbase_id, flybase_transcript_annotation_id


def gen_microcosm_human_predicted_targets():
    for target in gen_microcosm_targets(HUMAN):
        mirbase_id = target['mir']
        ensembl_transcript_id = target['transcript_id']
        yield mirbase_id, ensembl_transcript_id


def gen_microcosm_targets(species, version='v5'):
    '''
    Yield every row of the microcosm targets table as a dict containing 'mir'
    and 'transcript_id'.

    species: a microcosm species name, e.g. 'homo_sapiens' or 'drosophila_melanogaster'
    '''
    path = microcosm_species_targets_file(species, version)
    with open(path) as fh:
        for line in fh:
            # skip comments and blank lines
            if not line.strip() or line.strip().startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            assert len(fields) == 13
            mir = fields[1]
            transcript_id = fields[11]
            # For some reason, the human microcosm data has non-human mirs in
            # it.
            # Do not yield non-human mirs for Homo sapiens (or non dme mirs for
            # Drosophila melanogaster).
            if ((species == HUMAN and mir.startswith('hsa-')) 
                or (species == FLY and mir.startswith('dme-'))
                or (species != HUMAN and species != FLY)):
                yield {'mir': mir, 'transcript_id': transcript_id}


