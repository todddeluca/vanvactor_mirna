
import argparse
import collections
import contextlib
import logging
import os
import re
import sqlite3
import subprocess
import sys
import xml.etree.cElementTree as elementtree
import csv

import pandas
import numpy as np

import config
import temps


# species names used for db tables and microcosm files
mus = 'mus_musculus'
dro = 'drosophila_melanogaster'
cel = 'caenorhabditis_elegans'
homo = 'homo_sapiens'

# The seven mirs that were screened in muscle and CNS tissue to examine the
# differential gene expression effects.
screened_mirs = ['dme-miR-34', 'dme-miR-92b', 'dme-miR-137', 'dme-miR-190',
                 'dme-miR-219', 'dme-miR-276a', 'dme-miR-277']

# The twenty-seven mirs that were somehow "functionally validated" by McNeill
# and Van Vactor.
validated_mirs = ['dme-miR-8', 'dme-miR-13a', 'dme-miR-14', 'dme-miR-34',
                  'dme-miR-92a', 'dme-miR-92b', 'dme-miR-190', 'dme-miR-137',
                  'dme-miR-219', 'dme-miR-276a', 'dme-miR-277', 'dme-miR-279',
                  'dme-miR-287', 'dme-miR-304', 'dme-miR-308', 'dme-miR-313',
                  'dme-miR-314', 'dme-miR-316', 'dme-miR-932', 'dme-miR-953',
                  'dme-miR-969', 'dme-miR-970', 'dme-miR-978', 'dme-miR-979',
                  'dme-miR-982', 'dme-miR-999', 'dme-miR-1014']

FIVE_PRIME = '5_prime_sequence_homolog'
SEVENTY_PERCENT = '70_percent_full_sequence_homolog'


def example(message):

    print 'Example CLI function.'
    print 'message:', message


def call(cmd):
    '''
    Run a shell command using check_call.  Also print the command before
    running it.
    '''
    print cmd
    return subprocess.check_call(cmd, shell=True)


def makedirs(path, mode=0775):
    '''
    Idempotent function for making a directory.  If path is not a directory,
    make path and all of its non-existent parent directories.  Either way,
    return path.
    '''
    if not os.path.exists(path):
        os.makedirs(path, mode)
    return path


def download(url, dest, mode=0775):
    '''
    dest: the filename to download the url to.  Will create any parent
    directories of dest that do not already exist.
    '''
    makedirs(os.path.dirname(dest), mode)
    call('curl -o {} {}'.format(dest, url))
    return dest


def human_synaptic_mir_targets():
    pass


def fly_orthologs_of_human_synaptic_mir_targets():
    pass


def fly_human_mir_targets_and_synaptic_genes():
    pass


def conserved_synapse_genes_db_path():
    # dbpath = ':memory:' # in memory database for testing.
    dbpath = os.path.join(config.datadir, '20130330_dro_mus_homo_cel_mapping.db')
    return dbpath


@contextlib.contextmanager
def conserved_synapse_genes_db_cm():
    '''
    Yield an open connnection for use in a with statement.  Commit and close
    the connection when exiting the with statement.
    '''
    path = conserved_synapse_genes_db_path()
    conn = sqlite3.connect(path)
    try:
        yield conn
    except:
        conn.rollback()
        raise
    else:
        conn.commit()
    finally:
        conn.close()


def table_head(table, conn):
    '''
    Print the first few rows of table.
    '''
    sql = 'select * from {} limit 10'.format(table)
    print sql
    for row in conn.execute(sql):
        print row


def write_conserved_gene_and_mir_tables():

    path = conserved_synapse_genes_db_path()
    conn = sqlite3.connect(path)

    def write_table(sql, path, fields):
        print 'write', sql, 'to', path, 'with fields', fields
        with open(path, 'w') as fh:
            fh.write('# The first row is a header row.  All columns are tab-separated.\n')
            fh.write('\t'.join(fields) + '\n')
            for row in conn.execute(sql):
                fh.write('\t'.join(row) + '\n')

    # mir to gene for homo, mus, cel
    for s, gene_id_type in [
        (homo, 'ensembl_gene_id'), 
        (mus, 'ensembl_gene_id'), 
        (cel, 'wormbase_id'),
    ]:
        sql = ' '.join([
            'select m2t.mir_id, g2u.gene_id',
            'from {org}_gene_to_uniprot g2u',
            'join {org}_transcript_to_uniprot t2u',
            'join {org}_mir_to_trans m2t',
            'where t2u.uniprot_id = g2u.uniprot_id',
            'and m2t.trans_id = t2u.trans_id',
        ]).format(org=s)
        path = os.path.join(config.datadir, 
                            '20130330_{}_mir_to_gene.tsv'.format(s))
        write_table(sql, path, ['microcosm_mir_id', gene_id_type])

    # mir to gene for fly
    sql = ' '.join([
        'select m2t.mir_id, g2u.gene_id',
        'from {org}_gene_to_uniprot g2u',
        'join {org}_transcript_to_uniprot t2u',
        'join {org}_annotation_to_transcript a2t',
        'join {org}_mir_to_trans m2t',
        'where t2u.uniprot_id = g2u.uniprot_id',
        'and a2t.trans_id = t2u.trans_id',
        'and a2t.annotation_id = m2t.trans_id',
    ]).format(org=dro)
    path = os.path.join(config.datadir, 
                        '20130330_{}_mir_to_gene.tsv'.format(dro))
    write_table(sql, path, ['microcosm_mir_id', 'flybase_gene_id'])

    # gene to gene for human vs (fly or mouse or worm)
    for s, query_id_type in [
        (mus, 'ensembl_gene_id'),
        (dro, 'flybase_gene_id'),
        (cel, 'wormbase_gene_id'),
    ]:
        sql = ' '.join([
            'select query_g2u.gene_id, subject_g2u.gene_id',
            'from {org}_to_{homo}_orthologs ologs',
            'join {org}_gene_to_uniprot query_g2u',
            'join {homo}_gene_to_uniprot subject_g2u',
            'where query_g2u.uniprot_id = ologs.query_id',
            'and subject_g2u.uniprot_id = ologs.subject_id',
        ]).format(org=s, homo=homo)
        path = os.path.join(config.datadir, 
                            '20130330_{}_to_{}_gene_orthologs.tsv'.format(s, homo))
        write_table(sql, path, [query_id_type, 'ensembl_gene_id'])

    # human synapse genes
    sql = 'select id, symbol, desc from {homo}_synapse_genes'.format(homo=homo)
    path = os.path.join(config.datadir, 
                        '20130330_{}_synapsedb_genes.tsv'.format(homo))
    write_table(sql, path, ['ensembl_gene_id', 'gene_symbol', 'gene_description'])


def write_experimental_mir_targets_genes_lists():
    # select the genes that are differentially expressed in a tissue when a mir is perturbed.
    gene_sql = '''select distinct gene_id 
             from drosophila_melanogaster_experimental_mir_targets 
             join affy_probeset_to_flybase_gene 
             on affy_probeset_id = probeset_id 
             where tissue = ? and mir_id = ? '''
    # select the probeset ids that are differentially expressed in a tissue when a mir is perturbed.
    probe_sql = '''select distinct affy_probeset_id
             from drosophila_melanogaster_experimental_mir_targets 
             where tissue = ? and mir_id = ? '''
    
    with conserved_synapse_genes_db_cm() as conn:
        for tissue in ['muscle', 'CNS']:
            for mir in screened_mirs:
                fn = os.path.join(config.datadir, '20130523_{}_{}_diff_expr_flybase_genes.txt'.format(mir, tissue))
                probes = [row[0] for row in conn.execute(probe_sql, [tissue, mir])]
                genes = [row[0] for row in conn.execute(gene_sql, [tissue, mir])]
                print mir, tissue, len(probes), 'probes ->', len(genes), 'genes'
                with open(fn, 'w') as fh:
                    for gene in genes:
                        fh.write(gene)
                        fh.write(u'\n')


def load_2008_fly_human_mir_homologs_table():
    '''
    Create a db table containing the fly mirs homologous to human mirs, based
    on the paper "Sequence Relationships among C. elegans, D. melanogaster and
    Human microRNAs Highlight the Extensive Conservation of microRNAs in
    Biology" by Carolina Ibanez-Ventoso, Mehul Vora, and Monica Driscoll.

    The table contains a column for the drosophila mir id, the homo mir id, and
    the method used to assert the homologous relationship.  The two methods
    used in the paper were 5-prime sequence homology and greater than 70% full
    sequence homology.
    '''
    with conserved_synapse_genes_db_cm() as conn:
        table = '{}_{}_mir_homologs'.format(dro, homo)
        conn.execute('drop table if exists {}'.format(table))
        conn.execute('create table {} (dro_mir_id, homo_mir_id, method)'.format(table))
        params = [(dro_mir, homo_mir, method) for dro_mir, homo_mir, method in
                  gen_2008_fly_human_homologs()]
        conn.executemany('insert into {} values(?, ?, ?)'.format(table), params)
        table_head(table, conn)


def load_experimental_mir_targets_table():
    '''
    Create a db table containing the mirs, differentially expressed affy
    probeset ids, and tissue condition of the experiments done by Elizabeth
    McNeill and Davie Van Vactor.
    '''
    with conserved_synapse_genes_db_cm() as conn:
        table = '{}_experimental_mir_targets'.format(dro)
        conn.execute('drop table if exists {}'.format(table))
        conn.execute('create table {} (mir_id, tissue, affy_probeset_id)'.format(table))
        params = [(mir, tissue, probe) for mir, tissue, probe in 
                gen_drosophila_experimental_mir_targets()]
        conn.executemany('insert into {} values(?, ?, ?)'.format(table), params)
        table_head(table, conn)


def load_affy_probeset_to_flybase_gene_table():
    with conserved_synapse_genes_db_cm() as conn:
        table = 'affy_probeset_to_flybase_gene'
        conn.execute('drop table if exists {}'.format(table))
        conn.execute('create table {} (probeset_id, gene_id)'.format(table))
        params = [(probe, gene) for probe, gene in
                  gen_affymetrix_fly_probe_to_flybase_mapping()]
        conn.executemany('insert into {} values(?, ?)'.format(table), params)
        table_head(table, conn)


def load_conserved_synapse_genes_tables():
    '''
    For fly, mouse, human, and worm, create tables and load
    data for:

    - Microcosm mir to transcript targets.
    - Uniprot mapping of transcript ids to uniprot accessions.
    - Uniprot mapping of gene ids to uniprot accessions.
    - Roundup mappings of human uniprot ids to each of the other organisms.

    Also, since the Flybase Annotation IDs used as transcript ids by microcosm
    can not be translated by UniProt, create and load a table to map annotation
    ids to FlyBase Transcript Ids.

    Finally, create and load a table containing human synapse genes.

    For these tables, Ensembl ids are used for human and mouse genes and
    transcripts.  Flybase ids are used for fly genes and transcripts, except
    for the microcosm data.  And WormBase IDs are used for worm genes and
    transcripts.
    '''
    path = conserved_synapse_genes_db_path()
    if os.path.exists(path):
        os.remove(path)
    conn = sqlite3.connect(path)
    c = conn.cursor()

    specs = [dro, homo, mus, cel]
    # taxons and pairs used for roundup orthologs
    taxmap = {dro: '7227', homo: '9606', mus: '10090', cel: '6239'}
    ortholog_pairs = [(dro, homo), (mus, homo), (cel, homo)]

    # uniprot id mapping types for mapping from transcript to uniprotkb-ac (what I call uniprot ids)
    uniprot_mus_trans_type = 'Ensembl_TRS' # e.g. ENSMUST00000034287
    uniprot_homo_trans_type = 'Ensembl_TRS' # e.g. ENST00000295228
    uniprot_dro_trans_type = 'EnsemblGenome_TRS' # e.g. 'FBtr0082186', not 'CG11023-RA'
    uniprot_cel_trans_type = 'WormBase_TRS' # e.g. Y46E12A.4

    # mapping id types for mapping from gene to uniprotkb-ac (what I call uniprot ids)
    uniprot_mus_gene_type = 'Ensembl' # e.g. ENSMUSG00000047281
    uniprot_homo_gene_type = 'Ensembl' # e.g. ENSG00000166913
    uniprot_cel_gene_type = 'WormBase' # e.g. WBGene00017178
    uniprot_dro_gene_type = 'FlyBase' # e.g. FBgn0000008

    # table mapping flybase annotation ids to flybase transcript ids
    dro_annot_to_tran_table = '{}_annotation_to_transcript'.format(dro)
    # table mapping mir to transcript id
    def microcosm_table_name(species):
        return '{}_mir_to_trans'.format(species)

    def head(table):
        sql = 'select * from {} limit 10'.format(table)
        print sql
        for row in c.execute(sql):
            print row

    # create and load gene to uniprot tables
    for s, id_type, regex in [
        (mus, uniprot_mus_gene_type, re.compile(r'^ENSMUSG\d')),
        (homo, uniprot_homo_gene_type, re.compile(r'^ENSG\d')),
        (cel, uniprot_cel_gene_type, re.compile(r'.*')),
        (dro, uniprot_dro_gene_type, re.compile(r'.*')),
    ]:
        print 'loading', s, 'gene-to-uniprot table'
        table = '{}_gene_to_uniprot'.format(s)
        c.execute('create table {} (gene_id, uniprot_id)'.format(table))
        params = [(mapped, uniprot) for uniprot, id_type, mapped in
                  gen_uniprot_id_mappings(id_type=id_type) if
                  regex.search(mapped)]
        c.executemany('insert into {} values(?, ?)'.format(table), params)
        head(table)

    # create and load mir-to-transcript tables
    for s in specs:
        print 'loading', s, 'mir targets'
        table = microcosm_table_name(s)
        c.execute('create table {} (mir_id, trans_id)'.format(table))
        params = [(t['mir'], t['transcript_id']) for t in
                  gen_microcosm_targets(species=s)]
        c.executemany('insert into {} values (?, ?)'.format(table), params)
        head(table)

    # create and load fly annotation id transcript flybase transcript table
    print 'loading', dro, 'annotation id to transcript id table'
    table = dro_annot_to_tran_table
    c.execute('create table {} (annotation_id, trans_id)'.format(table))
    params = [(aid, tid) for org, tid, aid in
              gen_flybase_transcript_to_annotation_mapping() if org == 'Dmel']
    c.executemany('insert into {} values (?, ?)'.format(table), params)
    head(table)

    # create and load transcript-to-uniprot tables
    # human, mouse, and worm microcosm transcripts can all be translated
    # directly to uniprot ids from the transcript ids used by microcosm.
    # Microcosm uses flybase annotation ids, which need to be converted into
    # flybase transcript ids first.
    for s, id_type, trans_table, regex in [
        (mus, uniprot_mus_trans_type, microcosm_table_name(mus), re.compile(r'^ENSMUST\d')),
        (homo, uniprot_homo_trans_type, microcosm_table_name(homo), re.compile(r'^ENST\d')),
        (cel, uniprot_cel_trans_type, microcosm_table_name(cel), re.compile(r'.*')),
        (dro, uniprot_dro_trans_type, dro_annot_to_tran_table, re.compile(r'^FBtr\d')),
    ]:
        print 'loading', s, 'transcript-to-uniprot table'
        table = '{}_transcript_to_uniprot'.format(s)
        c.execute('create table {} (trans_id, uniprot_id)'.format(table))
        params = [(mapped, uniprot) for uniprot, id_type, mapped in
                  gen_uniprot_id_mappings(id_type=id_type) if
                  regex.search(mapped)]
        c.executemany('insert into {} values(?, ?)'.format(table), params)
        head(table)


    # create and load ortholog tables
    for pair in ortholog_pairs:
        print 'loading', pair, 'orthologs'
        table = '{}_to_{}_orthologs'.format(*pair)
        c.execute('create table {} (query_id, subject_id, distance)'.format(table))
        params = [(qid, sid, dist) for qid, sid, dist in
                  gen_orthologs(taxmap[pair[0]], taxmap[pair[1]])]
        c.executemany('insert into {} values (?, ?, ?)'.format(table), params)
        head(table)

    # create and load human synapse table
    print 'loading human synapse genes'
    table = '{}_synapse_genes'.format(homo)
    c.execute('create table {} (id primary key, symbol, desc)'.format(table))
    genes = parse_synaptomedb_all_genes()
    rows = {}
    for g in genes:
        for ensembl_g in g['ensembl_gene_ids']:
            row = (ensembl_g, g['gene_symbol'], g['gene_info'])
            if ensembl_g in rows:
                print 'Existing ensembl gene and duplicate:'
                print rows[ensembl_g]
                print row
            else:
                rows[ensembl_g] = row

    params = rows.values()
    c.executemany('insert into {} values (?, ?, ?)'.format(table), params)
    head(table)

    conn.commit()
    conn.close()


########################
# FLY-HUMAN MIR HOMOLOGS


def gen_2008_fly_human_homologs():
    '''
    From PLoS ONE. 2008; 3(7): e2818.
    Published online 2008 July 30. doi:  10.1371/journal.pone.0002818
    PMCID: PMC2486268
    Sequence Relationships among C. elegans, D. melanogaster and Human microRNAs Highlight the Extensive Conservation of microRNAs in Biology
    Carolina Ibanez-Ventoso, Mehul Vora, and Monica Driscoll

    Generate tuples of drosophila mir id, homo mir id, method) for the ">= 70%
    full sequence homology" and "5-prime sequence homology" methods for
    asserting homology between fly and human mirs from the paper.
    '''
    fn = os.path.join(config.datadir, '2008_Ibanez-Ventoso_drosophila_melanogaster_homo_sapiens_mir_orthologs.csv')
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
                yield current_dro_mir, full_homo_mir, SEVENTY_PERCENT


def print_2008_fly_human_homologs():
    '''
    See what the generated tuples of human-fly mir homologs look like to make
    sure they do not look wrong.
    '''
    for dro_mir, homo_mir, method in gen_2008_fly_human_homologs():
        print method, dro_mir, homo_mir



###################
# ROUNDUP FUNCTIONS


def gen_orthologs(query_taxon, subject_taxon, divergence='0.8', evalue='1e-5',
                  version='4'):
    d = os.path.join(config.datadir, 'roundup', 'v' + version)
    path = os.path.join(d, 'roundup-{}-orthologs_for_{}_{}_{}_{}.txt'.format(
        version, query_taxon, subject_taxon, divergence, evalue))
    with open(path) as fh:
        for line in fh:
            if line.startswith('OR'):
                kind, qid, sid, distance = line.rstrip('\n').split('\t')
                distance = float(distance)
                yield qid, sid, distance


###################
# UNIPROT FUNCTIONS


def uniprot_idmapping_path(version='2013_04', id_type=None):
    d = os.path.join(config.datadir, 'uniprot', version)
    if id_type is None:
        return os.path.join(d, 'idmapping.dat')
    else:
        return os.path.join(d, 'idmapping.{}.dat'.format(id_type))


def split_uniprot_idmapping_file():
    '''
    Split the idmapping.dat file into one file for each id_type.  This 
    allows more efficient searching for ids that having everything in one
    big file.  I really should load it into a database like sqlite or  berkeley
    db.
    '''
    fhs = {} # id type filehandles
    with open(uniprot_idmapping_path()) as fh:
        for i, line in enumerate(fh):
            if i % 10000000 == 0:
                print 'line', i
            id_type = line.split('\t')[1]
            if id_type not in fhs:
                fhs[id_type] = open(uniprot_idmapping_path(id_type=id_type), 'w')
            fhs[id_type].write(line)


def gen_uniprot_id_mappings(id_type=None):
    '''
    Return tuples of "UniProtKB-AC", "ID type", and "ID"
    '''
    path = uniprot_idmapping_path(id_type=id_type)
    with open(path) as fh:
        for line in fh:
            yield line.rstrip('\n').split('\t')


def guess_uniprot_mapping_idtype(identifer):
    '''
    Given an identifier of an unknown type, find the first (if any) id that
    matches the identifier in the uniprot idmapping file.  Return the type
    of that id.

    This can be useful for figuring out the exact name of the id type
    of your ids.
    '''
    for uniprot, id_type, mapped in gen_uniprot_id_mappings():
        if identifier == mapped:
            print id_type
            return id_type


def map_ids_to_uniprot(ids, id_type):
    '''
    Map each id in ids to a (possibly empty) list of UniProtKB-AC (a.k.a.
    uniprot accessions or uniprot ids) that are mapped to the id.

    ids: a list (or set) of ids.
    id_type: the type of the id in the uniprot id mapping file.
    '''
    idmap = {i: [] for i in ids}
    for uniprot, id_type, mapped in gen_uniprot_id_mappings(id_type=id_type):
        if mapped in idmap:
            idmap[mapped].append(uniprot)

    return idmap


#################################
# VAN VACTOR / MCNEILL EXPERIMENT


def gen_drosophila_experimental_mir_targets():
    '''
    Generate a tuple for every row in the CSV file containing the experimental
    results for affymetrix probesets that were differentially expressed in
    drosophila 'muscle' and 'CNS' tissues when 7 miRs were (individually)
    perturbed.
    '''
    fn = os.path.join(config.datadir, '20130523_vanvactor_fly_tissue_7_mir_targets.csv')
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



############
# AFFYMETRIX


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


###########
# MICROCOSM

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


def gen_microcosm_targets(species, version='v5'):
    '''
    Yield every row of the microcosm targets table as a dict containing 'mir'
    and 'transcript_id'.
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
            yield {'mir': mir, 'transcript_id': transcript_id}



###################
# FLYBASE FUNCTIONS


def flybase_transcript_to_annotation_id_table_path(version='FB2013_02'):
    d = flybase_release_dir(version)
    f = 'flybase_{}_transcript_to_annotation_id_table.txt'.format(version)
    return os.path.join(d, f)


def flybase_release_dir(version='FB2013_02'):
    return os.path.join(config.datadir, 'flybase', 'releases', version)


def gen_flybase_transcript_to_annotation_mapping(version='FB2013_02'):
    '''
    Yield tuples of (organism_abbreviation, transcript_id, annotation_id).
    For example:

        ('Dmel', 'FBtr0340478', 'CG3315-RB')
    '''
    with open(flybase_transcript_to_annotation_id_table_path(version)) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            organism, trans_id, annot_id = stripped.split('\t')
            yield organism, trans_id, annot_id


def make_flybase_transcript_mapping_file():
    '''
    Parse FBtr.xml file and write out a tab-separated map from FlyBase
    transcript ids (e.g. FBtr0005674) and FlyBase Annotation IDs (e.g.
    CG10890-RB).  Output lines like:

        Dmel    FBtr0005088     CG17291-RA
        Dmel    FBtr0005674     CG10890-RB
    '''
    # The basic structure to parse looks like:
    # <transcript>
        # <id acode="id">FBtr0005674</id>
        # <organism abbreviation="Dmel">...</organism>
        # <dbxref name="FlyBase Annotation IDs" is_current="1" description="" acode="N/A">CG10890-RB</dbxref>
        # ...
    # </transcript>
    version = 'FB2013_02'
    release_dir = flybase_release_dir(version)
    path = os.path.join(release_dir, 'reporting-xml', 'FBtr.xml')
    outfile = flybase_transcript_to_annotation_id_table_path(version)
    count = 0
    with open(outfile, 'w') as fh:
        fh.write('# Generated from FlyBase release {}.\n'.format(version))
        fh.write('# Derived from the reporting-xml FBtr.xml file.\n')
        fh.write('# Tab-separated fields: organism_abbreviation, flybase_transcript_id, flybase_annotation_id\n')

        for event, elem in elementtree.iterparse(path):
            if event == 'end' and elem.tag == "transcript":
                count += 1
                if count % 1000 == 0:
                    print 'Processing transcript #', count
                org = elem.find('organism')
                tid = elem.find('id')
                faids = elem.findall("./dbxref[@name='FlyBase Annotation IDs']")
                if len(faids) > 1:
                    print 'Multiple FlyBase Annotation IDs'
                if len(faids) == 0:
                    print 'No FlyBase Annotation IDs'
                for faid in faids:
                    fh.write('{}\t{}\t{}\n'.format(org.get('abbreviation'), tid.text, faid.text))
                # keep memory from getting too large
                elem.clear()



def flybase_genes(gene_conversion_table, outfile):
    '''
    Return a list of unique flybase gene ids from the gene_conversion table,
    skipping ids that did not convert.

    outfile: a file path where the genes will be written.  If None, genes will
    be written to standard output.


    '''
    genes = select_flybase_gene_ids(gene_conversion_table)
    with open(outfile, 'w') as fh:
        for gene in genes:
            fh.write('{}\n'.format(gene))


def gen_flybase_gene_conversion(gene_conversion_table):
    '''
    Given a conversion table that is converting "weird" transcript ids, e.g. CG10005-RA,
    to flybase transcript ids, e.g. FBtr0082507, yield tuples mapping the
    original ids to the transcript ids, e.g. ('CG10005-RA', 'FBtr0082507').
    Do not yield ids that fail to convert.
    '''
    for i, line in enumerate(open(gene_conversion_table)):
        # skip comments and blank lines
        if not line.strip() or line.strip().startswith('#'):
            continue

        # Example line with 5 columns with a flybase gene id
        # CG10005-RA      FBtr0082507             FBgn0037972     CG10005
        # Example line with 3 columns and no flybase gene id
        # CG6149-RA       unknown ID      -
        # Example line with 4 columns and no flybase gene id
        # CG6151-RC       FBtp0052133     -       -
        # Fields: Submitted ID, Current ID, mystery field, Converted ID, Related record
        splits = line.strip().split("\t")

        # skip values that flybase failed to convert to a gene id
        if splits[1] == 'unknown ID' or splits[3] == '-':
            continue

        # assuming this is a gene conversion table, then flybase converted the
        # submitted id into a flybase gene id.
        orig, gene = splits[0], splits[3]
        assert gene.startswith("FBgn")
        yield orig, gene


def select_flybase_gene_ids(gene_conversion_table):
    '''
    Return a list of unique flybase gene ids from the gene conversion table
    downloaded from flybase, skipping ids that did not convert.
    '''
    uniques = set()
    for orig, gene in gen_flybase_gene_conversion(gene_conversion_table):
        uniques.add(gene)

    return sorted(uniques)


#######################
# SYNAPTOMEDB FUNCTIONS


def synaptomedb_v1_all_genes_file():
    return os.path.join(config.datadir, 'synaptomedb', 'v1.06',
                        'synaptomedb-1.06-all_genes.csv')


def parse_synaptomedb_all_genes(filename=None):
    '''
    Return a list of genes.  Each gene is a dict containing 'gid' (a synaptome
    gene id?), 'gene_symbol', 'gene_info' (a description of the gene),
    'ensembl_gene_ids' (a list of ensembl gene ids), and 'ensembl_trans_ids' (a
    list of ensembl transcript ids).
    '''
    if filename is None:
        filename = synaptomedb_v1_all_genes_file()

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
    targets = list(gen_microcosm_targets('homo_sapiens'))
    # print 'targets', targets

    # quickly look up the mirs that target a transcript
    trans2mirs = collections.defaultdict(set)
    for target in targets:
        # transcripts are ensembl transcript ids
        trans2mirs[target['transcript_id']].add(target['mir'])

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


def get_uniprot_accs_from_idmapping_table(mapping_table):
    '''
    mapping table: a file like:

        # tab-separated fields: flybase_gene_id, uniprotkb_acc
        FBgn0000055     P00334
        FBgn0000083     A4V371
        FBgn0000083     P22464
        ...
        FBgn0263278     
    '''
    with open(mapping_table) as fh:
        uniprots = set()
        for line in fh:
            # skip blank and commment lines
            if not line.strip() or line.strip().startswith('#'):
                continue
            from_id, to_id = line.rstrip('\n').split('\t')
            if to_id:
                uniprots.add(to_id)

    # output the unique uniprot ids
    for uid in sorted(uniprots):
        print uid



def merge_uniprot_idmapping(from_id_type, to_id_type, mapping_table, not_mapped_list):
    '''
    mapping_table:  a file like:

        From    To
        FBgn0000083     A4V371
        FBgn0000083     P22464

    not_mapped_list: a file like:

        not mapped
        FBgn0058441
        FBgn0263278
    '''


    with open(mapping_table) as mapped_fh, open(not_mapped_list) as not_fh:
        print '# tab-separated fields: {}, {}'.format(from_id_type, to_id_type)
        for i, line in enumerate(mapped_fh):
            if i == 0:
                assert line == 'From\tTo\n', 'Unrecognized header line in uniprot id mapping table {}'.format(mapping_table)
            else:
                # make sure line parses
                fid, tid = line.strip().split('\t')
                print line.rstrip('\n')
        for i, line in enumerate(not_fh):
            if i == 0:
                assert line == 'not mapped\n', 'Unrecognized header line in uniprot id mapping not mapped list {}'.format(not_mapped_list)
            else:
                fid = line.strip()
                print fid + '\t'




def main():

    parser = argparse.ArgumentParser(description='')
    subparsers = parser.add_subparsers(dest='action', help='')

    subparser = subparsers.add_parser('example', help='Example help text')
    subparser.add_argument('message', help='print this')

    subparser = subparsers.add_parser('get_uniprot_accs_from_idmapping_table', help='')
    subparser.add_argument('mapping_table')

    subparser = subparsers.add_parser('merge_uniprot_idmapping', help='')
    subparser.add_argument('from_id_type')
    subparser.add_argument('to_id_type')
    subparser.add_argument('mapping_table')
    subparser.add_argument('not_mapped_list')

    subparser = subparsers.add_parser('load_2008_fly_human_mir_homologs_table')
    subparser = subparsers.add_parser('print_2008_fly_human_homologs')
    subparser = subparsers.add_parser('write_experimental_mir_targets_genes_lists')
    subparser = subparsers.add_parser('load_affy_probeset_to_flybase_gene_table')
    subparser = subparsers.add_parser('load_experimental_mir_targets_table')
    subparser = subparsers.add_parser('load_conserved_synapse_genes_tables', help='')

    subparser = subparsers.add_parser('make_flybase_transcript_map', help='')

    subparser = subparsers.add_parser('map_affymetrix_fly_probe_ids_to_flybase_genes', help='')

    subparser = subparsers.add_parser('write_conserved_gene_and_mir_tables', help='')

    subparser = subparsers.add_parser('download_microcosm_targets', help='')
    subparser.add_argument('species')

    subparser = subparsers.add_parser('flybase_genes', help='Extract flybase')
    subparser.add_argument('gene_conversion_table')
    subparser.add_argument('outfile')

    subparser = subparsers.add_parser('map_human_synaptic_genes_to_human_mirs', help='')

    # invoke a function named by action whose keyword parameters correspond to
    # cli arguments and whose parser name corresponds to function name
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


