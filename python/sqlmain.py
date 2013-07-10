
'''
Contains code for creating a sql database and querying it to answer Van Vactor
miRNA collaboration questions.
'''

import argparse
import collections
import contextlib
import csv
import itertools
import logging
import os
import re
import sqlite3
import subprocess
import xml.etree.cElementTree as elementtree

import config
import temps



# Structure of a Dmel Transcript Annotation ID
# http://flybase.org/static_pages/docs/nomenclature/nomenclature3.html#2.4.
# Annotation IDs are represented in a common way: a species-specific 2
# letter prefix followed by a four or five digit integer. For historical
# reasons, there are two 2-letter prefixes for D. melanogaster: CG for
# protein-coding genes and CR for non-protein-coding-genes.
# And from personal communication with Josh Goodman, "The '-RX' and '-PX'
# suffixes (X being one or more letters) are given to transcripts and
# polypeptides respectively."
#
# Notes: CG50-RF is an real annotation id that violates this pattern.
flybase_annotation_id_regex = re.compile(r'^(C(G|R)\d+-R[A-Z]+)$')

# species names used for db tables and microcosm files
mus = 'mus_musculus'
dro = 'drosophila_melanogaster'
cel = 'caenorhabditis_elegans'
homo = 'homo_sapiens'

# ncbi taxon ids for species
mus_taxon = '10090'
dro_taxon = '7227'
cel_taxon = '6239'
homo_taxon = '9606'

# The seven mirs that were screened in muscle and CNS tissue to examine the
# differential gene expression effects.
orginal_screened_mirs = ['dme-miR-34', 'dme-miR-92b', 'dme-miR-137',
                         'dme-miR-190', 'dme-miR-219', 'dme-miR-276a',
                         'dme-miR-277']
# The current mirbase ids for the seven mirs screened by McNeill and Van Vactor.
screened_mirs = ['dme-miR-34-5p', 'dme-miR-92b-3p', 'dme-miR-137-3p',
                 'dme-miR-190-5p', 'dme-miR-219-5p', 'dme-miR-276a-3p',
                 'dme-miR-277-3p']

muscle_tissue = 'muscle'
cns_tissue = 'CNS'
TISSUES = [muscle_tissue, cns_tissue]

# The twenty-six (formerly twenty-seven before the non-existent dme-miR-953 was
# dropped) mirs that were functionally validated by McNeill # and Van Vactor,
# by looking for morphological changes in the muscle phenotype.
original_validated_mirs = [
    'dme-miR-8', 'dme-miR-13a', 'dme-miR-14', 'dme-miR-34', 'dme-miR-92a',
    'dme-miR-92b', 'dme-miR-190', 'dme-miR-137', 'dme-miR-219', 'dme-miR-276a',
    'dme-miR-277', 'dme-miR-279', 'dme-miR-287', 'dme-miR-304', 'dme-miR-308',
    'dme-miR-313', 'dme-miR-314', 'dme-miR-316', 'dme-miR-932', 'dme-miR-969',
    'dme-miR-970', 'dme-miR-978', 'dme-miR-979', 'dme-miR-982', 'dme-miR-999',
    'dme-miR-1014']
# Twenty-six mirs have current miRBase ids.
validated_mirs = [
    u'dme-miR-8-3p', u'dme-miR-13a-3p', u'dme-miR-14-3p', u'dme-miR-34-5p',
    u'dme-miR-92a-3p', u'dme-miR-92b-3p', u'dme-miR-190-5p', u'dme-miR-137-3p',
    u'dme-miR-219-5p', u'dme-miR-276a-3p', u'dme-miR-277-3p',
    u'dme-miR-279-3p', u'dme-miR-287-3p', u'dme-miR-304-5p', u'dme-miR-308-3p',
    u'dme-miR-313-3p', u'dme-miR-314-3p', u'dme-miR-316-5p', u'dme-miR-932-5p',
    u'dme-miR-969-5p', u'dme-miR-970-3p', u'dme-miR-978-3p', u'dme-miR-979-3p',
    u'dme-miR-982-5p', u'dme-miR-999-3p', u'dme-miR-1014-3p']

# The two methods used to assess homology in the 2008 Ibanez-Ventoso paper.
FIVE_PRIME = '5_prime_sequence_homolog'
SEVENTY_PERCENT = '70_percent_full_sequence_homolog'

# taxons and pairs used for roundup orthologs
taxmap = {dro: '7227', homo: '9606', mus: '10090', cel: '6239'}
TAXON_TO_NAME = {cel_taxon: cel, dro_taxon: dro, homo_taxon: homo, mus_taxon: mus}
# Orthologs between human and the other 3 species.
ROUNDUP_PAIRS = [(dro_taxon, homo_taxon), (mus_taxon, homo_taxon), (cel_taxon, homo_taxon)]


###########################
# GENERAL UTILITY FUNCTIONS

def example(message):

    print 'Example CLI function.'
    print 'message:', message


def call(cmd):
    '''
    Run a shell command string using check_call.  Also print the command before
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


def download_and_unzip(url, dest_dir):
    '''
    download a url to a file in dest_dir named with the basename of the url.
    If the file is named *.zip, unzip the file in the directory.  Finally,
    return the filename with any '.zip' suffix returned.
    '''
    dest = os.path.join(dest_dir, os.path.basename(url))
    download(url, dest)
    if dest.endswith('.zip'):
        call('unzip -d {} {}'.format(dest_dir, dest))
    return dest.rstrip('.zip')



#######################################################
# CONSERVED SYNAPSE GENES RELATIONAL DATABASE FUNCTIONS

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


########################
# DATABASE DROP and LOAD

def drop_database():
    path = conserved_synapse_genes_db_path()
    if os.path.exists(path):
        os.remove(path)


def load_all_tables():
    drop_database()
    load_conserved_synapse_genes_tables()
    load_affy_probeset_to_flybase_gene_table()


def load_targetscan_etc():

    # load uniprot to ncbi gene id
    # load fly mir family to gene symbol
    # load fly gene symbol to flybase gene id
    # load fly mir family to mirbase id and mirbase acc

    # load human mir family to ncbi gene id
    # load human mir family to mirbase id and mirbase acc
    return


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


####################
# WRITE OUTPUT FILES

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
        for tissue in TISSUES:
            for mir in screened_mirs:
                fn = os.path.join(config.datadir, '20130523_{}_{}_diff_expr_flybase_genes.txt'.format(mir, tissue))
                probes = [row[0] for row in conn.execute(probe_sql, [tissue, mir])]
                genes = [row[0] for row in conn.execute(gene_sql, [tissue, mir])]
                print mir, tissue, len(probes), 'probes ->', len(genes), 'genes'
                with open(fn, 'w') as fh:
                    for gene in genes:
                        fh.write(gene)
                        fh.write(u'\n')


def write_overlap_between_screened_and_targetscan_predicted_fly_mir_targets():
    raise NotImplementedError()

    predicted_sql = '''select distinct dg2u.gene_id as fly_gene_id
    from drosophila_melanogaster_gene_to_uniprot dg2u
    join drosophila_melanogaster_transcript_to_uniprot dt2u
    join drosophila_melanogaster_annotation_to_transcript da2t
    join drosophila_melanogaster_mir_to_trans dm2t
    where 1
    and dg2u.uniprot_id = dt2u.uniprot_id
    and dt2u.trans_id = da2t.trans_id
    and da2t.annotation_id = dm2t.trans_id
    and dm2t.mir_id = ?
    '''

    screened_sql = '''select distinct ap2g.gene_id as fly_gene_id
    from drosophila_melanogaster_experimental_mir_targets as demt
    join affy_probeset_to_flybase_gene ap2g
    where demt.affy_probeset_id = ap2g.probeset_id
    and demt.mir_id = ?
    and demt.tissue = ?
    '''

    for mir in screened_mirs:
        for tissue in [cns_tissue]:
            fn = os.path.join(config.datadir, '{}_{}_overlap_between_screened_and_targetscan_predicted_targets.csv'.format(tissue, mir))
            with conserved_synapse_genes_db_cm() as conn:
                print mir, tissue
                predicted = set([row[0] for row in conn.execute(predicted_sql, [mir])])
                screened = set([row[0] for row in conn.execute(screened_sql, [mir, tissue])])
                write_set_overlap_file(screened, predicted, 'screened', 'predicted', fn)


def write_overlap_between_screened_and_microcosm_predicted_fly_mir_targets_old():

    predicted_sql = '''select distinct dg2u.gene_id as fly_gene_id
    from drosophila_melanogaster_gene_to_uniprot dg2u 
    join drosophila_melanogaster_transcript_to_uniprot dt2u
    join drosophila_melanogaster_annotation_to_transcript da2t
    join drosophila_melanogaster_mir_to_trans dm2t
    where 1
    and dg2u.uniprot_id = dt2u.uniprot_id
    and dt2u.trans_id = da2t.trans_id
    and da2t.annotation_id = dm2t.trans_id
    and dm2t.mir_id = ?
    '''

    screened_sql = '''select distinct ap2g.gene_id as fly_gene_id
    from drosophila_melanogaster_experimental_mir_targets as demt
    join affy_probeset_to_flybase_gene ap2g
    where demt.affy_probeset_id = ap2g.probeset_id
    and demt.mir_id = ?
    and demt.tissue = ?
    '''

    for mir in screened_mirs:
        for tissue in TISSUES:
            fn = os.path.join(config.datadir, '{}_{}_overlap_between_screened_and_microcosm_predicted_targets.csv'.format(tissue, mir))
            with conserved_synapse_genes_db_cm() as conn:
                print mir, tissue
                predicted = set([row[0] for row in conn.execute(predicted_sql, [mir])])
                screened = set([row[0] for row in conn.execute(screened_sql, [mir, tissue])])
                write_set_overlap_file(screened, predicted, 'screened', 'predicted', fn)


def write_set_overlap_file(set1, set2, name1, name2, filename):
    '''
    Given two sets, write out a csv file containing the set differences and 
    set intersection, placed in 3 columns, in order to be excel friendly.
    The first row is column headers based on name1 and name2.  The columns 
    are <name1>_not_<name2>, <name1>_and_<name2>, <name2>_not_<name1>.

    set1: a set.  Should have csv friendly values, since no quoting of any kind
    is done.
    set2: a set.
    name1: A name for set 1 used to build the column headers.  keep it simple.
    Avoid commas and special characters.
    name2: Like name 1, but for set 2.
    filename: where to write the csv file.
    '''
    one_not_two = sorted(set1 - set2)
    two_not_one = sorted(set2 - set1)
    one_and_two = sorted(set1 & set2)
    makedirs(os.path.dirname(filename))
    with open(filename, 'w') as fh:
        fh.write('{one}_not_{two},{one}_and_{two},{two}_not_{one}\n'.format(
            one=name1, two=name2))
        for row in itertools.izip_longest(one_not_two, one_and_two, two_not_one, fillvalue=''):
            fh.write(','.join([str(i) for i in row]) + '\n')


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


#########
# MIRBASE


def mirbase_aliases_url():
    return 'ftp://mirbase.org/pub/mirbase/19/aliases.txt.gz'


def mirbase_aliases_dest():
    return os.path.join(config.datadir, 'mirbase', '19',
                        os.path.basename(mirbase_aliases_url()))


def mirbase_aliases_filename():
    return mirbase_aliases_dest().rstrip('.gz')


def download_mirbase_aliases():

    url = mirbase_aliases_url()
    dest = mirbase_aliases_dest()
    download(url, dest)
    call('gunzip {}'.format(dest))


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
    Examples:
        P31946  Ensembl ENSG00000166913
        P36872  EnsemblGenome_TRS       FBtr0082186
        P61981  Ensembl_TRS     ENST00000307630
        P20905  FlyBase FBgn0004573
        Q9V785  UCSC    CG7761-RA
        Q2NKQ1  RefSeq_NT       NM_133454.2
    '''
    path = uniprot_idmapping_path(id_type=id_type)
    with open(path) as fh:
        for line in fh:
            yield line.rstrip('\n').split('\t')


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


def gen_microcosm_fly_predicted_targets():
    for target in gen_microcosm_targets(dro):
        mirbase_id = target['mir']
        flybase_transcript_annotation_id = target['transcript_id']
        yield mirbase_id, flybase_transcript_annotation_id


def gen_microcosm_human_predicted_targets():
    for target in gen_microcosm_targets(homo):
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
            if ((species == homo and mir.startswith('hsa-')) 
                or (species == dro and mir.startswith('dme-'))
                or (species != homo and species != dro)):
                yield {'mir': mir, 'transcript_id': transcript_id}


############
# TARGETSCAN


def targetscan_dir(version='6.2'):
    return os.path.join(targetscan_top_dir(), '6.2')


def targetscan_top_dir(version='6.2'):
    return os.path.join(config.datadir, 'targetscan')


def targetscan_fly_symbol_to_gene():
    symbol_to_gene = {}
    for gene_id, symbol in gen_targetscan_fly_gene_to_symbol():
        if symbol in symbol_to_gene:
            raise Exception('Symbol can map to only one gene', symbol, gene_id)
        else:
            symbol_to_gene[symbol] = gene_id

    return symbol_to_gene


def download_targetscan():
    fly_dir = os.path.join(targetscan_dir(), 'fly_12')
    # vert = vertebrates?
    vert_dir = os.path.join(targetscan_dir(), 'vert_61')

    url = 'http://www.targetscan.org/vert_61/vert_61_data_download/miR_Family_Info.txt.zip'
    download_and_unzip(url, vert_dir)
    url = 'http://www.targetscan.org/vert_61/vert_61_data_download/Predicted_Targets_Info.txt.zip'
    download_and_unzip(url, vert_dir)

    url = 'http://www.targetscan.org/fly_12/fly_12_data_download/flybase2symbol_fly12_targetscan.txt'
    download_and_unzip(url, fly_dir)
    url = 'http://www.targetscan.org/fly_12/fly_12_data_download/Conserved_Family_Conserved_Targets_Info.txt.zip'
    download_and_unzip(url, fly_dir)
    url = 'http://www.targetscan.org/fly_12/fly_12_data_download/miR_Family_Info.txt.zip'
    download_and_unzip(url, fly_dir)


def gen_targetscan_txt_fields(path):
    '''
    Yield the fields in each data line in path, skipping the header line.

    The structure of a targetscan txt file is pretty well formed.
    Every file starts with a header line giving the column names.
    Then all the other lines are tab-separated fields.
    '''
    with open(path) as fh:
        for i, line in enumerate(fh):
            # skip headers
            if i == 0:
                continue
            # lines have '\r\n' endings.  Damn you, Windows!
            fields = [f.strip() for f in line.split('\t')]
            yield fields


def gen_targetscan_fly_mir_family_info():
    '''
    There are some family rows with no Mirbase Acc:

    miR-10-3p/1006	AAAUUCG	dme-miR-10-3p	CAAAUUCGGUUCUAGAGAGGUUU	1	
    miR-281-2-5p	AGAGAGC	dme-miR-281-2-5p	AAGAGAGCUAUCCGUCGACAGUC	1	
    miR-2a-1/6/11/13/308	AUCACAG	dme-miR-2a-1	UAUCACAGCCAGCUUUGAUGAGCU	1	
    miR-2a-2/2c	CACAGCC	dme-miR-2a-2	UCACAGCCAGCUUUGAUGAGCUA	1	
    miR-10-5p	CCCUGUA	dme-miR-10-5p	ACCCUGUAGAUCCGAAUUUGUU	1	
    miR-281-1	GUCAUGG	dme-miR-281-1	UGUCAUGGAAUUGCUCUCUUUGU	1	
    miR-281-2-3p	UGUCAUG	dme-miR-281-2-3p	CUGUCAUGGAAUUGCUCUCUUUG	1	
    miR-210.2	UGUGCGU	dme-miR-210.2	UUGUGCGUGUGACAGCGGCUAU	1	
    miR-210.1	UUGUGCG	dme-miR-210.1	CUUGUGCGUGUGACAGCGGCUAU	1	
    miR-iab4as	UACGUAU	dme-miR-iab4as	UUACGUAUACUGAAGGUAUACCG	1

    Some, like 'dme-miR-210.1' are no where to be found in mirbase.  Others, like
    'dme-miR-iab4as' are not in mirbase, but a similar (previous) id, like
    'dme-miR-iab4as-5p', can be found.  Finally, some like 'dme-miR-10-3p' are
    in mirbase and so they should have a mirbase acc.
    '''
    path = os.path.join(targetscan_dir(), 'fly_12', 'miR_Family_Info.txt')
    # Fields: Family members  Seed+m8 MiRBase ID      Mature sequence Family Conservation?    MiRBase Accession
    for fields in gen_targetscan_txt_fields(path):
        family, mirbase_id, mirbase_acc = fields[0], fields[2], fields[5]
        if not mirbase_acc or not mirbase_id:
            print fields
            # Cannot raise exception b/c several rows are missing mirbase_acc.
            # raise Exception(fields)
        else:
            yield family, mirbase_id, mirbase_acc


def gen_targetscan_fly_predicted_targets():
    # Example lines
    # miR Family      Gene Symbol     Gene Tax ID     UTR start       UTR end MSA start       MSA end Seed match
    # miR-927 14-3-3epsilon   7260    91      97      154     161     1A
    # miR-927 14-3-3epsilon   7227    45      51      154     161     1A

    path = os.path.join(targetscan_dir(), 'fly_12',
                        'Conserved_Family_Conserved_Targets_Info.txt')
    for fields in gen_targetscan_txt_fields(path):
        family, gene_symbol, taxon = fields[0], fields[1], fields[2]
        if taxon == '7227':
            if not (family and gene_symbol and taxon):
                raise Exception(fields)
            else:
                yield family, gene_symbol


def gen_targetscan_fly_gene_to_symbol():
    # Example lines
    # gene_id symbol
    # FBgn0000008     a
    # FBgn0000011     ab
    path = os.path.join(targetscan_dir(), 'fly_12',
                        'flybase2symbol_fly12_targetscan.txt')
    for fields in gen_targetscan_txt_fields(path):
        flybase_gene_id, gene_symbol = fields[0], fields[1]
        if not flybase_gene_id or not gene_symbol:
            raise Exception(fields)
        else:
            yield flybase_gene_id, gene_symbol


def gen_targetscan_human_mir_family_info():
    path = os.path.join(targetscan_dir(), 'vert_61', 'miR_Family_Info.txt')
    # Fields: miR family      Seed+m8 Species ID      MiRBase ID      Mature sequence Family Conservation?    MiRBase Accession
    for fields in gen_targetscan_txt_fields(path):
        family, taxon, mirbase_id, mirbase_acc = fields[0], fields[2], fields[3], fields[6]
        if taxon == '9606':
            if not mirbase_acc or not mirbase_id:
                raise Exception(fields)
            else:
                yield family, mirbase_id, mirbase_acc


def gen_targetscan_human_predicted_targets():
    # Example lines
    # miR Family      Gene ID Gene Symbol     Transcript ID   Species ID      UTR start       UTR end MSA start       MSA end Seed match      PCT
    # let-7/98/4458/4500      29974   A1CF    NM_001198819    9913    3248    3255    4796    4809    8mer    0.00
    # let-7/98/4458/4500      29974   A1CF    NM_001198819    9598    3283    3290    4796    4809    8mer    0.00

    path = os.path.join(targetscan_dir(), 'vert_61',
                        'Predicted_Targets_Info.txt')
    for fields in gen_targetscan_txt_fields(path):
        family, refseq_transcript_id, taxon = fields[0], fields[3], fields[4]
        if taxon == '9606':
            if not (family and refseq_transcript_id and taxon):
                raise Exception(fields)
            else:
                yield family, refseq_transcript_id


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


def gen_flybase_transcript_mapping_file():
    '''
    Parse FBtr.xml file and yield FlyBase transcript ids (e.g. FBtr0005674) and
    FlyBase Annotation IDs (e.g.  CG10890-RB) for Drosophila melanogaster.
    Yield tuples like ('7227', 'FBtr0005088', 'CG17291-RA')
    '''
    # The basic structure to parse looks like:
    # <transcript>
        # <id acode="id">FBtr0005674</id>
        # <organism abbreviation="Dmel">...</organism>
        # <dbxref name="FlyBase Annotation IDs" is_current="1" description="" acode="N/A">CG10890-RB</dbxref>
        # ...
    # </transcript>
    version = 'FB2013_02'
    taxon = dro_taxon
    release_dir = flybase_release_dir(version)
    path = os.path.join(release_dir, 'reporting-xml', 'FBtr.xml')
    count = 0
    print path
    # iterate over complete elements
    for event, elem in elementtree.iterparse(path):
        if elem.tag == "transcript":
            # progress meter
            count += 1
            if count % 1000 == 0:
                print 'Processing transcript #', count

            # Organism
            org = elem.find('organism').get('abbreviation')
            # skip non Drosophila melanogaster sequences
            if org != 'Dmel':
                continue

            # Flybase Transcript Id
            tid = elem.find('id').text
            # confirm that id is a flybase transcript id
            if not tid.startswith('FBtr'):
                raise Exception('Expected FlyBase Transcript Id.', tid, org)

            # Flybase Transcript Annotation Ids
            faids = [e.text for e in elem.findall("./dbxref[@name='FlyBase Annotation IDs']")]
            if len(faids) > 1:
                # print 'Multiple FlyBase Annotation IDs.', tid, faids
                pass
            if len(faids) == 0:
                print 'No FlyBase Annotation IDs', tid
            for faid in faids:
                if faid.startswith('FBtr'):
                    # Report and skip instances where a transcript id is in place
                    # of an annotation id
                    msg = 'Flybase Transcript Id found where a Flybase '
                    msg += 'Annotation Id should be. org={}, transcript_id={},'
                    msg += 'annotation_id={}'
                    print msg.format(org, tid, faid)
                elif not flybase_annotation_id_regex.search(faid):
                    # confirm that faid is a Dmel flybase transcript annotation id
                    raise Exception('Expected Dmel Transcript Annotation Id. org={}, transcript_id={}, annotation_id={}'.format(org, tid, faid))

                yield (taxon, tid, faid)
            # try to keep memory from getting too large
            elem.clear()


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

    genes = []
    with open(filename) as fh:
        reader = csv.reader(fh)
        for i, row in enumerate(reader):
            # row 0 is a header row
            if i < 1:
                continue

            gid = row['gid']
            gene_symbol = row['gene_symbol']
            gene_info = row['gene_info']
            ensembl_gene_ids = row['ensembl_gene_id'].split('; ') if row['ensembl_gene_id'] else []
            ensembl_trans_ids = row['ensembl_trans_id'].split('; ') if row['ensembl_trans_id'] else []
            genes.append({'gid': gid, 'gene_symbol': gene_symbol,  'gene_info':
                          gene_info,  'ensembl_gene_ids': ensembl_gene_ids,
                          'ensembl_trans_ids': ensembl_trans_ids})

    # df = pandas.read_csv(filename)
    # for i, row in df.iterrows():
        # gid = row['gid']
        # gene_symbol = row['gene_symbol']
        # gene_info = row['gene_info']
        # ensembl_gene_ids = row['ensembl_gene_id'].split('; ') if row['ensembl_gene_id'] is not np.nan else []
        # ensembl_trans_ids = row['ensembl_trans_id'].split('; ') if row['ensembl_trans_id'] is not np.nan else []
        # genes.append({'gid': gid, 'gene_symbol': gene_symbol,  'gene_info': gene_info,  'ensembl_gene_ids': ensembl_gene_ids, 'ensembl_trans_ids': ensembl_trans_ids})

    return genes


############################
# NMJ_RNAI LOF/GOF GENE LIST

def nmj_rnai_set_path():
    return os.path.join(config.datadir, 'NMJ RNAi Search File.txt')


def get_nmj_rnai_genes():
    '''
    Return a list of flybase gene.

    The hits from three different screens (Aaron D'Antonio, Sanyal and
    Featherstone).  This contains FBgn IDs, which can be converted to gene
    symbols using flybase ID converter 
    '''
    path = nmj_rnai_set_path()
    print path
    with open(nmj_rnai_set_path()) as fh:
        genes = [line.strip() for i, line in enumerate(fh) if i > 0 and line.strip()]
        print genes
    return genes



#####################
# OTHER MAIN FUNCTION

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

    # database dropping, creating, and loading
    subparser = subparsers.add_parser('load_conserved_synapse_genes_tables', help='')
    subparser = subparsers.add_parser('load_affy_probeset_to_flybase_gene_table')
    subparser = subparsers.add_parser('load_experimental_mir_targets_table')

    subparser = subparsers.add_parser('write_targetscan_fly_mirs_targeting_conserved_synapse_genes')
    subparser = subparsers.add_parser('write_targetscan_human_mirs_targeting_conserved_synapse_genes')
    subparser = subparsers.add_parser('write_overlap_between_screened_and_targetscan_predicted_fly_mir_targets')
    subparser = subparsers.add_parser('write_overlap_between_validated_fly_mirs_and_targetscan_predicted_fly_mirs')
    subparser = subparsers.add_parser('write_overlap_between_microcosm_fly_mirs_and_targetscan_fly_mirs')
    subparser = subparsers.add_parser('write_overlap_between_microcosm_human_mirs_and_targetscan_human_mirs')

    subparser = subparsers.add_parser('write_overlap_of_validated_fly_mir_targets_and_nmj_rnai_genes')
    subparser = subparsers.add_parser('write_overlap_of_screened_fly_mir_targets_and_nmj_rnai_genes')
    subparser = subparsers.add_parser('write_overlap_between_screened_and_microcosm_and_targetscan_fly_mir_targets')

    subparser = subparsers.add_parser('write_microcosm_fly_mirs_targeting_conserved_synapse_genes')
    subparser = subparsers.add_parser('write_microcosm_human_mirs_targeting_conserved_synapse_genes')
    subparser = subparsers.add_parser('write_overlap_between_screened_and_microcosm_predicted_fly_mir_targets')
    subparser = subparsers.add_parser('write_overlap_between_validated_fly_mirs_and_microcosm_predicted_fly_mirs')
    subparser = subparsers.add_parser('print_ibanez_fly_human_homologs')
    subparser = subparsers.add_parser('write_experimental_mir_targets_genes_lists')

    subparser = subparsers.add_parser('make_flybase_transcript_map', help='')

    subparser = subparsers.add_parser('map_affymetrix_fly_probe_ids_to_flybase_genes', help='')

    subparser = subparsers.add_parser('write_conserved_gene_and_mir_tables', help='')

    subparser = subparsers.add_parser('download_minotar')
    subparser = subparsers.add_parser('download_targetscan')
    subparser = subparsers.add_parser('download_mirbase_aliases')
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


