
'''
Contains code for creating a semantic database and querying it to answer Van Vactor
miRNA collaboration questions.
'''

import argparse
import collections
import csv
import itertools
import json
import logging
import os
import re
import subprocess
import urllib
import xml.etree.cElementTree as elementtree
import StringIO

import rdflib
import requests

import config
import secrets
import temps
import virtuoso
# import stardog
import sparql


virt7 = virtuoso.Virtuoso7(secrets.virtuoso_dba_password,
                           config.virtuoso_load_dir)
sparq = sparql.Sparql(virt7.sparql_endpoint()) # used for querying


# RDF URIs and URI generators
# Real Entities
dro_iri = rdflib.URIRef('http://purl.uniprot.org/taxonomy/7227')
homo_iri = rdflib.URIRef('http://purl.uniprot.org/taxonomy/9606')
muscle_tissue_iri = rdflib.URIRef('http://purl.uniprot.org/tissues/642')
cns_tissue_iri = rdflib.URIRef('http://purl.uniprot.org/tissues/150')
# Fake Entities
synaptic_iri = rdflib.URIRef('http://purl.synaptomedb.org/owl/synaptic')
five_prime_method_iri = rdflib.URIRef('http://purl.example.com/ibanez/five_prime_method')
full_seq_method_iri = rdflib.URIRef('http://purl.example.com/ibanez/full_sequence_method')

# Fake Named Graph URIs used to:
# - select which tissue or ibanez homology method to use in queries
# - organize loading and reloading of data into the database.
affymetrix_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_affymetrix')
flybase_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_flybase')
microcosm_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_microcosm')
mirbase_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_mirbase')
roundup_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_roundup')
synaptomedb_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_synaptomedb')
targetscan_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_targetscan')
uniprot_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_uniprot')
ibanez_five_prime_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_ibanez_five_prime')
ibanez_full_sequence_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_ibanez_full_sequence')
mcneill_cns_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_mcneill_cns')
mcneill_muscle_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_mcneill_muscle')

# Fake Classes
mature_mirna_cl = rdflib.URIRef('http://purl.mirbase.org/owl/mature_mirna')
# A real class for mature microRNAs exists in Sequence Ontology.  I am not
# using it because it is not readable (in a SPARQL query).  I did ask on the
# SO mailing list and OBO discuss list and Andy from flybase says that yes
# they use it to annotate mirbase mature mirnas.
# http://purl.obolibrary.org/obo/SO_0000276

# Real Predicates
orthologous_to_pred = rdflib.URIRef('http://purl.org/obo/owl/obo#orthologous_to')
homologous_to_pred = rdflib.URIRef('http://purl.org/obo/owl/obo#homologous_to')
organism_pred = rdflib.URIRef('http://purl.uniprot.org/core/organism')
type_pred = rdflib.URIRef('http://www.w3.org/1999/02/22-rdf-syntax-ns#type')
see_also_pred = rdflib.URIRef('http://www.w3.org/2000/01/rdf-schema#seeAlso')
classified_pred = rdflib.URIRef('http://purl.uniprot.org/core/classifiedWith')
db_pred = rdflib.URIRef('http://purl.uniprot.org/core/database')
# Fake Predicates
is_mature_pred = rdflib.URIRef('http://purl.mirbase.org/owl/is_mature')
current_id_pred = rdflib.URIRef('http://purl.mirbase.org/owl/current_id')
has_id_pred = rdflib.URIRef('http://purl.mirbase.org/owl/has_id')
targets_pred = rdflib.URIRef('http://purl.mirbase.org/owl/targets')
regulates_pred = rdflib.URIRef('http://purl.example.com/owl/regulates') # a regulates the gene expression of a gene.

# Fake Database IRIs.  Databases represent the origin of a sequence/entity, the
# scope of an id.  There is probably a better way to represent this, since this
# seems to conflate the organization and id/sequence type.
targetscan_mir_family_db = rdflib.URIRef('http://purl.example.com/database/targetscan_mir_family')
refseq_db = rdflib.URIRef('http://purl.example.com/database/ncbi_refseq')
flybase_annotation_db = rdflib.URIRef('http://purl.example.com/database/flybase_annotation')
flybase_gene_db = rdflib.URIRef('http://purl.example.com/database/flybase_gene')
flybase_transcript_db = rdflib.URIRef('http://purl.example.com/database/flybase_transcript')
mirbase_acc_db = rdflib.URIRef('http://purl.example.com/database/mirbase_acc')
mirbase_id_db = rdflib.URIRef('http://purl.example.com/database/mirbase_id')
ensembl_gene_db = rdflib.URIRef('http://purl.example.com/database/ensembl_gene')
ensembl_transcript_db = rdflib.URIRef('http://purl.example.com/database/ensembl_transcript')
affymetrix_probeset_db = rdflib.URIRef('http://purl.example.com/database/affymetrix_probeset')
uniprot_db = rdflib.URIRef('http://purl.example.com/database/uniprot')


def targetscan_mir_family_iri(family):
    return rdflib.URIRef('http://purl.targetscan.org/mir_family/{}'.format(urllib.quote(family, safe='')))

def refseq_iri(refseq):
    return rdflib.URIRef('http://purl.ncbi.nlm.nih.gov/refseq/{}'.format(refseq ))

def flybase_annotation_iri(annotation):
    return rdflib.URIRef('http://purl.flybase.org/annotation_id/{}'.format(annotation))

def flybase_gene_iri(gene):
    return rdflib.URIRef('http://purl.flybase.org/gene_id/{}'.format(gene))

def flybase_transcript_iri(transcript):
    return rdflib.URIRef('http://purl.flybase.org/transcript_id/{}'.format(transcript))

def mirbase_acc_iri(acc):
    return rdflib.URIRef('http://purl.mirbase.org/mirna_acc/{}'.format(acc))

def mirbase_id_iri(mirna_id):
    return rdflib.URIRef('http://purl.mirbase.org/mirna_id/{}'.format(mirna_id))

def ensembl_gene_iri(gene):
    return rdflib.URIRef('http://purl.ensembl.org/gene/{}'.format(gene))

def ensembl_transcript_iri(transcript):
    return rdflib.URIRef('http://purl.ensembl.org/transcript/{}'.format(transcript))

def uniprot_iri(acc):
    return rdflib.URIRef('http://purl.uniprot.org/uniprot/{}'.format(acc))

def affymetrix_probeset_iri(probeset):
    return rdflib.URIRef('http://purl.affymetrix.com/probeset/{}'.format(probeset))


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


def makefiledirs(path, mode=0775):
    '''
    Idempotent function for making any non-existent directories for the
    file path `path`.  Returns path.
    '''
    makedirs(os.path.dirname(path), mode=mode)
    return path

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


##############################
# RDF GRAPH DATABASE FUNCTIONS


def sparql_json_query(query):
    print 'sparql_json_query'
    print query
    return sparq.query(query, accept='application/sparql-results+json')


def query_for_ids(query, binding):
    '''
    Wrap a sparql query into a call to stardog, and parse out the id value of 
    the results for the given binding.  The binding values should be IRIs of
    the form "http://example.com/path/to/id".  Return a list of ids.
    '''
    result = sparql_json_query(query)
    # IRIs are like u'http://purl.targetscan.org/mir_family/miR-33'
    # or u'http://purl.targetscan.org/mir_family/miR-279%2F286%2F996'
    iris = [b[binding]['value'] for b in result['results']['bindings']]
    return iris_to_ids(iris)



def iris_to_ids(iris, token='/'):
    '''
    Convert a list of URIs/IRIs into plain old identifiers.  The ids are also
    URL unquoted, which should be safe, considering that, as IRIs, they should
    be quoted, right?  Or is that just a URL thing?
    iris: a list of IRIs in the form: 'http://example.com/path/to/identifier'
    token: a character (or string) to split the identifier from the rest of
    the IRI.  In RDF, it is typically a slash ('/') or a hash ('#').
    '''
    return [urllib.unquote(iri.rsplit(token, 1)[-1]) for iri in iris]


def prefixes():
    '''
    Return a string containing several common PREFIX clauses used in the
    mirna rdf database.
    '''
    return '''
    PREFIX syndb:<http://purl.synaptomedb.org/owl/>
    PREFIX mb:<http://purl.mirbase.org/owl/>
    PREFIX ts:<http://purl.targetscan.org/owl/>
    PREFIX ex:<http://purl.example.com/owl/>
    PREFIX obo:<http://purl.org/obo/owl/obo#>
    PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
    PREFIX up:<http://purl.uniprot.org/core/>
    PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
    PREFIX db:<http://purl.example.com/database/>
    '''


def stardog_database_name():

    return 'mirna'


def drop_stardog_database():
    cmd = 'time stardog-admin db drop --username {} --passwd {}'.format(
        secrets.stardog_admin_user, secrets.stardog_admin_password)
    cmd += ' ' + stardog_database_name()
    call(cmd)


def load_stardog_database():
    '''
    Data is loaded when the database is created b/c loading is much faster
    that way, according to the stardog docs.
    '''
    rdf_paths = all_rdf_paths()
    cmd = 'time stardog-admin db create --username {} --passwd {}'.format(
        secrets.stardog_admin_user, secrets.stardog_admin_password)
    cmd += ' --name {} {}'.format(stardog_database_name(), ' '.join(rdf_paths))
    call(cmd)


def all_rdf_files_and_graphs():
    return [
        (affymetrix_fly_annotations_rdf_path(), affymetrix_ng),
        (flybase_rdf_path(), flybase_ng),
        (ibanez_5prime_mir_homologs_rdf_path(), ibanez_five_prime_ng),
        (ibanez_fullseq_mir_homologs_rdf_path(), ibanez_full_sequence_ng),
        (mcneill_cns_screen_rdf_path(), mcneill_cns_ng),
        (mcneill_muscle_screen_rdf_path(), mcneill_muscle_ng),
        (microcosm_rdf_path(), microcosm_ng),
        (mirbase_rdf_path(), mirbase_ng),
        (roundup_rdf_path(), roundup_ng),
        (synaptomedb_v1_rdf_path(), synaptomedb_ng),
        (targetscan_rdf_path(), targetscan_ng),
        (uniprot_rdf_path(), uniprot_ng),
    ]


def all_rdf_paths():
    return [filename for filename, graph in all_rdf_files_and_graphs()]


def load_virtuoso_database():
    for filename, graph in all_rdf_files_and_graphs():
        # drop any existing graph
        if sparq.exists_graph(graph):
            print 'dropping graph {}'.format(graph)
            sparq.drop_graph(graph)
        # load files into graph, one at a time
        print 'loading graph {}'.format(graph)
        virt7.load_graph(graph, filename)


################
# SPARQL QUERIES

def update_mirbase_ids(mirbase_ids):
    '''
    Take a list of mirbase ids for MATURE mirbase accessions, and return a dict
    that maps each id to the current mirbase id for the accession the mirbase
    id maps to.  Possible data quirks.  If the mirbase id does not map, it will
    not be in the returned dict.  If the mirbase id maps to multiple mature
    mirbase accs (I do not *think* this happens), the returned dict will only
    map the original to one of the current mirbase ids for those mirbase accs.
    Finally, if the original mirbase id could be the same as the current
    mirbase id if it is up-to-date.
    '''
    # Generate URIs for the mirbase ids
    mids = [mirbase_id_iri(mirbase_id) for mirbase_id in mirbase_ids + ['foobar']]
    # Query for a mapping from original id to current id
    query = prefixes() + '''
    SELECT DISTINCT ?dm_old ?dm
    WHERE {
    VALUES ?dm_old { ''' + ' '.join(['<{}>'.format(mid) for mid in mids]) + ''' }
    # convert old mirbase ids to mirbase accs
    ?dma mb:has_id ?dm_old .
    ?dm_old up:database db:mirbase_id .
    # convert mirbase accs to current mirbase ids
    ?dma up:database db:mirbase_acc .
    ?dma a mb:mature_mirna .
    ?dma mb:current_id ?dm .
    }
    "
    '''
    result = sparql_json_query(query)
    # iris are like u'http://purl.targetscan.org/mir_family/miR-33'
    # or u'http://purl.targetscan.org/mir_family/miR-279%2F286%2F996'
    lookup = dict((os.path.basename(b['dm_old']['value']), 
                   os.path.basename(b['dm']['value'])) for b in
                  result['results']['bindings'])
    return lookup # map old id to new id


def update_validated_mirs_to_current_mirbase_ids():
    lookup = update_mirbase_ids(original_validated_mirs)
    for mir in original_validated_mirs:
        print mir, lookup.get(mir)

    current_mirs = [lookup[mir] for mir in original_validated_mirs if mir in lookup]
    print current_mirs
    print len(current_mirs)
    return current_mirs


def targetscan_human_mirs_targeting_conserved_synapse_genes():
    query = prefixes() + '''
    SELECT DISTINCT ?hm
    WHERE {
    # human synaptic genes
    ?hg up:classifiedWith syndb:synaptic .
    ?hg up:organism taxon:9606 .
    ?hu rdfs:seeAlso ?hg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .

    # human refseq transcripts
    ?hu rdfs:seeAlso ?ht .
    ?ht up:database db:ncbi_refseq .
    ?ht up:organism taxon:9606 .

    # targetscan mirs targeting human refseq transcripts
    ?hmf mb:targets ?ht .
    ?hmf up:database db:targetscan_mir_family .

    # convert targetscan mir family into current mirbase accs
    ?hmf rdfs:seeAlso ?hma .
    ?hma up:organism taxon:9606 .

    # convert mirbase accs into current mirbase ids
    ?hma up:database db:mirbase_acc .
    ?hma a mb:mature_mirna .
    ?hma mb:current_id ?hm .
    }
    '''
    return query_for_ids(query, 'hm')

def targetscan_fly_mirs_targeting_conserved_synapse_genes():
    '''
    fly targetscan family targeting flybase gene linked to uniprot orthologous
    to human uniprot linked to human ensembl gene classified with synaptic.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?dm
    WHERE {
    # human synaptic ensembl genes
    ?hg up:classifiedWith syndb:synaptic .
    ?hg up:organism taxon:9606 .
    ?hg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?hg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .

    # fly genes
    ?du rdfs:seeAlso ?dg .
    ?dg up:database db:flybase_gene .
    # targetscan mir family targeting fly genes
    ?dmf mb:targets ?dg .
    ?dmf up:database db:targetscan_mir_family .
    # convert targetscan mir family to mirbase accs
    ?dmf rdfs:seeAlso ?dma .
    ?dma up:organism taxon:7227 .
    # convert mirbase accs to current mirbase ids
    ?dma up:database db:mirbase_acc .
    ?dma a mb:mature_mirna .
    ?dma mb:current_id ?dm .
    }
    '''
    return query_for_ids(query, 'dm')


def microcosm_human_mirs_targeting_conserved_synapse_genes():
    '''
    Return a list of human mirs that are predicted to target genes that are
    synaptic genes (according the synapsedb) and orthologous to fly genes.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?hm
    WHERE {

    # synaptic human ensembl genes
    ?hg up:classifiedWith syndb:synaptic .
    ?hg up:organism taxon:9606 .
    ?hg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?hg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .

    # human ensembl transcripts
    ?hu rdfs:seeAlso ?ht .
    ?ht up:database db:ensembl_transcript.
    ?ht up:organism taxon:9606 .

    # microcosm mirbase ids targeting transcripts
    ?hm_old mb:targets ?ht .

    # convert microcosm mirbase ids into mirbase accs
    ?hma mb:has_id ?hm_old .
    ?hma up:organism taxon:9606 .

    # convert mirbase accs into current mirbase ids
    ?hma up:database db:mirbase_acc .
    ?hma a mb:mature_mirna .
    ?hma mb:current_id ?hm .
    }
    '''
    return query_for_ids(query, 'hm')



def microcosm_fly_mirs_targeting_conserved_synapse_genes():
    '''
    fly targetscan family targeting flybase gene linked to uniprot orthologous
    to human uniprot linked to human ensembl gene classified with synaptic.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?dm
    WHERE {
    # human synaptic ensembl genes
    ?hg up:classifiedWith syndb:synaptic .
    ?hg up:organism taxon:9606 .
    ?hg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?hg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .

    # fly transcrips
    ?du rdfs:seeAlso ?dt .
    ?dt up:database db:flybase_transcript .

    # fly annotations
    ?dt rdfs:seeAlso ?da .
    ?da up:database db:flybase_annotation .

    # microcosm mirbase ids targeting fly annotation ids
    ?dm_old mb:targets ?da .

    # convert old mirbase ids to mirbase accs
    ?dma mb:has_id ?dm_old .
    ?dma up:database db:mirbase_acc .
    ?dma a mb:mature_mirna .

    # convert mirbase accs to current mirbase ids
    ?dma mb:current_id ?dm .
    ?dm up:database db:mirbase_id .
    }
    '''
    return query_for_ids(query, 'dm')


def mircocosm_predicted_human_mir_targets(mirbase_id):
    '''
    Given a mirbase id, go to mature mirbase accession, to mirbase ids
    associated with that accession, to ensembl transcript targets to uniprot
    to ensembl genes.
    '''
    mir = mirbase_id_iri(mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?gene
    WHERE {{
    # original mirbase id to mature mirbase accession
    ?macc mb:has_id <{mid_orig}> .
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .

    # mirbase acc to all mirbase ids
    ?macc mb:has_id ?mid .
    ?mid up:database db:mirbase_id .

    # mirbase id to targeted ensembl transcripts
    ?mid mb:targets ?enst .
    ?enst up:database db:ensembl_transcript .

    # ensembl transcript to uniprot
    ?u rdfs:seeAlso ?enst .
    ?u up:database db:uniprot .

    # uniprot to human ensembl gene .
    ?u rdfs:seeAlso ?gene .
    ?gene up:database db:ensembl_gene .
    ?gene up:organism taxon:9606 .
    }}
    '''.format(mid_orig=mir)
    return query_for_ids(query, 'gene')


def mircocosm_predicted_fly_mir_targets(mirbase_id):
    '''
    Given a mirbase id, go to a mature mirbase accession and back to mirbase
    ids (since we are not sure if the input mirbase_id is the same as one of
    the mirbase_ids used in microcosm), then to flybase accession targets and eventually to flybase genes.
    '''
    mir = mirbase_id_iri(mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?gene
    WHERE {{
    # original mirbase id to mature mirbase accession
    ?macc mb:has_id <{mid_orig}> .
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .

    # mirbase acc to all mirbase ids
    ?macc mb:has_id ?mid .
    ?mid up:database db:mirbase_id .

    # mirbase id to targeted flybase annotation ids
    ?mid mb:targets ?fba .
    ?fba up:database db:flybase_annotation .
    ?fba up:organism taxon:7227 .

    # flybase annotation id to flybase transcript
    ?fbt rdfs:seeAlso ?fba .
    ?fbt up:database db:flybase_transcript .

    # flybase transcript to uniprot
    ?u rdfs:seeAlso ?fbt .
    ?u up:database db:uniprot .

    # uniprot to flybase gene .
    ?u rdfs:seeAlso ?gene .
    ?gene up:database db:flybase_gene .
    ?gene up:organism taxon:7227 .
    }}
    '''.format(mid_orig=mir)
    return query_for_ids(query, 'gene')


def targetscan_predicted_human_mir_targets(mirbase_id):
    '''
    Given a mirbase id, go to a mature mirbase accession, to a targetscan mir family,
    to predicted targeted refseq to uniprot to ensembl gene.
    '''
    mir = mirbase_id_iri(mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?gene
    WHERE {{
    # original mirbase id to mature mirbase accession
    ?macc mb:has_id <{mid_orig}> .
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .

    # mirbase acc to targetscan mir family
    ?tmf rdfs:seeAlso ?macc .
    ?tmf up:database db:targetscan_mir_family .

    # targetscan mir family to targeted refseq transcripts
    ?tmf mb:targets ?ref .
    ?ref up:database db:ncbi_refseq .

    # refseq to uniprot
    ?u rdfs:seeAlso ?ref .
    ?u up:database db:uniprot .

    # uniprot to human ensembl gene .
    ?u rdfs:seeAlso ?gene .
    ?gene up:database db:ensembl_gene .
    ?gene up:organism taxon:9606 .
    }}
    '''.format(mid_orig=mir)
    return query_for_ids(query, 'gene')


def targetscan_predicted_fly_mir_targets(mirbase_id):
    '''
    Given a mirbase id, go to a mature mirbase accession, to a targetscan mir family,
    to predicted targeted flybase genes.
    '''
    mir = mirbase_id_iri(mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?gene
    WHERE {{
    # original mirbase id to mature mirbase accession
    ?macc mb:has_id <{mid_orig}> .
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .

    # mirbase acc to targetscan mir family
    ?tmf rdfs:seeAlso ?macc .
    ?tmf up:database db:targetscan_mir_family .

    # targetscan mir family to targeted flybase genes
    ?tmf mb:targets ?gene .
    ?gene up:database db:flybase_gene .
    ?gene up:organism taxon:7227 .
    }}
    '''.format(mid_orig=mir)
    return query_for_ids(query, 'gene')


def screened_fly_mir_targets_in_tissue(mirbase_id, tissue):

    # mcneill screen edges
    # graph.add((mir, regulates_pred, probeset))
    # graph.add((mir, db_pred, mirbase_id_db))
    # graph.add((mir, organism_pred, dro_iri))
    # graph.add((probeset, db_pred, affymetrix_probeset_db))
    # graph.add((probeset, organism_pred, dro_iri))

    # affy edges
    # graph.add((probeset, see_also_pred, gene))
    # graph.add((probeset, db_pred, affymetrix_probeset_db))
    # graph.add((probeset, organism_pred, dro_iri))
    # graph.add((gene, db_pred, flybase_gene_db))
    # graph.add((gene, organism_pred, dro_iri))

    if tissue == muscle_tissue:
        tissue_ng = mcneill_muscle_ng
    elif tissue == cns_tissue:
        tissue_ng = mcneill_cns_ng
    else:
        raise Exception('Unknown tissue {}. mirbase_id {}.'.format(
            tissue, mirbase_id))

    mir = mirbase_id_iri(mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?gene
    WHERE {{
    # original mirbase id to mature mirbase accession
    ?macc mb:has_id <{mid_orig}> .
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .

    # mirbase acc to all mirbase ids
    ?macc mb:has_id ?mid .
    ?mid up:database db:mirbase_id .

    # mirbase id regulates an affymetrix probe
    GRAPH <{tissue_ng}> {{
        ?mid ex:regulates ?affy .
        ?affy up:database db:affymetrix_probeset .
    }}

    # affy probe to flybase gene
    ?affy rdfs:seeAlso ?gene .
    ?gene up:database db:flybase_gene .
    ?gene up:organism taxon:7227 .
    }}
    '''.format(tissue_ng=tissue_ng, mid_orig=mir)
    return query_for_ids(query, 'gene')


####################
# WRITE OUTPUT FILES




# PHASE I
# HUMAN and FLY, MICROCOSM and TARGETSCAN, mirs targeting conserved synapse genes.
# ALSO VALIDATED FLY MIRS

def write_targetscan_human_mirs_targeting_conserved_synapse_genes():
    mirs = sorted(targetscan_human_mirs_targeting_conserved_synapse_genes())
    fn = os.path.join(config.datadir, 'predicted_targetscan_human_mirs_targeting_conserved_synapse_genes.txt')
    with open(fn, 'w') as fh:
        fh.write(''.join([mir + '\n' for mir in mirs]))


def write_microcosm_human_mirs_targeting_conserved_synapse_genes():
    mirs = sorted(microcosm_human_mirs_targeting_conserved_synapse_genes())
    fn = os.path.join(config.datadir, 'predicted_microcosm_human_mirs_targeting_conserved_synapse_genes.txt')
    with open(fn, 'w') as fh:
        fh.write(''.join([mir + '\n' for mir in mirs]))


def write_overlap_between_microcosm_human_mirs_and_targetscan_human_mirs():
    fn = os.path.join(config.datadir, 'overlap_between_microcosm_and_targetscan_human_mirs.csv')
    micro = set(microcosm_human_mirs_targeting_conserved_synapse_genes())
    targ = set(targetscan_human_mirs_targeting_conserved_synapse_genes())
    write_set_overlap_file(micro, targ, 'microcosm', 'targetscan', fn)


def write_targetscan_fly_mirs_targeting_conserved_synapse_genes():
    fn = os.path.join(config.datadir, 'predicted_targetscan_fly_mirs_targeting_conserved_synapse_genes.txt')
    mirs = sorted(targetscan_fly_mirs_targeting_conserved_synapse_genes())
    with open(fn, 'w') as fh:
        fh.write(''.join([mir + '\n' for mir in mirs]))


def write_microcosm_fly_mirs_targeting_conserved_synapse_genes():
    '''
    Write a list of fly mirs that are predicted to target genes that are
    orthologous to human synaptic genes (according the synapsedb).
    '''
    fn = os.path.join(config.datadir, 'predicted_microcosm_fly_mirs_targeting_conserved_synapse_genes.txt')
    mirs = sorted(microcosm_fly_mirs_targeting_conserved_synapse_genes())
    with open(fn, 'w') as fh:
        fh.write(''.join([mir + '\n' for mir in mirs]))


def write_overlap_between_validated_fly_mirs_and_microcosm_predicted_fly_mirs():
    fn = os.path.join(config.datadir, 'overlap_between_validated_and_predicted_microcosm_fly_mirs.csv')
    predicted = set(microcosm_fly_mirs_targeting_conserved_synapse_genes())
    validated = set(validated_mirs)
    write_set_overlap_file(validated, predicted, 'validated', 'predicted', fn)


def write_overlap_between_validated_fly_mirs_and_targetscan_predicted_fly_mirs():
    fn = os.path.join(config.datadir, 'overlap_between_validated_and_predicted_targetscan_fly_mirs.csv')
    predicted = set(targetscan_fly_mirs_targeting_conserved_synapse_genes())
    validated = set(validated_mirs)
    write_set_overlap_file(validated, predicted, 'validated', 'predicted', fn)


def write_overlap_between_microcosm_fly_mirs_and_targetscan_fly_mirs():
    fn = os.path.join(config.datadir, 'overlap_between_microcosm_and_targetscan_fly_mirs.csv')
    micro = set(microcosm_fly_mirs_targeting_conserved_synapse_genes())
    targ = set(targetscan_fly_mirs_targeting_conserved_synapse_genes())
    write_set_overlap_file(micro, targ, 'microcosm', 'targetscan', fn)


# PHASE II
# SETS and OVERLAPS for
# For EACH SCREENED MIR
#     MICROCOSM TARGETED FLY GENES (FBgn)
#     TARGETSCAN TARGETED FLY GENES (FBgn)
#     MICROCOSM TARGETED HUMAN GENES (ENSG)
#     TARGETSCAN TARGETED HUMAN GENES (ENSG)
#     For EACH TISSUE TYPE
#         SCREENED REGULATED FLY GENES (FBgn)
# - TARGETS OF 7 SCREENED MIRS in 2 TISSUES VS PREDICTED TARGETS of MICROCOSM and TARGETSCAN

def write_overlap_between_screened_and_microcosm_and_targetscan_fly_mir_targets():
    for mir in screened_mirs:
        microcosm_targets = set(mircocosm_predicted_fly_mir_targets(mir))
        targetscan_targets = set(targetscan_predicted_fly_mir_targets(mir))
        mt_fn = os.path.join(config.datadir, '{}_overlap_between_microcosm_and_targetscan_targets.csv'.format(mir))
        write_set_overlap_file(microcosm_targets, targetscan_targets, 'microcosm', 'targetscan', mt_fn)
        for tissue in TISSUES:
            screened_targets = set(screened_fly_mir_targets_in_tissue(mir, tissue))
            sm_fn = os.path.join(config.datadir, '{}_{}_overlap_between_screened_and_microcosm_targets.csv'.format(mir, tissue))
            write_set_overlap_file(screened_targets, microcosm_targets, 'screened', 'microcosm', sm_fn)
            st_fn = os.path.join(config.datadir, '{}_{}_overlap_between_screened_and_targetscan_targets.csv'.format(mir, tissue))
            write_set_overlap_file(screened_targets, targetscan_targets, 'screened', 'targetscan', st_fn)


def write_overlap_of_screened_fly_mir_targets_and_nmj_rnai_genes():
    nmj_rnai = set(get_nmj_rnai_genes())
    for mir in screened_mirs:
        for tissue in TISSUES:
            screened_targets = set(screened_fly_mir_targets_in_tissue(mir, tissue))
            fn = os.path.join(config.datadir, 'overlap_of_nmj_rnai_genes_and_screened_fly_targets',
                              '{}_{}_overlap_between_nmj_rnai_genes_and_screened_fly_targets.csv'.format(mir, tissue))
            write_set_overlap_file(nmj_rnai, screened_targets, 'nmj_rnai_genes', 'screened', fn)


def write_overlap_of_validated_fly_mir_targets_and_nmj_rnai_genes():
    nmj_rnai = set(get_nmj_rnai_genes())
    for mir in validated_mirs:
        microcosm_targets = set(mircocosm_predicted_fly_mir_targets(mir))
        fn = os.path.join(config.datadir, 'overlap_of_nmj_rnai_genes_and_predicted_fly_targets',
                          '{}_overlap_of_nmj_rnai_genes_and_microcosm_fly_targets.csv'.format(mir))
        write_set_overlap_file(nmj_rnai, microcosm_targets, 'nmj_rnai_genes', 'microcosm', fn)
        targetscan_targets = set(targetscan_predicted_fly_mir_targets(mir))
        fn = os.path.join(config.datadir, 'overlap_of_nmj_rnai_genes_and_predicted_fly_targets',
                          '{}_overlap_of_nmj_rnai_genes_and_targetscan_fly_targets.csv'.format(mir))
        write_set_overlap_file(nmj_rnai, targetscan_targets, 'nmj_rnai_genes', 'targetscan', fn)


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


###############################
# IBANEZ FLY-HUMAN MIR HOMOLOGS


def ibanez_5prime_mir_homologs_rdf_path():
    return os.path.join(config.datadir, 'ibanez', '2008', 'ibanez-2008-5prime-mir-homologs.nt')


def ibanez_fullseq_mir_homologs_rdf_path():
    return os.path.join(config.datadir, 'ibanez', '2008', 'ibanez-2008-fullseq-mir-homologs.nt')


def write_ibanez_5prime_mir_homologs_rdf():
    print 'write_ibanez_5prime_mir_homologs_rdf'
    write_ibanez_method_mir_homologs_rdf(FIVE_PRIME, ibanez_5prime_mir_homologs_rdf_path())


def write_ibanez_fullseq_mir_homologs_rdf():
    print 'write_ibanez_fullseq_mir_homologs_rdf'
    write_ibanez_method_mir_homologs_rdf(SEVENTY_PERCENT, ibanez_fullseq_mir_homologs_rdf_path())


def write_ibanez_method_mir_homologs_rdf(method, filename):
    '''
    Write the homologs for a specific method to a file as n-triples rdf.
    method: '5_prime_sequence_homolog' or '70_percent_full_sequence_homolog'
    '''
    print 'Populating graph for Ibanez 2008 miR fly-human homologs for method {}'.format(method)
    graph = rdflib.Graph()
    fly_mirs = set()
    human_mirs = set()

    for fly_mir_id, human_mir_id, meth in gen_ibanez_fly_human_homologs():
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


def download_ibanez_fly_human_homologs():
    raise NotImplementedError('''
Since there is no data file on the internet available for download, the Ibanez
data was hand-copied from the paper into a CSV file.  Upload it from my laptop
to the server location $datadir/ibanez/2008/2008_Ibanez-Ventoso_drosophila_melanogaster_homo_sapiens_mir_orthologs.csv
''')


def gen_ibanez_fly_human_homologs():
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
                yield current_dro_mir, full_homo_mir, SEVENTY_PERCENT


def print_ibanez_fly_human_homologs():
    '''
    See what the generated tuples of human-fly mir homologs look like to make
    sure they do not look wrong.
    '''
    for dro_mir, homo_mir, method in gen_ibanez_fly_human_homologs():
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


def roundup_rdf_path(divergence='0.8', evalue='1e-5', version='4'):
    d = os.path.join(config.datadir, 'roundup', 'v' + version)
    return os.path.join(d, 'roundup-{}-orthologs_{}_{}.nt'.format(
        version, divergence, evalue))


def download_roundup_orthologs():
    raise NotImplementedError('''Go to http://roundup.hms.harvard.edu/download.
Download Homo sapiens (9606) and Drosophila melanogaster (7227) 
with a 0.8 divergence and 1e-5 evalue thresholds.
This should be implemented by making an http request to the right download
url but the Orchestra filesystem is melting down right now (2013/07/07).
''')


def write_roundup_orthologs_rdf():
    '''
    For the input file for the given query taxon, subject taxon, divergence,
    evalue, and roundup version, write out a triple for each query id, subject
    id pair, the symmetrical triple, and then two triples to say that the query
    id is from the query_taxon organism and the subject id is from the subject
    taxon too.

    Note the loss of divergence, evalue, roundup version, and ortholog distance
    score.
    '''
    print 'write_roundup_orthologs_rdf'
    # parameters
    version = '4'
    divergence = '0.8'
    evalue = '1e-5'

    done_ids = set()
    graph = rdflib.Graph()
    for qid, sid, distance in gen_orthologs(dro_taxon, homo_taxon,
                                            divergence, evalue, version):
        q = uniprot_iri(qid)
        s = uniprot_iri(sid)
        graph.add((q, orthologous_to_pred, s)) # link query to subject
        graph.add((s, orthologous_to_pred, q)) # link query to subject
        if qid not in done_ids:
            graph.add((q, organism_pred, dro_iri)) # link query to taxon
            done_ids.add(qid)
        if sid not in done_ids:
            graph.add((s, organism_pred, homo_iri)) # link subject to taxon
            done_ids.add(sid)

    # n-triples file extension
    with open(roundup_rdf_path(), 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))


#########
# MIRBASE


def mirbase_aliases_url():
    return 'ftp://mirbase.org/pub/mirbase/19/aliases.txt.gz'


def mirbase_aliases_dest():
    return os.path.join(config.datadir, 'mirbase', '19',
                        os.path.basename(mirbase_aliases_url()))


def mirbase_aliases_filename():
    return mirbase_aliases_dest().rstrip('.gz')


def mirbase_rdf_path():
    return os.path.join(config.datadir, 'mirbase', '19', 'mirbase.nt')


def download_mirbase_aliases():

    url = mirbase_aliases_url()
    dest = mirbase_aliases_dest()
    download(url, dest)
    call('gunzip {}'.format(dest))


def write_mirbase_rdf():
    print 'write_mirbase_rdf'
    graph = rdflib.Graph()
    human_accs = set()
    human_ids = set()
    fly_accs = set()
    fly_ids = set()

    mature_re = re.compile(r'MIMAT\d+$')
    # link mirbase accessions to mirbase ids and collect the accession and ids
    with open(mirbase_aliases_filename()) as fh:
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

    with open(mirbase_rdf_path(), 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))


###################
# UNIPROT FUNCTIONS



def get_uniprot_version_of_sparql_endpoint():
    # Get the uniprot version.  Right now this is scraped from the html.
    # I hope it will be a VoID property someday.
    r = requests.get("http://beta.sparql.uniprot.org")
    version = re.search(r"UniProt release (\d\d\d\d_\d\d)", r.text).group(1)
    return version


def uniprot_version_path():
    '''
    Where the "active" uniprot version/release is saved.  This tracks which
    version of UniProt we have downloaded and generated RDF for.
    '''
    return os.path.join(config.datadir, 'uniprot', 'active_version.txt')


def set_uniprot_version(version=None):
    '''
    version: if falsy, use the current version of the uniprot sparql endpoint.
    '''
    version = version or get_uniprot_version_of_sparql_endpoint()
    fn = uniprot_version_path()
    with open(fn, 'w') as fh:
        fh.write(version + '\n')


def get_uniprot_version():
    fn = uniprot_version_path()
    with open(fn) as fh:
        return fh.read().strip()


def uniprot_fly_rdf_path():
    return os.path.join(uniprot_data_dir(get_uniprot_version()),
                        'uniprot_fly.nt')


def uniprot_homo_rdf_path():
    return os.path.join(uniprot_data_dir(get_uniprot_version()),
                        'uniprot_human.nt')


def uniprot_flybase_annotation_rdf_path():
    return os.path.join(uniprot_data_dir(get_uniprot_version()),
                        'uniprot_flybase_annotation.nt')


def uniprot_flybase_gene_rdf_path():
    return os.path.join(uniprot_data_dir(get_uniprot_version()),
                        'uniprot_flybase_gene.nt')


def uniprot_flybase_transcript_rdf_path():
    return os.path.join(uniprot_data_dir(get_uniprot_version()),
                        'uniprot_flybase_transcript.nt')


def uniprot_homo_refseq_rdf_path():
    return os.path.join(uniprot_data_dir(get_uniprot_version()),
                        'uniprot_homo_refseq.nt')


def uniprot_homo_ensembl_gene_rdf_path():
    return os.path.join(uniprot_data_dir(get_uniprot_version()),
                        'uniprot_homo_ensembl_gene.nt')


def uniprot_homo_ensembl_transcript_rdf_path():
    return os.path.join(uniprot_data_dir(get_uniprot_version()),
                        'uniprot_homo_ensembl_transcript.nt')


def download_uniprot_fly_rdf():
    # Construct triples for every fly uniprot entry.
    # Say that they are proteins, from fly, and from uniprot.
    qry = '''
        PREFIX up:<http://purl.uniprot.org/core/> 
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX skos:<http://www.w3.org/2004/02/skos/core#> 
        PREFIX owl:<http://www.w3.org/2002/07/owl#> 
        PREFIX bibo:<http://purl.org/ontology/bibo/> 
        PREFIX dc:<http://purl.org/dc/terms/> 
        PREFIX xsd:<http://www.w3.org/2001/XMLSchema#> 
        PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
        PREFIX db:<http://purl.uniprot.org/database/>
        CONSTRUCT { 
          ?protein a up:Protein .
          ?protein up:organism taxon:7227 .
          ?protein up:database <''' + uniprot_db + '''> .
        }
        WHERE {
          ?protein a up:Protein .
          ?protein up:organism taxon:7227 .
        }
    '''
    fn = uniprot_fly_rdf_path()
    count = download_uniprot_xrefs(qry, fn)
    # As of 2013/07 there were 37828 fly uniprots.
    if count < 30000 * 3:
        raise Exception('Expected more triples.', count)


def download_uniprot_homo_rdf():
    # Construct triples for every homo uniprot entry.
    # Say that they are proteins, from homo, and from uniprot.
    qry = '''
        PREFIX up:<http://purl.uniprot.org/core/> 
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX skos:<http://www.w3.org/2004/02/skos/core#> 
        PREFIX owl:<http://www.w3.org/2002/07/owl#> 
        PREFIX bibo:<http://purl.org/ontology/bibo/> 
        PREFIX dc:<http://purl.org/dc/terms/> 
        PREFIX xsd:<http://www.w3.org/2001/XMLSchema#> 
        PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
        PREFIX db:<http://purl.uniprot.org/database/>
        CONSTRUCT { 
          ?protein a up:Protein .
          ?protein up:organism taxon:9606 .
          ?protein up:database <''' + uniprot_db + '''> .
        }
        WHERE {
          ?protein a up:Protein .
          ?protein up:organism taxon:9606 .
        }
    '''
    fn = uniprot_homo_rdf_path()
    count = download_uniprot_xrefs(qry, fn)
    # As of 2013/07 there were 134288 homo uniprots.
    if count < 130000 * 3:
        raise Exception('Expected more triples.', count)


def download_uniprot_flybase_annotation_rdf():
    qry = '''
        PREFIX up:<http://purl.uniprot.org/core/> 
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX skos:<http://www.w3.org/2004/02/skos/core#> 
        PREFIX owl:<http://www.w3.org/2002/07/owl#> 
        PREFIX bibo:<http://purl.org/ontology/bibo/> 
        PREFIX dc:<http://purl.org/dc/terms/> 
        PREFIX xsd:<http://www.w3.org/2001/XMLSchema#> 
        PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
        PREFIX db:<http://purl.uniprot.org/database/>
        CONSTRUCT { 
          ?protein rdfs:seeAlso ?xref .
          ?xref up:organism taxon:7227 .
          ?xref up:database <''' + flybase_annotation_db + '''> .
        }
        WHERE {
          ?protein a up:Protein .
          ?protein up:organism taxon:7227 .
          ?protein rdfs:seeAlso ?xref . 
          ?xref up:database db:UCSC .
        }
    '''
    fn = uniprot_flybase_annotation_rdf_path()
    count = download_uniprot_xrefs(qry, fn)
    # As of 2013/07 there were 14953 fly flybase annotations.
    if count < 14000 * 3:
        raise Exception('Expected more triples.', count)


def download_uniprot_flybase_gene_rdf():
    qry = '''
        PREFIX up:<http://purl.uniprot.org/core/> 
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX skos:<http://www.w3.org/2004/02/skos/core#> 
        PREFIX owl:<http://www.w3.org/2002/07/owl#> 
        PREFIX bibo:<http://purl.org/ontology/bibo/> 
        PREFIX dc:<http://purl.org/dc/terms/> 
        PREFIX xsd:<http://www.w3.org/2001/XMLSchema#> 
        PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
        PREFIX db:<http://purl.uniprot.org/database/>
        CONSTRUCT { 
          ?protein rdfs:seeAlso ?xref .
          ?xref up:organism taxon:7227 .
          ?xref up:database <''' + flybase_gene_db + '''> .
        }
        WHERE {
          ?protein a up:Protein .
          ?protein up:organism taxon:7227 .
          ?protein rdfs:seeAlso ?xref .
          ?xref up:database db:FlyBase .
        }
    '''
    fn = uniprot_flybase_gene_rdf_path()
    count = download_uniprot_xrefs(qry, fn)
    # As of 2013/07 there were 23348 fly protein + flybase gene annotations.
    # And 13k unique flybase genes in those annotation
    # Wrote 50988 triples.
    if count < 50000:
        raise Exception('Expected more triples.', count)


def download_uniprot_flybase_transcript_rdf():
    # UniProt annotates proteins with flybase transcripts (FBtr* ids) using
    # the database EnsemblMetazoa.
    qry = '''
        PREFIX up:<http://purl.uniprot.org/core/> 
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX skos:<http://www.w3.org/2004/02/skos/core#> 
        PREFIX owl:<http://www.w3.org/2002/07/owl#> 
        PREFIX bibo:<http://purl.org/ontology/bibo/> 
        PREFIX dc:<http://purl.org/dc/terms/> 
        PREFIX xsd:<http://www.w3.org/2001/XMLSchema#> 
        PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
        PREFIX db:<http://purl.uniprot.org/database/>
        CONSTRUCT { 
          ?protein rdfs:seeAlso ?xref .
          ?xref up:organism taxon:7227 .
          ?xref up:database <''' + flybase_transcript_db + '''> .
        }
        WHERE {
          ?protein a up:Protein .
          ?protein up:organism taxon:7227 .
          ?protein rdfs:seeAlso ?xref .
          ?xref up:database db:EnsemblMetazoa .
        }
    '''
    fn = uniprot_flybase_transcript_rdf_path()
    count = download_uniprot_xrefs(qry, fn)
    # As of 2013/07 there were 21819 fly flybase transcript annotations.
    # Wrote 65457 triples.
    if count < 60000:
        raise Exception('Expected more triples.', count)


def download_uniprot_homo_refseq_rdf():
    # Construct human protein refseq annotation triples.
    # The NP_* ids are URIs.  The NM_* ids are comment strings on the
    # NP_* URIs.  Weird.  Also, strip the version number off the refseq ids,
    # since the TargetScan transcript ids do not have the version number
    # suffix.
    qry = '''
        PREFIX up:<http://purl.uniprot.org/core/> 
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX skos:<http://www.w3.org/2004/02/skos/core#> 
        PREFIX owl:<http://www.w3.org/2002/07/owl#> 
        PREFIX bibo:<http://purl.org/ontology/bibo/> 
        PREFIX dc:<http://purl.org/dc/terms/> 
        PREFIX xsd:<http://www.w3.org/2001/XMLSchema#> 
        PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
        PREFIX db:<http://purl.uniprot.org/database/>
        CONSTRUCT { 
          ?protein rdfs:seeAlso ?xrefuri .
          ?protein rdfs:seeAlso ?xref2uri .
          ?xrefuri up:organism taxon:9606 .
          ?xref2uri up:organism taxon:9606 .
          ?xrefuri up:database <''' + refseq_db + '''> .
          ?xref2uri up:database <''' + refseq_db + '''> .
        }
        WHERE {
          ?protein a up:Protein .
          ?protein up:organism taxon:9606 .
          ?protein rdfs:seeAlso ?xref . 
          ?xref up:database db:RefSeq .
          ?xref rdfs:comment ?xref2 .
          BIND (URI(REPLACE(str(?xref), 
            '(\\\\d+)\\\\.\\\\d+$', '$1')) as ?xrefuri)
          BIND (URI(CONCAT('http://purl.uniprot.org/refseq/', 
            REPLACE(?xref2, '(\\\\d+)\\\\.\\\\d+$', '$1'))) as ?xref2uri)
        }
    '''
    fn = uniprot_homo_refseq_rdf_path()
    count = download_uniprot_xrefs(qry, fn)
    # As of 2013/07 there were 83574 homo refseq annotation.
    # Wrote 250652 triples.
    if count < 80000 * 3:
        raise Exception('Expected more triples.', count)


def download_uniprot_homo_ensembl_transcript_rdf():
    qry = '''
        PREFIX up:<http://purl.uniprot.org/core/> 
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX skos:<http://www.w3.org/2004/02/skos/core#> 
        PREFIX owl:<http://www.w3.org/2002/07/owl#> 
        PREFIX bibo:<http://purl.org/ontology/bibo/> 
        PREFIX dc:<http://purl.org/dc/terms/> 
        PREFIX xsd:<http://www.w3.org/2001/XMLSchema#> 
        PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
        PREFIX db:<http://purl.uniprot.org/database/>
        CONSTRUCT { 
          ?protein rdfs:seeAlso ?xref .
          ?xref up:organism taxon:9606 .
          ?xref up:database <''' + ensembl_transcript_db + '''> .
        }
        WHERE {
          ?protein a up:Protein .
          ?protein up:organism taxon:9606 .
          ?protein rdfs:seeAlso ?xref . 
          ?xref up:database db:Ensembl .
          ?xref a up:Transcript_Resource .
        }
    '''
    fn = uniprot_homo_ensembl_transcript_rdf_path()
    count = download_uniprot_xrefs(qry, fn)
    # As of 2013/07 there were 98088 homo ensembl transcript annotations.
    # Wrote 294078 triples.
    if count < 90000 * 3:
        raise Exception('Expected more triples.', count)

def download_uniprot_homo_ensembl_gene_rdf():
    qry = '''
        PREFIX up:<http://purl.uniprot.org/core/> 
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX skos:<http://www.w3.org/2004/02/skos/core#> 
        PREFIX owl:<http://www.w3.org/2002/07/owl#> 
        PREFIX bibo:<http://purl.org/ontology/bibo/> 
        PREFIX dc:<http://purl.org/dc/terms/> 
        PREFIX xsd:<http://www.w3.org/2001/XMLSchema#> 
        PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
        PREFIX db:<http://purl.uniprot.org/database/>
        CONSTRUCT { 
          ?protein rdfs:seeAlso ?xref .
          ?xref up:organism taxon:9606 .
          ?xref up:database <''' + ensembl_gene_db + '''> .
        }
        WHERE {
          ?protein a up:Protein .
          ?protein up:organism taxon:9606 .
          ?protein rdfs:seeAlso ?enst .
          ?enst up:database db:Ensembl .
          ?enst a up:Transcript_Resource .
          ?enst up:transcribedFrom ?xref .
        }
    '''
    fn = uniprot_homo_ensembl_gene_rdf_path()
    count = download_uniprot_xrefs(qry, fn)
    # As of 2013/07 there were 72427 homo ensembl gene annotations.
    # Wrote 216873 triples.  Note that this is odd, b/c there are
    # 71237 distinct pairs of ?protein and ?xref, and 22092 distinct
    # ?xrefs in the query above.  So we would expect 71237 + 2 * 22092 triples
    # (= 115421).  There are that many distinct triples, so the construct
    # query returned non-distinct triples.  Weird.
    if count < 70000 * 3:
        raise Exception('Expected more triples.', count)


def download_uniprot_xrefs(qry, filename, validate=None):
    '''
    Return the count of lines (triples) returned from the query.
    qry: a SPARQL construct query.
    filename: where to save the constructed triples (as N-Triples).
    validate: a function to validate each line.
    '''
    if validate == None:
        def validate(line):
            pass
    version = get_uniprot_version_of_sparql_endpoint()
    print 'Downloading Xrefs for UniProt version:', version
    print 'Query:', qry
    sp = sparql.Sparql('http://beta.sparql.uniprot.org')
    try:
        # Use text/plain for N-triples
        nt = sp.query(qry, accept='text/plain')
    except Exception as e:
        print dir(e)
        print e.response.text
        print e.response.headers
        raise

    print 'Writing triples to file:', filename
    count = 0
    with open(filename, 'w') as fh:
        for line in StringIO.StringIO(nt):
            count += 1
            validate(line)
            fh.write(line)

    print 'Wrote {} triples.'.format(count)
    return count


def download_uniprot_idmapping():
    raise NotImplementedError('''
The idmapping file can only be downloaded for the current uniprot release.
Since we want release 2013_04, for the sake of consistency, the uniprot
idmapping file must be uploaded from my laptop.

The idmapping file is big (2.1 GB compressed).  It would be cooler, and 
and perhaps faster to download idmappings just for human and drosophila
from the uniprot sparql endpoint, http://beta.sparql.uniprot.org/.
''')


def write_uniprot_rdf():
    '''
    '''
    print 'write_uniprot_rdf'
    # id dbs, id kind, uniprot id mapping file
    # mrna refseq ids (targetscan) (human) YES, idmapping.RefSeq_NT.dat (but need to strip version off the end)
    # flybase annotation ids (microcosm, flybase) (fly) YES/MAYBE, idmapping.UCSC.dat
    # flybase gene ids (targetscan) (fly) YES, idmapping.FlyBase.dat 
    # flybase transcript ids (flybase) (fly), YES, idmapping.EnsemblGenome_TRS.dat
    # ensembl transcript ids (microcosm) (human), YES, idmapping.Ensembl_TRS.dat
    # ensembl gene ids (synaptomedb) (human), YES, idmapping.Ensembl.dat
    xref_data = [
        (refseq_iri, 'RefSeq_NT', re.compile(r'^(NM_\d+)\.\d+$'), homo_iri, refseq_db),
        (flybase_annotation_iri, 'UCSC', flybase_annotation_id_regex, dro_iri, flybase_annotation_db),
        (flybase_gene_iri, 'FlyBase', re.compile(r'(FBgn\d+)$'), dro_iri, flybase_gene_db),
        (flybase_transcript_iri, 'EnsemblGenome_TRS', re.compile(r'(FBtr\d+)$'), dro_iri, flybase_transcript_db),
        (ensembl_transcript_iri, 'Ensembl_TRS', re.compile(r'(ENST\d+)$'), homo_iri, ensembl_transcript_db),
        (ensembl_gene_iri, 'Ensembl', re.compile(r'(ENSG\d+)$'), homo_iri, ensembl_gene_db),
    ]

    graph = rdflib.Graph()
    human_uniprots = set()
    fly_uniprots = set()
    uniprots_map = {homo_iri: human_uniprots, dro_iri: fly_uniprots}

    # id_iri makes an IRI/URI for the id.
    # id_type specifies what file has the uniprot to id mapping.
    # id_regex matches the specific ids in the file that we want, since some
    # files have data for many kinds of ids or species.  Also used to select
    # only part of an id.
    # id_taxon is the IRI for the taxon of the ids we want from the file.
    # id_database is the IRI for the database of origin of the ids (e.g
    # flybase genes use 'http://purl.example.com/database/flybase_gene').
    for id_iri, id_type, id_regex, id_taxon, id_database in xref_data:
        print 'Doing {} for {} in {}'.format(id_type, id_taxon, id_database)
        for uniprot, id_kind, mapped_id in gen_uniprot_id_mappings(id_type=id_type):
            m = id_regex.search(mapped_id)
            if not m:
                # skip ids not matching the specific type and species we want
                continue
            else:
                # for refseq ids, we trim off the version suffix.
                trimmed_id = m.group(1)

            uni = uniprot_iri(uniprot)
            mapped = id_iri(trimmed_id)
            graph.add((uni, see_also_pred, mapped))
            graph.add((mapped, db_pred, id_database))
            graph.add((mapped, organism_pred, id_taxon))
            uniprots_map[id_taxon].add(uni)

    for taxon, uniprots in uniprots_map.items():
        for uni in uniprots:
            graph.add((uni, organism_pred, taxon))
            graph.add((uni, db_pred, uniprot_db))

    with open(uniprot_rdf_path(), 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))


def uniprot_data_dir(version):
    return makedirs(os.path.join(config.datadir, 'uniprot', version))


def uniprot_rdf_path(version='2013_04'):
    return os.path.join(uniprot_data_dir(version),
                        'uniprot-{}-idmapping.nt'.format(version))


def uniprot_idmapping_path(version='2013_04', id_type=None):
    d = uniprot_data_dir(version)
    if id_type is None:
        return os.path.join(d, 'idmapping.dat')
    else:
        return os.path.join(d, 'idmapping.{}.dat'.format(id_type))


def split_uniprot_idmapping_file():
    '''
    Split the idmapping.dat file into one file for each id_type.  This 
    allows more efficient searching for ids than having everything in one
    big file.  I really should load it into a database like sqlite or
    berkeleydb.
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


def guess_uniprot_mapping_idtype(identifier):
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
# MCNEILL SCREEN

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
    write_mcneill_tissue_screen_rdf(muscle_tissue,
                                    mcneill_muscle_screen_rdf_path())


def write_mcneill_cns_screen_rdf():
    '''
    Write an rdf file for CNS tissue results from Elizabeth McNeill's
    expression array screen of 7 fly miRs.
    '''
    print 'write_mcneill_cns_screen_rdf'
    write_mcneill_tissue_screen_rdf(cns_tissue,
                                    mcneill_cns_screen_rdf_path())


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


#########
# MinoTar


def minotar_dir():
    # version is made up b/c I can not find a version on the minotar site.
    return os.path.join(config.datadir, 'minotar', 'v1')


def minotar_fly_dir():
    return os.path.join(minotar_dir(), 'drosophila_melanogaster')


def minotar_human_dir():
    return os.path.join(minotar_dir(), 'homo_sapiens')


def minotar_fly_csv_path():
    return os.path.join(minotar_fly_dir(), 'DrosConservedSeedsConservedTargs.csv')


def minotar_human_csv_path():
    return os.path.join(minotar_human_dir(), 'HumanConservedSeedsConservedTargs.csv')


def download_minotar():
    '''
    For fly:
    The columns in the target prediction files are:
        Gene ID
        Transcript ID
        Symbol
        CG identifier
        Conserved 8Mers
        Conserved 7Mers-m8
        Conserved 7Mers-1A
        Non-Conserved 8Mers
        Non-Conesrved 7Mers-m8
        Non-Conserved 7Mers-1A
        Probability Conserved Targeting

    For human:
    The columns in the target prediction files are:
        RefSeq ID
        Gene Symbol
        Conserved 8Mers
        Non-Conserved 8Mers
        Non-Conserved 7Mers-m8
        NonConserved 7Mers-1A
        Probability Conserved Targeting
    '''
    # Minotar (v1) was generated using mir families from fly_52orfs and vert_50
    # according to personal correspondence with Bonnie Berger.
    # Download the mir family info for those datasets in addition to the
    # minotar data.


    print 'Downloadding MinoTar supporting data from TargetScanVert version vert_50.'
    vert50_dir = os.path.join(targetscan_top_dir(), '5.2', 'vert_50')
    url = 'http://www.targetscan.org//vert_50/vert_50_data_download/miR_Family_Info.txt.zip'
    download_and_unzip(url, vert50_dir)

    print 'Downloading MinoTar supporting data from TargetScanFly version fly_52orfs.'
    orfs_dir = os.path.join(targetscan_top_dir(), '6.2', 'fly_52orfs')
    url = 'http://www.targetscan.org/fly_52orfs/fly_52orfs_data_download/miR_Family_Info.txt.zip'
    download_and_unzip(url, orfs_dir)
    url = 'http://www.targetscan.org/fly_52orfs/fly_52orfs_data_download/ORF_Sequences.txt.zip'
    download_and_unzip(url, orfs_dir)
    url = 'http://www.targetscan.org/fly_52orfs/fly_52orfs_data_download/Predicted_Targets_Info.txt.zip'
    download_and_unzip(url, orfs_dir)

    print 'Downloading MinoTar data files.'
    url = 'http://www.flyrnai.org/supplement/HumanUniqueMirnas.txt'
    download_and_unzip(url, minotar_human_dir())
    url = 'http://www.flyrnai.org/supplement/HumanConservedSeedsConservedTargs.zip'
    download_and_unzip(url, minotar_human_dir())
    url = 'http://www.flyrnai.org/supplement/DrosUniqueMirnas.txt'
    download_and_unzip(url, minotar_fly_dir())
    url = 'http://www.flyrnai.org/supplement/DrosConservedSeedsConservedTargs.zip'
    download_and_unzip(url, minotar_fly_dir())


def concatenate_minotar_fly():
    '''
    Concatenate the dozens of files downloaded from minotar, each of which
    encodes the targetscan mir family in the filename, into a single csv file
    where the (unencoded) mir family is prepended as the first field of each
    line.
    '''
    targetfiles_dir = os.path.join(minotar_fly_dir(), 'ConservedSeedsConservedTargs')
    def family_parser(fn):
        return os.path.splitext(os.path.basename(fn).replace('_', '/'))[0]

    concatenate_minotar_sub(targetfiles_dir, family_parser, minotar_fly_csv_path())


def concatenate_minotar_human():
    '''
    Concatenate the dozens of files downloaded from minotar, each of which
    encodes the targetscan mir family in the filename, into a single csv file
    where the (unencoded) mir family is prepended as the first field of each
    line.
    '''
    targetfiles_dir = os.path.join(minotar_human_dir(), 'ConservedSeedsConservedTargs')
    def family_parser(fn):
        return os.path.splitext(os.path.basename(fn).replace(':', '/'))[0]

    concatenate_minotar_sub(targetfiles_dir, family_parser, minotar_human_csv_path())


def concatenate_minotar_sub(files_dir, family_parser, outfile):
    with open(outfile, 'w') as fh:
        for fn in os.listdir(files_dir):
            if not fn.endswith('.txt'):
                print 'Warning: file not ending in ".txt" dir={}, file={}'.format(files_dir, fn)
                continue
            mir_family = family_parser(fn)
            print mir_family
            with open(os.path.join(files_dir, fn)) as infh:
                for line in infh:
                    if not line.strip():
                        raise Exception('Unexpected blank line', fn, line)
                    if line.strip().startswith('#'):
                        raise Exception('Unexpected comment line', fn, line)
                    fh.write('{},{}'.format(mir_family, line))






############
# TARGETSCAN


def targetscan_dir(version='6.2'):
    return os.path.join(targetscan_top_dir(), '6.2')


def targetscan_top_dir(version='6.2'):
    return os.path.join(config.datadir, 'targetscan')


def targetscan_rdf_path(version='6.2'):
    return os.path.join(targetscan_dir(version),
                        'targetscan-{}.nt'.format(version))


def targetscan_fly_symbol_to_gene():
    symbol_to_gene = {}
    for gene_id, symbol in gen_targetscan_fly_gene_to_symbol():
        if symbol in symbol_to_gene:
            raise Exception('Symbol can map to only one gene', symbol, gene_id)
        else:
            symbol_to_gene[symbol] = gene_id

    return symbol_to_gene


def write_targetscan_rdf():
    print 'write_targetscan_rdf'
    graph = rdflib.Graph()
    families = set()
    fly_mirbase_accs = set()
    fly_genes = set()
    human_mirbase_accs = set()
    human_refseqs = set()

    # link family to fly gene target
    # converting fly symbol to flybase gene id
    symbol_to_gene = targetscan_fly_symbol_to_gene()
    for family, symbol in gen_targetscan_fly_predicted_targets():
        fam = targetscan_mir_family_iri(family)
        gene = flybase_gene_iri(symbol_to_gene[symbol])
        graph.add((fam, targets_pred, gene))
        families.add(fam)
        fly_genes.add(gene)

    # link mir_family to fly mirbase_acc
    for family, mirbase_id, mirbase_acc in gen_targetscan_fly_mir_family_info():
        fam = targetscan_mir_family_iri(family)
        acc = mirbase_acc_iri(mirbase_acc)
        graph.add((fam, see_also_pred, acc))
        families.add(fam)
        fly_mirbase_accs.add(acc)

    # link family to human refseq transcript
    for family, refseq_transcript in gen_targetscan_human_predicted_targets():
        fam = targetscan_mir_family_iri(family)
        ref = refseq_iri(refseq_transcript)
        graph.add((fam, targets_pred, ref))
        families.add(fam)
        human_refseqs.add(ref)

    # link mir_family to human mirbase_acc
    for family, mirbase_id, mirbase_acc in gen_targetscan_human_mir_family_info():
        fam = targetscan_mir_family_iri(family)
        acc = mirbase_acc_iri(mirbase_acc)
        graph.add((fam, see_also_pred, acc))
        families.add(fam)
        human_mirbase_accs.add(acc)

    # indicate that families are from the targetscan db
    for fam in families:
        graph.add((fam, db_pred, targetscan_mir_family_db))

    # fly genes are from the flybase database and drosophila organism
    for gene in fly_genes:
        graph.add((gene, db_pred, flybase_gene_db))
        graph.add((gene, organism_pred, dro_iri))

    # human refseqs are from refseq database and human organism
    for ref in human_refseqs:
        graph.add((ref, db_pred, refseq_db))
        graph.add((ref, organism_pred, homo_iri))

    # fly mirbase accs are from the mirbase database and drosophila organism
    for acc in fly_mirbase_accs:
        graph.add((acc, db_pred, mirbase_acc_db))
        graph.add((acc, organism_pred, dro_iri))

    # human mirbase accs are from mirbase database and human organism
    for acc in human_mirbase_accs:
        graph.add((acc, db_pred, mirbase_acc_db))
        graph.add((acc, organism_pred, homo_iri))

    with open(targetscan_rdf_path(), 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))


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


def flybase_rdf_path(version='FB2013_02'):
    return os.path.join(flybase_release_dir(version), 'flybase_{}.nt'.format(version))


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


def write_flybase_rdf():
    '''
    '''
    print 'write_flybase_rdf'
    graph = rdflib.Graph()
    for taxon, trans_id, anno_id in gen_flybase_transcript_mapping_file():
        transcript = flybase_transcript_iri(trans_id)
        annotation = flybase_annotation_iri(anno_id)
        graph.add((transcript, see_also_pred, annotation))
        graph.add((transcript, db_pred, flybase_transcript_db))
        graph.add((annotation, db_pred, flybase_annotation_db))

    outfile = flybase_rdf_path()
    with open(outfile, 'w') as outfh:
        outfh.write(graph.serialize(format='nt'))


def download_flybase_transcript_reporting_xml():
    raise NotImplementedError()


def flybase_transcript_reporting_xml_path():
    release_dir = flybase_release_dir()
    return os.path.join(release_dir, 'reporting-xml', 'FBtr.xml')


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
    taxon = dro_taxon
    path = flybase_transcript_reporting_xml_path()
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
    path = flybase_transcript_reporting_xml_path()
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



###############
# MAIN WORKFLOW


def cli():

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
    subparser = subparsers.add_parser('drop_stardog_database')
    subparser = subparsers.add_parser('load_stardog_database')

    subparser = subparsers.add_parser('write_targetscan_fly_mirs_targeting_conserved_synapse_genes')
    subparser = subparsers.add_parser('write_targetscan_human_mirs_targeting_conserved_synapse_genes')
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

    subparser = subparsers.add_parser('make_flybase_transcript_map', help='')

    subparser = subparsers.add_parser('map_affymetrix_fly_probe_ids_to_flybase_genes', help='')

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


def workflow():

    # Clear out a graph loaded into the wrong graph URI
    # Slow for big graphs.
    # virt.clear_graph('http://purl.roundup.hms.harvard.edu/graph/uniprot_taxonomy')

    if False: # Done tasks
        pass

        download_microcosm_targets()
        download_minotar()
        download_mirbase_aliases()
        download_targetscan()
        # Not Implemented b/c not available online or not available anymore for
        # the version we want (uniprot 2013_04) or just not done yet (Roundup,
        # Flybase)
        # download_affymetrix_fly_annotations_file()
        # download_ibanez_fly_human_homologs()
        # download_mcneill_screen()
        # download_roundup_orthologs()
        # download_synaptome_v1()
        # download_flybase_transcript_reporting_xml()

        write_affymetrix_fly_annotations_rdf()
        write_flybase_rdf()
        write_ibanez_5prime_mir_homologs_rdf()
        write_ibanez_fullseq_mir_homologs_rdf()
        write_mcneill_cns_screen_rdf()
        write_mcneill_muscle_screen_rdf()
        write_microcosm_rdf()
        write_mirbase_rdf()
        write_roundup_orthologs_rdf()
        write_synaptomedb_rdf()
        write_targetscan_rdf()

        # Get UniProt ID Mapping RDF from beta.sparql.uniprot.org.
        set_uniprot_version()
        download_uniprot_fly_rdf()
        download_uniprot_homo_rdf()
        download_uniprot_flybase_annotation_rdf()
        download_uniprot_flybase_gene_rdf()
        download_uniprot_flybase_transcript_rdf()
        download_uniprot_homo_refseq_rdf()
        download_uniprot_homo_ensembl_transcript_rdf()
        download_uniprot_homo_ensembl_gene_rdf()
        # Or get UniProt ID Mapping RDF from idmapping.dat (UniProt ftp file.)
        # write_uniprot_rdf()

        load_virtuoso_database()
    else:

        # download_uniprot_xrefs()
        # load_affymetrix_fly_annotations_rdf()
        # load_flybase_rdf()
        # load_ibanez_5prime_mir_homologs_rdf()
        # load_ibanez_fullseq_mir_homologs_rdf()
        # load_mcneill_cns_screen_rdf()
        # load_mcneill_muscle_screen_rdf()
        # load_microcosm_rdf()
        # load_mirbase_rdf()
        # load_roundup_orthologs_rdf()
        # load_synaptomedb_rdf()
        # load_targetscan_rdf()
        # load_uniprot_rdf()



        pass




def main():
    workflow()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception('')
        raise


#################
# DEPRECATED CODE


def stardog_json_query(query):
    raise DeprecationWarning()
    cmd = 'stardog --username {} --passwd {} query'.format(
        secrets.stardog_user, secrets.stardog_password)
    cmd += ' --format JSON'
    cmd += ' "mirna;reasoning=QL"'
    cmd += ' "' + query + '"'
    print 'stardog_json_query()'
    print cmd
    out = subprocess.check_output(cmd, shell=True)
    return json.loads(out)


def ibanez_mir_homologs_rdf_path():
    raise DeprecationWarning()
    return os.path.join(config.datadir, 'ibanez', '2008', 'ibanez-2008-mir-homologs.trix')


def write_ibanez_mir_homologs_rdf():
    raise DeprecationWarning()
    ds = rdflib.Dataset()

    five_prime_graph = ds.graph(ibanez_five_prime_ng)
    make_ibanez_mir_homologs_method_graph(FIVE_PRIME, five_prime_graph)

    full_seq_graph = ds.graph(ibanez_full_sequence_ng)
    make_ibanez_mir_homologs_method_graph(SEVENTY_PERCENT, full_seq_graph)

    with open(ibanez_mir_homologs_rdf_path(), 'w') as outfh:
        outfh.write(ds.serialize(format='trix'))


def make_ibanez_mir_homologs_method_graph(method, graph):
    '''
    Write the homologs for a specific method to a
    method: '5_prime_sequence_homolog' or '70_percent_full_sequence_homolog'
    '''
    raise DeprecationWarning()

    print 'Populating graph for Ibanez 2008 miR fly-human homologs for method {}'.format(method)
    fly_mirs = set()
    human_mirs = set()

    for fly_mir_id, human_mir_id, meth in gen_ibanez_fly_human_homologs():
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




