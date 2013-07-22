
'''
Contains code for creating a semantic database and querying it to answer Van Vactor
miRNA collaboration questions.
'''

import csv
import itertools
import json
import logging
import os
import subprocess
import urllib
import sys
import datetime
import functools

import rdflib

import config
import secrets
import virtuoso
# import stardog
import sparql
from mirna import mcneill, ibanez
from mirna.util import makedirs, map_flatten_set_list, call
from mirna.core import (mirbase_id_iri, ensembl_gene_iri, get_nmj_rnai_genes)
from mirna.uniprot import (uniprot_fly_rdf_path,
                           uniprot_flybase_annotation_rdf_path,
                           uniprot_flybase_gene_rdf_path,
                           uniprot_flybase_transcript_rdf_path,
                           uniprot_homo_ensembl_gene_rdf_path,
                           uniprot_homo_ensembl_transcript_rdf_path,
                           uniprot_homo_rdf_path, uniprot_homo_refseq_rdf_path,
                           download_uniprot_fly_rdf,
                           download_uniprot_flybase_annotation_rdf,
                           download_uniprot_flybase_gene_rdf,
                           download_uniprot_flybase_transcript_rdf,
                           download_uniprot_homo_ensembl_gene_rdf,
                           download_uniprot_homo_ensembl_transcript_rdf,
                           download_uniprot_homo_rdf,
                           download_uniprot_homo_refseq_rdf,
                           set_uniprot_version)
from mirna.microcosm import (download_microcosm_targets, microcosm_rdf_path,
                             write_microcosm_rdf)
from mirna.mcneill import (write_mcneill_cns_screen_rdf,
                           write_mcneill_muscle_screen_rdf,
                           mcneill_cns_screen_rdf_path,
                           mcneill_muscle_screen_rdf_path)
from mirna.mirbase import (download_mirbase_aliases, mirbase_rdf_path,
                           write_mirbase_rdf)
from mirna.targetscan import (download_targetscan, targetscan_rdf_path,
                              write_targetscan_rdf)
from mirna.affymetrix import (affymetrix_fly_annotations_rdf_path,
                              write_affymetrix_fly_annotations_rdf)
from mirna.flybase import (flybase_rdf_path, write_flybase_rdf)
from mirna.ibanez import (ibanez_five_prime_mir_homologs_rdf_path,
                          ibanez_fullseq_mir_homologs_rdf_path,
                          write_ibanez_five_prime_mir_homologs_rdf,
                          write_ibanez_fullseq_mir_homologs_rdf)
from mirna.roundup import (download_roundup_from_sparql, roundup_rdf_path,
                           roundup_sparql_rdf_path,
                           write_roundup_orthologs_rdf)
from mirna.synaptomedb import (synaptomedb_v1_rdf_path, write_synaptomedb_rdf)



VIRT7_DB = 'virt7'
STARDOG_DB = 'stardog'
# DB_TYPE = VIRT7_DB
DB_TYPE = STARDOG_DB
virt7 = virtuoso.Virtuoso7(secrets.virtuoso_dba_password,
                           config.virtuoso_load_dir)
sparq = sparql.Sparql(virt7.sparql_endpoint()) # used for querying


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

# species names used for db tables and microcosm files
mus = 'mus_musculus'
dro = 'drosophila_melanogaster'
cel = 'caenorhabditis_elegans'
homo = 'homo_sapiens'

# mirna target prediction databases
MICROCOSM = 'microcosm'
TARGETSCAN = 'targetscan'


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

muscle_tissue = mcneill.MUSCLE
cns_tissue = mcneill.CNS
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
mcneill_cns_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_mcneill_cns')
mcneill_muscle_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_mcneill_muscle')
ibanez_five_prime_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_ibanez_five_prime')
ibanez_full_sequence_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_ibanez_full_sequence')


# taxons and pairs used for roundup orthologs
taxmap = {dro: '7227', homo: '9606', mus: '10090', cel: '6239'}
TAXON_TO_NAME = {cel_taxon: cel, dro_taxon: dro, homo_taxon: homo, mus_taxon: mus}
# Orthologs between human and the other 3 species.
ROUNDUP_PAIRS = [(dro_taxon, homo_taxon), (mus_taxon, homo_taxon), (cel_taxon, homo_taxon)]

IBANEZ_METHOD_TO_NG = {ibanez.FIVE_PRIME: ibanez_five_prime_ng,
                          ibanez.FULL_SEQ: ibanez_full_sequence_ng,
                         }


##############################
# RDF GRAPH DATABASE FUNCTIONS


def stop_stardog():
    user = secrets.stardog_admin_user
    password = secrets.stardog_admin_password
    print 'To stop stardog, run the following command:'
    print 'stardog-admin server stop --username {} --passwd {}'.format(user, password)


def start_stardog():
    # In order to use Stardog in federated queries, the authentication needs to
    # be disabled for the default anonymous user.
    # http://www.stardog.com/docs/faq/.
    # Disable authentication for the `anonymous` user
    user = secrets.stardog_admin_user
    password = secrets.stardog_admin_password
    print 'To start stardog, run the following command:'
    print 'stardog-admin --disable-security server start --username {} --passwd {}'.format(user, password)


def stardog_json_query(query):
    db_name = stardog_database_name()
    # Include reasoning if inferring symmetrical edges.  Avoid it if
    # symmetrical ortholog edges were loaded.
    # db_name = '{};reasoning=QL'.format(stardog_database_name()),

    cmd = ['stardog', 'query', 'execute', '--username', secrets.stardog_user,
           '--passwd', secrets.stardog_password, '--format', 'JSON', db_name,
           query]
    print 'stardog_json_query()'
    print query
    out = subprocess.check_output(cmd)
    try:
        data = json.loads(out)
    except ValueError:
        sys.stdout.write(out)
        raise
    return data


def virtuoso_json_query(query):
    print 'virtuoso_json_query'
    print query
    return sparq.query(query, accept='application/sparql-results+json')


def sparql_json_query(query):
    if DB_TYPE == VIRT7_DB:
        return virtuoso_json_query(query)
    elif DB_TYPE == STARDOG_DB:
        return stardog_json_query(query)
    else:
        raise Exception('Unrecognized DB_TYPE', DB_TYPE)


def query_for_id_tuples(query, bindings):
    '''
    Make a sparql query, parse out the bindings (with no type conversion) from
    the results for the given bindings.  The binding values should be IRIs of
    the form "http://example.com/path/to/id".  Return a list of tuples of ids
    (not IRIs), where the ids are ordered according to `bindings`.
    '''
    result = sparql_json_query(query)
    # IRIs are like u'http://purl.targetscan.org/mir_family/miR-33'
    # or u'http://purl.targetscan.org/mir_family/miR-279%2F286%2F996'
    id_tuples = []
    for row in result['results']['bindings']:
        iris = [row[binding]['value'] for binding in bindings]
        ids = tuple(iri_to_id(iri) for iri in iris)
        id_tuples.append(ids)
    return id_tuples


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
    return [iri_to_id(iri, token) for iri in iris]

def iri_to_id(iri, token='/'):
    '''
    Convert a URIs/IRIs into a plain old identifier.  The id are also
    URL unquoted, which should be safe, considering that, as IRIs, they should
    be quoted, right?  Or is that just a URL thing?
    iri: str.  An IRI in the form: 'http://example.com/path/to/identifier'
    token: a character (or string) to split the identifier from the rest of
    the IRI.  In RDF, it is typically a slash ('/') or a hash ('#').
    '''
    return urllib.unquote(iri.rsplit(token, 1)[-1])


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
    '''
    The name of the database for the vanvactor mirna collaboration
    '''
    return 'mirna'


def drop_stardog_database():
    cmd = 'time stardog-admin db drop --username {} --passwd {} {}'.format(
        secrets.stardog_admin_user, secrets.stardog_admin_password,
        stardog_database_name())
    call(cmd)


def load_stardog_database():
    '''
    Data is loaded when the database is created b/c loading is much faster
    that way, according to the stardog docs.
    '''
    cmd = 'time stardog-admin db create --username {} --passwd {}'.format(
        secrets.stardog_admin_user, secrets.stardog_admin_password)
    cmd += ' --name {}'.format(stardog_database_name())
    # quickly bulk load data at database creation time
    # cmd += ' ' + ' '.join(all_rdf_paths())
    call(cmd)

    for filename, graph in all_rdf_files_and_graphs():
        cmd = 'time stardog data add --username {} --passwd {}'.format(
            secrets.stardog_admin_user, secrets.stardog_admin_password)
        cmd += ' --named-graph {} {} {}'.format(
            graph, stardog_database_name(), filename)
        call(cmd)


def stardog_default_graph_rdf_files_and_graphs():
    return [
        (affymetrix_fly_annotations_rdf_path(), affymetrix_ng),
        (flybase_rdf_path(), flybase_ng),
        (microcosm_rdf_path(), microcosm_ng),
        (mirbase_rdf_path(), mirbase_ng),
        (synaptomedb_v1_rdf_path(), synaptomedb_ng),
        (targetscan_rdf_path(), targetscan_ng),
        # UniProt from the idmapping.dat file
        # (uniprot_rdf_path(), uniprot_ng),
        # UniProt from beta.sparql.uniprot.org
        (uniprot_homo_ensembl_gene_rdf_path(), uniprot_ng),
        (uniprot_homo_ensembl_transcript_rdf_path(), uniprot_ng),
        (uniprot_homo_refseq_rdf_path(), uniprot_ng),
        (uniprot_flybase_transcript_rdf_path(), uniprot_ng),
        (uniprot_flybase_gene_rdf_path(), uniprot_ng),
        (uniprot_flybase_annotation_rdf_path(), uniprot_ng),
        (uniprot_homo_rdf_path(), uniprot_ng),
        (uniprot_fly_rdf_path(), uniprot_ng),
        # Roundup from the downloaded txt files
        (roundup_rdf_path(), roundup_ng),
        # Roundup from sparql.roundup.hms.harvard.edu
        (roundup_sparql_rdf_path(), roundup_ng),
    ]


def stardog_named_graph_rdf_files_and_graphs():
    return [
        (ibanez_five_prime_mir_homologs_rdf_path(), ibanez_five_prime_ng),
        (ibanez_fullseq_mir_homologs_rdf_path(), ibanez_full_sequence_ng),
        (mcneill_cns_screen_rdf_path(), mcneill_cns_ng),
        (mcneill_muscle_screen_rdf_path(), mcneill_muscle_ng),
    ]


def all_rdf_files_and_graphs():
    return [
        (affymetrix_fly_annotations_rdf_path(), affymetrix_ng),
        (flybase_rdf_path(), flybase_ng),
        (ibanez_five_prime_mir_homologs_rdf_path(), ibanez_five_prime_ng),
        (ibanez_fullseq_mir_homologs_rdf_path(), ibanez_full_sequence_ng),
        (mcneill_cns_screen_rdf_path(), mcneill_cns_ng),
        (mcneill_muscle_screen_rdf_path(), mcneill_muscle_ng),
        (microcosm_rdf_path(), microcosm_ng),
        (mirbase_rdf_path(), mirbase_ng),
        (synaptomedb_v1_rdf_path(), synaptomedb_ng),
        (targetscan_rdf_path(), targetscan_ng),
        # UniProt from the idmapping.dat file
        # (uniprot_rdf_path(), uniprot_ng),
        # UniProt from beta.sparql.uniprot.org
        (uniprot_homo_ensembl_gene_rdf_path(), uniprot_ng),
        (uniprot_homo_ensembl_transcript_rdf_path(), uniprot_ng),
        (uniprot_homo_refseq_rdf_path(), uniprot_ng),
        (uniprot_flybase_transcript_rdf_path(), uniprot_ng),
        (uniprot_flybase_gene_rdf_path(), uniprot_ng),
        (uniprot_flybase_annotation_rdf_path(), uniprot_ng),
        (uniprot_homo_rdf_path(), uniprot_ng),
        (uniprot_fly_rdf_path(), uniprot_ng),
        # Roundup from the downloaded txt files
        (roundup_rdf_path(), roundup_ng),
        # Roundup from sparql.roundup.hms.harvard.edu
        (roundup_sparql_rdf_path(), roundup_ng),
    ]


def all_rdf_paths():
    return [filename for filename, graph in all_rdf_files_and_graphs()]


def load_virtuoso_database():
    for filename, graph in all_rdf_files_and_graphs():
        # drop any existing graph
        virt7.drop_graph(graph, silent=True)
        # load files into graph, one at a time
        print 'loading graph {}'.format(graph)
        virt7.load_graph(graph, filename)


def load_rdf_database():
    if DB_TYPE == VIRT7_DB:
        load_virtuoso_database()
    elif DB_TYPE == STARDOG_DB:
        drop_stardog_database()
        load_stardog_database()
    else:
        raise Exception('Unrecognized DB_TYPE', DB_TYPE)


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
    '''
    Used to update a list of mirs from Elizabeth to their current mirbase ids,
    so they can be more easily used in the database queries
    '''
    lookup = update_mirbase_ids(original_validated_mirs)
    for mir in original_validated_mirs:
        print mir, lookup.get(mir)

    current_mirs = [lookup[mir] for mir in original_validated_mirs if mir in lookup]
    print current_mirs
    print len(current_mirs)
    return current_mirs


def human_to_fly_conserved_synapse_genes_table():
    '''
    Return a list of tuples of human ensembl gene and flybase gene.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?heg ?fbg
    WHERE {
    # human synaptic ensembl genes
    ?heg up:classifiedWith syndb:synaptic .
    ?heg up:organism taxon:9606 .
    ?heg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?heg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # flybase genes
    ?du rdfs:seeAlso ?fbg .
    ?fbg up:database db:flybase_gene .
    }
    '''
    return query_for_id_tuples(query, ['heg', 'fbg'])


def human_conserved_synapse_genes():
    '''
    Synapse genes are human genes from SynapseDB.
    Conserved synapse genes are the subset that have fly orthologs.
    Return a list of human ensembl synaptic genes with fly orthologs.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?heg
    WHERE {
    # human synaptic ensembl genes
    ?heg up:classifiedWith syndb:synaptic .
    ?heg up:organism taxon:9606 .
    ?heg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?heg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # flybase genes
    ?du rdfs:seeAlso ?dg .
    ?dg up:database db:flybase_gene .
    }
    '''
    return query_for_ids(query, 'heg')

def fly_conserved_synapse_genes():
    '''
    Synapse genes are human genes from SynapseDB.
    Conserved synapse genes are the subset that have fly orthologs.
    Return a list of fly genes that are orthologous to human synaptic genes.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?dg
    WHERE {
    # human synaptic ensembl genes
    ?heg up:classifiedWith syndb:synaptic .
    ?heg up:organism taxon:9606 .
    ?heg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?heg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # flybase genes
    ?du rdfs:seeAlso ?dg .
    ?dg up:database db:flybase_gene .
    }
    '''
    return query_for_ids(query, 'dg')

def targetscan_human_mirs_targeting_conserved_synapse_genes():
    '''
    Synapse genes are human genes from SynapseDB.
    Conserved synapse genes are the subset that have fly orthologs.
    Return the human mirs predicted by targetscan to target one or more of
    the conserved synapse genes.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?hm
    WHERE {
    # human synaptic genes
    ?hg up:classifiedWith syndb:synaptic .
    ?hg up:organism taxon:9606 .
    ?hu rdfs:seeAlso ?hg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

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
    Synapse genes are human genes from SynapseDB.
    Conserved synapse genes are the subset that have fly orthologs.
    Return the fly mirs that map to targetscan mir families predicted by
    targetscan to target one or more of fly orthologs of the conserved synapse
    genes.

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
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

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
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

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
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

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


def microcosm_predicted_human_mir_targets(mirbase_id):
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


def microcosm_predicted_fly_mir_targets(mirbase_id):
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
    Given a mirbase id, go to a mature mirbase accession, to a targetscan mir
    family, to predicted targeted refseq to uniprot to ensembl gene.
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


def human_homologs_of_fly_mir(mirbase_id, ibanez_method):
    mir = mirbase_id_iri(mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?hmid
    WHERE {{
    # original mirbase id to mature mirbase accession
    ?macc mb:has_id <{mid_orig}> .
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .

    # mirbase acc to all fly mirbase ids
    ?macc mb:has_id ?mid .
    ?mid up:database db:mirbase_id .
    ?mid up:organism taxon:7227 .

    GRAPH <{ibanez_method}> {{
        ?mid obo:homologous_to ?hmid_old .
    }}

    # human mirbase id
    ?hmid_old up:organism taxon:9606 .
    ?hmid_old up:database db:mirbase_id .

    # update mirbase id to current mirbase id
    ?hmacc mb:has_id ?hmid_old .
    ?hmacc up:database db:mirbase_acc .
    ?hmacc a mb:mature_mirna .
    ?hmacc mb:current_id ?hmid .
    ?hmid up:database db:mirbase_id .
    ?hmid up:organism taxon:9606 .
    }}
    '''.format(mid_orig=mir,
               ibanez_method=IBANEZ_METHOD_TO_NG[ibanez_method])
    return query_for_ids(query, 'hmid')


def targets_of_mirs(mirbase_ids, species, db):
    '''
    Return a list of the union of all the targets of the mirnas in mirbase_ids,
    for the given species and database.

    db: MICROCOSM or TARGETSCAN
    species: homo or mus
    mirbase_ids: a list.
    '''
    func_map = {(homo, MICROCOSM): microcosm_predicted_human_mir_targets,
                (homo, TARGETSCAN): targetscan_predicted_human_mir_targets,
                (dro, MICROCOSM): microcosm_predicted_fly_mir_targets,
                (dro, TARGETSCAN): targetscan_predicted_fly_mir_targets,}
    return map_flatten_set_list(func_map[(species, db)], mirbase_ids)


def fly_orthologs_of_targetscan_targets_of_human_mirs_homologous_to_fly_mir(
    mirbase_id, ibanez_method):
    '''
    '''
    mir = mirbase_id_iri(mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?fbg
    WHERE {{
    # original mirbase id to mature mirbase accession
    ?dmacc mb:has_id <{dmid_orig}> .
    ?dmacc a mb:mature_mirna .
    ?dmacc up:database db:mirbase_acc .

    # mirbase acc to all fly mirbase ids
    ?dmacc mb:has_id ?dmid .
    ?dmid up:database db:mirbase_id .
    ?dmid up:organism taxon:7227 .

    GRAPH <{ibanez_method}> {{
        ?dmid obo:homologous_to ?hmid .
    }}

    # human mirbase id
    ?hmid up:organism taxon:9606 .
    ?hmid up:database db:mirbase_id .

    # original mirbase id to mature mirbase accession
    ?hmacc mb:has_id ?hmid .
    ?hmacc a mb:mature_mirna .
    ?hmacc up:database db:mirbase_acc .

    # mirbase acc to targetscan mir family
    ?tmf rdfs:seeAlso ?hmacc .
    ?tmf up:database db:targetscan_mir_family .

    # targetscan mir family to targeted refseq transcripts
    ?tmf mb:targets ?ref .
    ?ref up:database db:ncbi_refseq .

    # refseq to uniprot
    ?hu rdfs:seeAlso ?ref .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # uniprot to flybase gene .
    ?du rdfs:seeAlso ?fbg .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    }}
    '''.format(dmid_orig=mir, ibanez_method=ibanez_method)
    return query_for_ids(query, 'gene')



def targets_of_human_homologs_of_fly_mir(mirbase_id, ibanez_method, target_db):
    '''
    Get the human homologs of a fly mir and then get the targets of those
    homologs.  Return the set of targets of the homologs, as a list.
    '''
    homolog_mirs = human_homologs_of_fly_mir(mirbase_id, ibanez_method)
    return targets_of_mirs(homolog_mirs, homo, target_db)

def targets_of_human_homologs_of_fly_mirs(mirbase_ids, ibanez_method,
                                          target_db):
    '''
    Combine the targets of each homolog of each fly mirbase id into one
    giant set, and retur the set as a list.
    '''
    func = functools.partial(targets_of_human_homologs_of_fly_mir,
                             ibanez_method=ibanez_method, target_db=target_db)
    return map_flatten_set_list(func, mirbase_ids)


def slow_fly_orthologs_of_human_genes(genes):
    return map_flatten_set_list(fly_orthologs_of_human_gene, genes)


def fly_orthologs_of_human_genes(genes):
    if not genes:
        return []

    gene_iris = [ensembl_gene_iri(gene) for gene in genes]
    query = prefixes() + '''
    SELECT DISTINCT ?fbg
    WHERE {{
    # human ensembl gene to uniprot
    ?heg up:organism taxon:9606 .
    ?heg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?heg .

    FILTER ?heg IN ({hegs})

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # uniprot to flybase gene .
    ?du rdfs:seeAlso ?fbg .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .

    }}
    '''.format(hegs=', '.join('<{}>'.format(iri) for iri in gene_iris))
    return query_for_ids(query, 'fbg')




def fly_orthologs_of_human_gene(gene):
    '''
    gene: a human ensembl gene id.
    returns: a list of flybase genes.
    '''
    gene_iri = ensembl_gene_iri(gene)
    query = prefixes() + '''
    SELECT DISTINCT ?fbg
    WHERE {{
    # human ensembl gene to uniprot
    <{heg}> up:organism taxon:9606 .
    <{heg}> up:database db:ensembl_gene .
    ?hu rdfs:seeAlso <{heg}> .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # uniprot to flybase gene .
    ?du rdfs:seeAlso ?fbg .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    }}
    '''.format(heg=gene_iri)
    return query_for_ids(query, 'fbg')


def affymetrix_to_flybase_table():
    query = prefixes() + '''
    SELECT DISTINCT ?affy ?fbg
    WHERE {
        ?affy up:database db:affymetrix_probeset .
        ?affy up:organism taxon:7227 .
        ?affy rdfs:seeAlso ?fbg .
        ?fbg up:database db:flybase_gene .
        ?fbg up:organism taxon:7227 .
    }
    '''
    return query_for_id_tuples(query, ['affy', 'fbg'])


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



def results_dir():
    dn = datetime.datetime.now().strftime('%Y%m%d') + '_results'
    return makedirs(os.path.join(config.datadir, dn))


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


def write_list_file(items, filename):
    makedirs(os.path.dirname(filename))
    with open(filename, 'w') as fh:
        fh.write(''.join([item + '\n' for item in items]))


def write_csv_file(rows, filename, headers=None):
    '''
    :param rows: list of rows, where each row is a list of value.
    :param headers: list of strings, used as a header row.  If None, no header
    row will be output.
    '''
    with open(filename, 'w') as fh:
        writer = csv.writer(fh)
        if headers:
            writer.writerow(headers)
        for row in rows:
            writer.writerow(row)


# PHASE I
# HUMAN and FLY, MICROCOSM and TARGETSCAN, mirs targeting conserved synapse genes.
# ALSO VALIDATED FLY MIRS

def write_targetscan_human_mirs_targeting_conserved_synapse_genes():
    dn = os.path.join(results_dir(), 'phase1')
    mirs = sorted(targetscan_human_mirs_targeting_conserved_synapse_genes())
    fn = os.path.join(dn, 'predicted_targetscan_human_mirs_targeting_conserved_synapse_genes.txt')
    write_list_file(mirs, fn)


def write_microcosm_human_mirs_targeting_conserved_synapse_genes():
    dn = os.path.join(results_dir(), 'phase1')
    mirs = sorted(microcosm_human_mirs_targeting_conserved_synapse_genes())
    fn = os.path.join(dn, 'predicted_microcosm_human_mirs_targeting_conserved_synapse_genes.txt')
    write_list_file(mirs, fn)


def write_overlap_between_microcosm_human_mirs_and_targetscan_human_mirs():
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'overlap_between_microcosm_and_targetscan_human_mirs.csv')
    micro = set(microcosm_human_mirs_targeting_conserved_synapse_genes())
    targ = set(targetscan_human_mirs_targeting_conserved_synapse_genes())
    write_set_overlap_file(micro, targ, 'microcosm', 'targetscan', fn)


def write_targetscan_fly_mirs_targeting_conserved_synapse_genes():
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'predicted_targetscan_fly_mirs_targeting_conserved_synapse_genes.txt')
    mirs = sorted(targetscan_fly_mirs_targeting_conserved_synapse_genes())
    write_list_file(mirs, fn)


def write_microcosm_fly_mirs_targeting_conserved_synapse_genes():
    '''
    Write a list of fly mirs that are predicted to target genes that are
    orthologous to human synaptic genes (according the synapsedb).
    '''
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'predicted_microcosm_fly_mirs_targeting_conserved_synapse_genes.txt')
    mirs = sorted(microcosm_fly_mirs_targeting_conserved_synapse_genes())
    write_list_file(mirs, fn)


def write_overlap_between_microcosm_fly_mirs_and_targetscan_fly_mirs():
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'overlap_between_microcosm_and_targetscan_fly_mirs.csv')
    micro = set(microcosm_fly_mirs_targeting_conserved_synapse_genes())
    targ = set(targetscan_fly_mirs_targeting_conserved_synapse_genes())
    write_set_overlap_file(micro, targ, 'microcosm', 'targetscan', fn)


def write_overlap_between_validated_fly_mirs_and_microcosm_predicted_fly_mirs():
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'overlap_between_validated_and_predicted_microcosm_fly_mirs.csv')
    predicted = set(microcosm_fly_mirs_targeting_conserved_synapse_genes())
    validated = set(validated_mirs)
    write_set_overlap_file(validated, predicted, 'validated', 'predicted', fn)


def write_overlap_between_validated_fly_mirs_and_targetscan_predicted_fly_mirs():
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'overlap_between_validated_and_predicted_targetscan_fly_mirs.csv')
    predicted = set(targetscan_fly_mirs_targeting_conserved_synapse_genes())
    validated = set(validated_mirs)
    write_set_overlap_file(validated, predicted, 'validated', 'predicted', fn)


def write_conserved_synaptic_genes():
    dn = os.path.join(results_dir(), 'conserved_synapse_genes')
    fn = os.path.join(dn, 'fly_conserved_synapse_genes.txt')
    genes = fly_conserved_synapse_genes()
    write_list_file(genes, fn)
    fn = os.path.join(dn, 'human_conserved_synapse_genes.txt')
    genes = human_conserved_synapse_genes()
    write_list_file(genes, fn)


def write_affymetrix_to_flybase_table():
    fn = os.path.join(results_dir(), 'affymetrix_to_flybase_genes.csv')
    write_csv_file(affymetrix_to_flybase_table(), fn,
                   headers=('affymetrix_probeset_id', 'flybase_gene'))

def write_human_to_fly_conserved_synapse_genes_table():
    dn = os.path.join(results_dir(), 'conserved_synapse_genes')
    fn = os.path.join(dn, 'human_to_fly_conserved_synapse_genes.csv')
    write_csv_file(human_to_fly_conserved_synapse_genes_table(), fn,
                   headers=('ensembl_gene', 'flybase_gene'))


def write_overlap_of_validated_mir_targets_and_conserved_targets_of_human_mir_homologs():
    dn = os.path.join(results_dir(), 'overlap_of_validated_mir_targets_and_conserved_targets_of_human_mir_homologs')
    # filename template
    fnt = 'overlap_of_{}_{}_targets_and_conserved_human_mir_{}_homolog_targets'
    for target_db in [MICROCOSM, TARGETSCAN]:
        for ibanez_method in [ibanez.FIVE_PRIME, ibanez.FULL_SEQ]:
            for mir in validated_mirs[:1]:
                print target_db, ibanez_method, mir
                fn = os.path.join(dn, fnt.format(mir, target_db, ibanez_method))
                fly_targets = set(targets_of_mirs([mir], dro, target_db))
                fly_ortholog_targets = set(fly_orthologs_of_targetscan_targets_of_human_mirs_homologous_to_fly_mir(
                    mir, ibanez_method))
                write_set_overlap_file(fly_targets, fly_ortholog_targets,
                                       'fly_targets', 'human_target_orthologs',
                                       fn)
                # homo_targets = targets_of_human_homologs_of_fly_mir(
                    # mir, ibanez_method, target_db)
                # fly_ortholog_targets = fly_orthologs_of_human_genes(homo_targets)


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
#
# FOR EACH VALIDATED FLY MIR:
#     PREDICTED TARGETS OF THE FLY MIR COMPARED TO (UNION OF PREDICTED TARGETS OF THE HUMAN HOMOLOGS OF THE FLY MIR)

def write_overlap_between_screened_and_microcosm_and_targetscan_fly_mir_targets():
    dn = os.path.join(results_dir(), 'overlap_of_screened_and_predicted_target_genes')
    for mir in screened_mirs:
        microcosm_targets = set(microcosm_predicted_fly_mir_targets(mir))
        targetscan_targets = set(targetscan_predicted_fly_mir_targets(mir))
        mt_fn = os.path.join(dn, '{}_overlap_between_microcosm_and_targetscan_targets.csv'.format(mir))
        write_set_overlap_file(microcosm_targets, targetscan_targets, 'microcosm', 'targetscan', mt_fn)
        for tissue in TISSUES:
            screened_targets = set(screened_fly_mir_targets_in_tissue(mir, tissue))
            sm_fn = os.path.join(dn, '{}_{}_overlap_between_screened_and_microcosm_targets.csv'.format(mir, tissue))
            write_set_overlap_file(screened_targets, microcosm_targets, 'screened', 'microcosm', sm_fn)
            st_fn = os.path.join(dn, '{}_{}_overlap_between_screened_and_targetscan_targets.csv'.format(mir, tissue))
            write_set_overlap_file(screened_targets, targetscan_targets, 'screened', 'targetscan', st_fn)


def write_overlap_of_screened_fly_mir_targets_and_nmj_rnai_genes():
    dn = os.path.join(results_dir(), 'overlap_of_nmj_rnai_genes_and_screened_fly_targets')
    nmj_rnai = set(get_nmj_rnai_genes())
    for mir in screened_mirs:
        for tissue in TISSUES:
            screened_targets = set(screened_fly_mir_targets_in_tissue(mir, tissue))
            fn = os.path.join(dn, '{}_{}_overlap_between_nmj_rnai_genes_and_screened_fly_targets.csv'.format(mir, tissue))
            write_set_overlap_file(nmj_rnai, screened_targets, 'nmj_rnai_genes', 'screened', fn)


def write_overlap_of_validated_fly_mir_targets_and_nmj_rnai_genes():
    dn = os.path.join(results_dir(), 'overlap_of_nmj_rnai_genes_and_predicted_fly_targets')
    nmj_rnai = set(get_nmj_rnai_genes())
    for mir in validated_mirs:
        microcosm_targets = set(microcosm_predicted_fly_mir_targets(mir))
        fn = os.path.join(dn, '{}_overlap_of_nmj_rnai_genes_and_microcosm_fly_targets.csv'.format(mir))
        write_set_overlap_file(nmj_rnai, microcosm_targets, 'nmj_rnai_genes', 'microcosm', fn)
        targetscan_targets = set(targetscan_predicted_fly_mir_targets(mir))
        fn = os.path.join(dn, '{}_overlap_of_nmj_rnai_genes_and_targetscan_fly_targets.csv'.format(mir))
        write_set_overlap_file(nmj_rnai, targetscan_targets, 'nmj_rnai_genes', 'targetscan', fn)



###############
# MAIN WORKFLOW


def workflow():

    # Clear out a graph loaded into the wrong graph URI
    # Slow for big graphs.
    # virt.clear_graph('http://purl.roundup.hms.harvard.edu/graph/uniprot_taxonomy')

    if False: # Done tasks
        pass

        download_microcosm_targets()
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
        write_ibanez_five_prime_mir_homologs_rdf()
        write_ibanez_fullseq_mir_homologs_rdf()
        write_mcneill_cns_screen_rdf()
        write_mcneill_muscle_screen_rdf()
        write_microcosm_rdf()
        write_mirbase_rdf()
        write_synaptomedb_rdf()
        write_targetscan_rdf()

        # Construct Roundup Orthologs from sparql.roundup.hms.harvard.edu
        download_roundup_from_sparql()
        # Or get them from the roundup text files from a raw data download
        write_roundup_orthologs_rdf()

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

        load_rdf_database()

        # PHASE I
        write_targetscan_human_mirs_targeting_conserved_synapse_genes()
        write_microcosm_human_mirs_targeting_conserved_synapse_genes()
        write_overlap_between_microcosm_human_mirs_and_targetscan_human_mirs()
        write_targetscan_fly_mirs_targeting_conserved_synapse_genes()
        write_microcosm_fly_mirs_targeting_conserved_synapse_genes()
        write_overlap_between_validated_fly_mirs_and_microcosm_predicted_fly_mirs()
        write_overlap_between_validated_fly_mirs_and_targetscan_predicted_fly_mirs()
        write_overlap_between_microcosm_fly_mirs_and_targetscan_fly_mirs()

        # PHASE II
        write_overlap_between_screened_and_microcosm_and_targetscan_fly_mir_targets()
        write_overlap_of_screened_fly_mir_targets_and_nmj_rnai_genes()
        write_overlap_of_validated_fly_mir_targets_and_nmj_rnai_genes()

        write_conserved_synaptic_genes()
        write_affymetrix_to_flybase_table()
        write_human_to_fly_conserved_synapse_genes_table()

    else:

        write_overlap_of_validated_mir_targets_and_conserved_targets_of_human_mir_homologs()

        pass



def main():
    workflow()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception('')
        raise


