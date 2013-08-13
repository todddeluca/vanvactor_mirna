

import os

import rdflib

import config
from mirna.util import makedirs, download_sparql_construct_rdf
from mirna.core import (dro_iri, dro_taxon, homo_iri, homo_taxon,
                        organism_pred, orthologous_to_pred, uniprot_iri)


###################
# ROUNDUP FUNCTIONS

def download_roundup_from_sparql():
    qry = '''
    PREFIX up:<http://purl.uniprot.org/core/> 
    PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
    PREFIX obo:<http://purl.org/obo/owl/obo#>
    CONSTRUCT { 
        ?gene1 obo:orthologous_to ?gene2 .
        ?gene2 obo:orthologous_to ?gene1 .
        ?gene1 up:organism taxon:7227 .
        ?gene2 up:organism taxon:9606 .
    }
    WHERE {
        ?gene1 obo:orthologous_to ?gene2 .
        ?gene1 up:organism taxon:7227 .
        ?gene2 up:organism taxon:9606 .
    }
    '''
    url = 'http://sparql.roundup.hms.harvard.edu:8890/sparql'
    filename = roundup_sparql_rdf_path()
    makedirs(os.path.dirname(filename))
    return download_sparql_construct_rdf(qry, filename, endpoint=url)


def gen_roundup_orthologs(query_taxon, subject_taxon, divergence='0.8', evalue='1e-5',
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


def roundup_sparql_rdf_path(taxon1='7227', taxon2='9606', version='qfo_2013_04'):
    # lexicographically normalize taxons
    taxon1, taxon2 = sorted([taxon1, taxon2])
    d = os.path.join(config.datadir, 'roundup', version)
    return os.path.join(d, 'roundup-{}_{}-orthologs.nt'.format(taxon1, taxon2))


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
    print 'write_orthologs_rdf'
    # parameters
    version = '4'
    divergence = '0.8'
    evalue = '1e-5'

    done_ids = set()
    graph = rdflib.Graph()
    for qid, sid, distance in gen_roundup_orthologs(dro_taxon, homo_taxon,
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


