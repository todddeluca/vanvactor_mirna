
import os
import re

import requests
import rdflib

import config
from mirna.core import (db_pred, ensembl_gene_db, ensembl_gene_iri,
                        ensembl_transcript_db, ensembl_transcript_iri,
                        flybase_annotation_db, flybase_annotation_iri,
                        flybase_gene_db, flybase_gene_iri,
                        flybase_transcript_db, flybase_transcript_iri,
                        refseq_db, refseq_iri, see_also_pred, uniprot_db,
                        uniprot_iri, )
from mirna.flybase import flybase_annotation_id_regex
from mirna.util import (makedirs, download_sparql_ntriples,
                        count_and_write_ntriples)

uniprot_sparql_url = 'http://beta.sparql.uniprot.org'



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


def download_uniprot_sparql_rdf(construct_phrase, where_phrase, taxon,
                                filename, min_count, substitution = None):
    '''
    construct_phrase: str.  A query fragment for what triples to create.
    taxon: e.g '7227' (fly) or '9606' (human)
    where_phrase: str. A query fragment for what pattern to match.  It is
    appended to patterns to look for up:Protein triples with up:organism
    taxon:<taxon> triples.
    substitution: if not None, a tuple of pattern and replacement to be used by
    an re.sub() call to alter each triple line.  E.g.:
        (r'http://purl.uniprot.org/ucsc/', r'http://purl.flybase.org/annotation_id/')
    min_count: int.  Minor sanity check and bias.  If fewer than min_count
    triples are returned, an exception will be raised.
    '''
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
        CONSTRUCT {{
        {construct}
        }}
        WHERE {{
          ?protein a up:Protein .
          ?protein up:organism taxon:{taxon} .
          {where}
        }}
    '''.format(construct=construct_phrase, taxon=taxon, where=where_phrase)
    nt = download_sparql_ntriples(qry, uniprot_sparql_url)
    if substitution:
        # Change UniProt URIs to my URIs
        def fix(line):
            return re.sub(substitution[0], substitution[1], line)
        nt = ''.join(fix(line) for line in nt.splitlines(True))
    count = count_and_write_ntriples(nt, filename)
    # As of 2013/07 there were 14953 fly flybase annotations.
    if count < min_count:
        raise Exception('Expected >= {} triples.'.format(min_count), count)


def download_uniprot_fly_rdf():
    fn = uniprot_fly_rdf_path()
    taxon = '7227'
    construct = '''
        ?protein a up:Protein .
        ?protein up:organism taxon:{taxon} .
        ?protein up:database <{db}> .
        '''.format(taxon=taxon, db=uniprot_db)
    where = ''
    # As of 2013/07 there were 37828 fly uniprots.
    min_count = 30000 * 3
    download_uniprot_sparql_rdf(construct, where, taxon, fn, min_count)


def download_uniprot_homo_rdf():
    fn = uniprot_homo_rdf_path()
    taxon = '9606'
    construct = '''
        ?protein a up:Protein .
        ?protein up:organism taxon:{taxon} .
        ?protein up:database <{db}> .
        '''.format(taxon=taxon, db=uniprot_db)
    where = ''
    # As of 2013/07 there were 134288 homo uniprots.
    min_count = 130000 * 3
    download_uniprot_sparql_rdf(construct, where, taxon, fn, min_count)


def download_uniprot_flybase_annotation_rdf():
    fn = uniprot_flybase_annotation_rdf_path()
    taxon = '7227'
    construct = '''
        ?protein rdfs:seeAlso ?xref .
        ?xref up:organism taxon:{taxon} .
        ?xref up:database <{db}> .
        '''.format(taxon=taxon, db=flybase_annotation_db)
    where = ''' 
        ?protein rdfs:seeAlso ?xref .
        ?xref up:database db:UCSC .
        '''
    # As of 2013/07 there were 134288 homo uniprots.
    min_count = 14000 * 3
    sub = (r'http://purl.uniprot.org/ucsc/',
           r'http://purl.flybase.org/annotation_id/')
    download_uniprot_sparql_rdf(construct, where, taxon, fn, min_count, sub)


def download_uniprot_flybase_gene_rdf():
    fn = uniprot_flybase_gene_rdf_path()
    taxon = '7227'
    construct = '''
        ?protein rdfs:seeAlso ?xref .
        ?xref up:organism taxon:{taxon} .
        ?xref up:database <{db}> .
        '''.format(taxon=taxon, db=flybase_gene_db)
    where = ''' 
        ?protein rdfs:seeAlso ?xref .
        ?xref up:database db:FlyBase .
        '''
    # As of 2013/07 there were 23348 fly protein + flybase gene annotations.
    # And 13k unique flybase genes in those annotation
    # Wrote 50988 triples.
    min_count = 50000
    sub = (r'http://purl.uniprot.org/flybase/',
           r'http://purl.flybase.org/gene_id/')
    download_uniprot_sparql_rdf(construct, where, taxon, fn, min_count, sub)


def download_uniprot_flybase_transcript_rdf():
    # UniProt annotates proteins with flybase transcripts (FBtr* ids) using
    # the database EnsemblMetazoa.
    fn = uniprot_flybase_transcript_rdf_path()
    taxon = '7227'
    construct = '''
        ?protein rdfs:seeAlso ?xref .
        ?xref up:organism taxon:{taxon} .
        ?xref up:database <{db}> .
        '''.format(taxon=taxon, db=flybase_transcript_db)
    where = ''' 
        ?protein rdfs:seeAlso ?xref .
        ?xref up:database db:EnsemblMetazoa .
        '''
    # As of 2013/07 there were 21819 fly flybase transcript annotations.
    # Wrote 65457 triples.
    min_count = 60000
    sub = (r'http://purl.uniprot.org/ensemblmetazoa/',
           r'http://purl.flybase.org/transcript_id/')
    download_uniprot_sparql_rdf(construct, where, taxon, fn, min_count, sub)


def download_uniprot_homo_refseq_rdf():
    # Construct human protein refseq annotation triples.
    # The NP_* ids are URIs.  The NM_* ids are comment strings on the
    # NP_* URIs.  Weird.  Also, strip the version number off the refseq ids,
    # since the TargetScan transcript ids do not have the version number
    # suffix.
    fn = uniprot_homo_refseq_rdf_path()
    taxon = '9606'
    construct = '''
        ?protein rdfs:seeAlso ?xrefuri .
        ?protein rdfs:seeAlso ?xref2uri .
        ?xrefuri up:organism taxon:{taxon} .
        ?xref2uri up:organism taxon:{taxon} .
        ?xrefuri up:database <{db}> .
        ?xref2uri up:database <{db}> .
        '''.format(taxon=taxon, db=refseq_db)
    where = ''' 
        ?protein rdfs:seeAlso ?xref . 
        ?xref up:database db:RefSeq .
        ?xref rdfs:comment ?xref2 .
        BIND (URI(REPLACE(str(?xref), 
          '(\\\\d+)\\\\.\\\\d+$', '$1')) as ?xrefuri)
        BIND (URI(CONCAT('http://purl.uniprot.org/refseq/', 
          REPLACE(?xref2, '(\\\\d+)\\\\.\\\\d+$', '$1'))) as ?xref2uri)
        '''
    # As of 2013/07 there were 21819 fly flybase transcript annotations.
    # Wrote 65457 triples.
    min_count = 80000 * 3
    sub = (r'http://purl.uniprot.org/refseq/',
           r'http://purl.ncbi.nlm.nih.gov/refseq/')
    download_uniprot_sparql_rdf(construct, where, taxon, fn, min_count, sub)


def download_uniprot_homo_ensembl_transcript_rdf():
    fn = uniprot_homo_ensembl_transcript_rdf_path()
    taxon = '9606'
    construct = '''
        ?protein rdfs:seeAlso ?xref .
        ?xref up:organism taxon:{taxon} .
        ?xref up:database <{db}> .
        '''.format(taxon=taxon, db=ensembl_transcript_db)
    where = ''' 
        ?protein rdfs:seeAlso ?xref .
        ?xref up:database db:Ensembl .
        ?xref a up:Transcript_Resource .
        '''
    # As of 2013/07 there were 98088 homo ensembl transcript annotations.
    # Wrote 294078 triples.
    min_count = 90000 * 3
    sub = (r'http://purl.uniprot.org/ensembl/',
           r'http://purl.ensembl.org/transcript/')
    download_uniprot_sparql_rdf(construct, where, taxon, fn, min_count, sub)


def download_uniprot_homo_ensembl_gene_rdf():
    fn = uniprot_homo_ensembl_gene_rdf_path()
    taxon = '9606'
    construct = '''
        ?protein rdfs:seeAlso ?xref .
        ?xref up:organism taxon:{taxon} .
        ?xref up:database <{db}> .
        '''.format(taxon=taxon, db=ensembl_gene_db)
    where = ''' 
        ?protein rdfs:seeAlso ?enst .
        ?enst up:database db:Ensembl .
        ?enst a up:Transcript_Resource .
        ?enst up:transcribedFrom ?xref .
        '''
    # As of 2013/07 there were 72427 homo ensembl gene annotations.
    # Wrote 216873 triples.  Note that this is odd, b/c there are
    # 71237 distinct pairs of ?protein and ?xref, and 22092 distinct
    # ?xrefs in the query above.  So we would expect 71237 + 2 * 22092 triples
    # (= 115421).  There are that many distinct triples, so the construct
    # query returned non-distinct triples.  Weird.
    min_count = 70000 * 3
    sub = (r'http://purl.uniprot.org/ensembl/',
           r'http://purl.ensembl.org/gene/')
    download_uniprot_sparql_rdf(construct, where, taxon, fn, min_count, sub)


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
    Generate triples mapping uniprot ids to other ids that we are interested 
    in, like RefSeq mRNA ids, FlyBase gene ids, FlyBase transcript ids, 
    Human Ensembl transcript ids, Human Ensembl gene ids, and FlyBase
    annotation ids.  For each mapping, generate:

    - a "seeAlso" triple between the uniprot id and the mapped id
    - a "database" triple for the mapped id

    Also, for each unique uniprot id, generate:

    - a "database" triple for the uniprot id

    Note that "organism" triples are not generated, because with the exception
    of the Ensembl ids, it is not possible to know specifically the taxon the
    ids come from.
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
        (refseq_iri, 'RefSeq_NT', re.compile(r'^(NM_\d+)\.\d+$'), refseq_db),
        (flybase_annotation_iri, 'UCSC', flybase_annotation_id_regex, flybase_annotation_db),
        (flybase_gene_iri, 'FlyBase', re.compile(r'(FBgn\d+)$'), flybase_gene_db),
        (flybase_transcript_iri, 'EnsemblGenome_TRS', re.compile(r'(FBtr\d+)$'), flybase_transcript_db),
        (ensembl_transcript_iri, 'Ensembl_TRS', re.compile(r'(ENST\d+)$'), ensembl_transcript_db),
        (ensembl_gene_iri, 'Ensembl', re.compile(r'(ENSG\d+)$'), ensembl_gene_db),
    ]

    graph = rdflib.Graph()
    uniprots = set()

    # id_iri makes an IRI/URI for the id.
    # id_type specifies what file has the uniprot to id mapping.
    # id_regex matches the specific ids in the file that we want, since some
    #   files have data for many kinds of ids or species.  Also used to select
    #   only part of an id.
    # id_database is the IRI for the database of origin of the ids (e.g
    #   flybase genes use 'http://purl.example.com/database/flybase_gene').
    for id_iri, id_type, id_regex, id_database in xref_data:
        print 'Doing {} in {}'.format(id_type, id_database)
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
            uniprots.add(uni)

    for uni in uniprots:
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


