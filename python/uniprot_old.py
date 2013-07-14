
'''
This code has been deprecated by extracting id mappings directly from
http://beta.sparql.uniprot.org/.

Download idmapping.dat file, split it into several smaller files (so it does
not take as long to traverse), filter the ids that I want and write them to
an RDF file.

These are the ids I want, the data source they are used to map to, the organism
I want them for, whether or not they are in the uniprot data and what file/id
type they are found in.

- mrna refseq ids (targetscan) (human) YES, idmapping.RefSeq_NT.dat (but need to strip version off the end)
- flybase annotation ids (microcosm, flybase) (fly) YES/MAYBE, idmapping.UCSC.dat
- flybase gene ids (targetscan) (fly) YES, idmapping.FlyBase.dat 
- flybase transcript ids (flybase) (fly), YES, idmapping.EnsemblGenome_TRS.dat
- ensembl transcript ids (microcosm) (human), YES, idmapping.Ensembl_TRS.dat
- ensembl gene ids (synaptomedb) (human), YES, idmapping.Ensembl.dat

The UniProt world differs from this project in that it sees the source of
FBtr* ids as EnsemblMetazoa/EnsemblGenome_TRS whereas I see them as FlyBase
Transcript ids.  This differing on names / database sources for ids exists
for multiple id types.

'''

###################
# UNIPROT FUNCTIONS

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


def uniprot_rdf_path(version='2013_04'):
    return os.path.join(config.datadir, 'uniprot', version,
                        'uniprot-{}-idmapping.nt'.format(version))


def uniprot_idmapping_path(version='2013_04', id_type=None):
    d = os.path.join(config.datadir, 'uniprot', version)
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


