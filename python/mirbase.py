
import os
import re

import rdflib

import downloadutil # contains download, gunzip, call


# Real Entities
dro_iri = rdflib.URIRef('http://purl.uniprot.org/taxonomy/7227')
homo_iri = rdflib.URIRef('http://purl.uniprot.org/taxonomy/9606')

# Real Predicates
organism_pred = rdflib.URIRef('http://purl.uniprot.org/core/organism')
type_pred = rdflib.URIRef('http://www.w3.org/1999/02/22-rdf-syntax-ns#type')
db_pred = rdflib.URIRef('http://purl.uniprot.org/core/database')

# Fake Classes
mature_mirna_cl = rdflib.URIRef('http://purl.mirbase.org/owl/mature_mirna')
# A real class for mature microRNAs exists in Sequence Ontology.  I am not
# using it because it is not readable (in a SPARQL query).  I did ask on the
# SO mailing list and OBO discuss list and Andy from flybase says that yes
# they use it to annotate mirbase mature mirnas.
# http://purl.obolibrary.org/obo/SO_0000276

# Fake Predicates
is_mature_pred = rdflib.URIRef('http://purl.mirbase.org/owl/is_mature')
current_id_pred = rdflib.URIRef('http://purl.mirbase.org/owl/current_id')
has_id_pred = rdflib.URIRef('http://purl.mirbase.org/owl/has_id')
targets_pred = rdflib.URIRef('http://purl.mirbase.org/owl/targets')

# Fake Database IRIs.  Databases represent the origin of a sequence/entity, the
# scope of an id.  There is probably a better way to represent this, since this
# seems to conflate the organization and id/sequence type.
mirbase_acc_db = rdflib.URIRef('http://purl.example.com/database/mirbase_acc')
mirbase_id_db = rdflib.URIRef('http://purl.example.com/database/mirbase_id')

def mirbase_acc_iri(acc):
    return rdflib.URIRef('http://purl.mirbase.org/mirna_acc/{}'.format(acc))

def mirbase_id_iri(mirna_id):
    return rdflib.URIRef('http://purl.mirbase.org/mirna_id/{}'.format(mirna_id))


# The URL containing mirbase predicates and classes
owl_prefix = 'http://purl.mirbase.org/owl/'


class Mirbase:

    def __init__(self, version, datadir, name='mirbase'):
        self.version = version
        self.name = name
        self.datadir = datadir

    def mirbase_aliases_url(self):
        return 'ftp://mirbase.org/pub/mirbase/{}/aliases.txt.gz'.format(
            self.version)


    def mirbase_aliases_dest(self):
        return os.path.join(self.datadir, 'mirbase', self.version,
                            os.path.basename(self.mirbase_aliases_url()))


    def mirbase_aliases_filename(self):
        return self.mirbase_aliases_dest().rstrip('.gz')


    def mirbase_rdf_path(self):
        return os.path.join(self.datadir, 'mirbase', self.version,
                            'mirbase.nt')


    def download_mirbase_aliases(self):

        url = self.mirbase_aliases_url()
        dest = self.mirbase_aliases_dest()
        downloadutil.download(url, dest)
        downloadutil.gunzip(dest)

    def mirbase_rdf_mediatype(self):
        # N-Triples: 'text/plain'
        return 'text/plain'

    def write_mirbase_rdf(self):
        graph = rdflib.Graph()
        human_accs = set()
        human_ids = set()
        fly_accs = set()
        fly_ids = set()

        mature_re = re.compile(r'MIMAT\d+$')
        # link mirbase accessions to mirbase ids and collect the accession and ids
        with open(self.mirbase_aliases_filename()) as fh:
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

        with open(self.mirbase_rdf_path(), 'w') as outfh:
            outfh.write(graph.serialize(format='nt'))

