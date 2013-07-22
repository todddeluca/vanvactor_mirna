'''
Shared RDF URIs, utility functions, etc.
'''

import urllib

import rdflib


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


