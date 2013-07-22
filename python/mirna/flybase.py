

import elementtree
import os
import re

import rdflib

import config
from mirna.core import (db_pred, dro_taxon, flybase_annotation_db,
                        flybase_annotation_iri, flybase_transcript_db,
                        flybase_transcript_iri, see_also_pred)


###################
# FLYBASE FUNCTIONS


# Notes: CG50-RF is an real annotation id that violates this pattern.
flybase_annotation_id_regex = re.compile(r'^(C(G|R)\d+-R[A-Z]+)$')


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


