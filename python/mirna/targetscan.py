
import os

import rdflib

import config
from mirna.core import (db_pred, dro_iri, flybase_gene_db, flybase_gene_iri,
                        homo_iri, mirbase_acc_db, mirbase_acc_iri,
                        organism_pred, refseq_db, refseq_iri, see_also_pred,
                        targets_pred, targetscan_mir_family_db,
                        targetscan_mir_family_iri)
from mirna.util import (download_and_unzip)



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

    # converting fly symbol to flybase gene id
    fly_symbol_to_gene = targetscan_fly_symbol_to_gene()

    # link family to fly gene target
    for family, symbol in gen_targetscan_fly_predicted_targets():
        fam = targetscan_mir_family_iri(family)
        gene = flybase_gene_iri(fly_symbol_to_gene[symbol])
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


