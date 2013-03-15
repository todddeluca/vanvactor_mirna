


import argparse
import logging
import sys



def select_flybase_gene_ids(gene_conversion_table):
    '''
    Return a list of unique flybase gene ids from the gene conversion table
    downloaded from flybase, skipping ids that did not convert.
    '''
    uniques = set()
    for i, line in enumerate(open(gene_conversion_table)):
        # skip values that flybase failed to convert
        if "unknown ID" in line:
            continue

        # Example line
        # CG10005-RA      FBtr0082507             FBgn0037972     CG10005
        # Fields: Submitted ID, Current ID, Converted ID, Related record
        splits = line.strip().split("\t")
        submitted, current, mystery_field, converted, related = splits

        # assuming this is a gene conversion table, then flybase converted the
        # submitted id into a flybase gene id.
        gene = converted
        assert gene.startswith("FBgn")
        uniques.add(gene)

    genes = sorted(uniques)
    return genes



