
import os

import config


############################
# NMJ_RNAI LOF/GOF GENE LIST

def nmj_rnai_set_path():
    return os.path.join(config.datadir, 'NMJ RNAi Search File.txt')


def nmj_rnai_gain_of_function_set_path():
    return os.path.join(config.datadir, 'NMJ_RNAi_gain_of_function_flybase_ids.txt')


def get_nmj_rnai_genes():
    '''
    Return a list of flybase gene ids.

    The hits from three different screens (Aaron D'Antonio, Sanyal and
    Featherstone).  This contains FBgn IDs, which can be converted to gene
    symbols using flybase ID converter 
    '''
    path = nmj_rnai_set_path()
    print path
    with open(path) as fh:
        # Skip first line, the header
        genes = [line.strip() for i, line in enumerate(fh) if i > 0 and line.strip()]
    return genes


def get_nmj_rnai_gain_of_function_genes():
    '''
    Return a list of flybase gene ids.

    The gain of function genes should be a curated subset of the NMJ RNAi
    genes.  They were defined in a file Elizabeth McNeill sent, 
    "NMJ RNAi Gainoffunctionscreens.xlsx".

    The hits from three different screens (Aaron D'Antonio, Sanyal and
    Featherstone).  This contains FBgn IDs, which can be converted to gene
    symbols using flybase ID converter.
    '''
    path = nmj_rnai_gain_of_function_set_path()
    print path
    with open(path) as fh:
        # Skip first line, the header
        genes = [line.strip() for i, line in enumerate(fh) if i > 0 and line.strip()]
    return genes




