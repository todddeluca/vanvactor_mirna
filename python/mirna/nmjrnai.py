
import os

import config


############################
# NMJ_RNAI LOF/GOF GENE LIST

def nmj_rnai_set_path():
    return os.path.join(config.datadir, 'NMJ RNAi Search File.txt')


def get_nmj_rnai_genes():
    '''
    Return a list of flybase gene.

    The hits from three different screens (Aaron D'Antonio, Sanyal and
    Featherstone).  This contains FBgn IDs, which can be converted to gene
    symbols using flybase ID converter 
    '''
    path = nmj_rnai_set_path()
    print path
    with open(nmj_rnai_set_path()) as fh:
        genes = [line.strip() for i, line in enumerate(fh) if i > 0 and line.strip()]
        print genes
    return genes



