

import os

from deployenv import datadir, virtuoso_load_dir
# a temporary directory used to load data into virtuoso
# files are copied here, loaded, then deleted.
virtuoso_load_dir

microcosm_dir = os.path.join(datadir, 'microcosm')


