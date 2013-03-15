import codecs
import datetime
import os

import requests


def download(destdir, kind, version, mir, genome, datestamp=None):
    '''
    :param kind: The format of the file to download, 'txt' or 'gff'.
    :param version: The version of the microcosm database.  e.g. 'v5'
    :param mir: e.g. dme-miR-13b
    :param genome: A microcosm genome id.  Why didn't they use ncbi taxon ids
    like the rest of the world?  e.g. 2188 for drosophila melanogaster.  
    :param datestamp: added to the filename for additional versioning.  e.g. "20130214"
    '''
    if datestamp is None:
        datestamp = datetime.datetime.today().strftime("%Y%m%d")
    url = "http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/{version}/download_formatter.pl?format={kind}&genome_id={genome}&mirna_id={mir}"
    path = os.path.join(destdir, "{version}/microcosm-{version}-{mir}-targets-{datestamp}.{kind}")
    url = url.format(kind=kind, version=version, mir=mir, genome=genome)
    path = path.format(kind=kind, version=version, mir=mir, genome=genome,
                       datestamp=datestamp)
    print "saving {} to {}".format(url, path)
    r = requests.get(url)
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
        pass
    with codecs.open(path, "w", "utf-8") as fh:
        fh.write(r.text)
        pass
    return path


def main():
    import config

    for format in ['txt', 'gff']:
        for mir in ["dme-miR-2b", "dme-miR-13b"]:
            download(destdir=config.microcosm_dir, kind=format, version='v5',
                     mir=mir, genome='2188')


if __name__ == '__main__':
    main()



