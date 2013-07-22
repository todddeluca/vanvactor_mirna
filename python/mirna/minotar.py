
import os

import config
from mirna.util import download_and_unzip
from mirna.targetscan import targetscan_top_dir


#########
# MinoTar


def minotar_dir():
    # version is made up b/c I can not find a version on the minotar site.
    return os.path.join(config.datadir, 'minotar', 'v1')


def minotar_fly_dir():
    return os.path.join(minotar_dir(), 'drosophila_melanogaster')


def minotar_human_dir():
    return os.path.join(minotar_dir(), 'homo_sapiens')


def minotar_fly_csv_path():
    return os.path.join(minotar_fly_dir(), 'DrosConservedSeedsConservedTargs.csv')


def minotar_human_csv_path():
    return os.path.join(minotar_human_dir(), 'HumanConservedSeedsConservedTargs.csv')


def download_minotar():
    '''
    For fly:
    The columns in the target prediction files are:
        Gene ID
        Transcript ID
        Symbol
        CG identifier
        Conserved 8Mers
        Conserved 7Mers-m8
        Conserved 7Mers-1A
        Non-Conserved 8Mers
        Non-Conesrved 7Mers-m8
        Non-Conserved 7Mers-1A
        Probability Conserved Targeting

    For human:
    The columns in the target prediction files are:
        RefSeq ID
        Gene Symbol
        Conserved 8Mers
        Non-Conserved 8Mers
        Non-Conserved 7Mers-m8
        NonConserved 7Mers-1A
        Probability Conserved Targeting
    '''
    # Minotar (v1) was generated using mir families from fly_52orfs and vert_50
    # according to personal correspondence with Bonnie Berger.
    # Download the mir family info for those datasets in addition to the
    # minotar data.


    print 'Downloadding MinoTar supporting data from TargetScanVert version vert_50.'
    vert50_dir = os.path.join(targetscan_top_dir(), '5.2', 'vert_50')
    url = 'http://www.targetscan.org//vert_50/vert_50_data_download/miR_Family_Info.txt.zip'
    download_and_unzip(url, vert50_dir)

    print 'Downloading MinoTar supporting data from TargetScanFly version fly_52orfs.'
    orfs_dir = os.path.join(targetscan_top_dir(), '6.2', 'fly_52orfs')
    url = 'http://www.targetscan.org/fly_52orfs/fly_52orfs_data_download/miR_Family_Info.txt.zip'
    download_and_unzip(url, orfs_dir)
    url = 'http://www.targetscan.org/fly_52orfs/fly_52orfs_data_download/ORF_Sequences.txt.zip'
    download_and_unzip(url, orfs_dir)
    url = 'http://www.targetscan.org/fly_52orfs/fly_52orfs_data_download/Predicted_Targets_Info.txt.zip'
    download_and_unzip(url, orfs_dir)

    print 'Downloading MinoTar data files.'
    url = 'http://www.flyrnai.org/supplement/HumanUniqueMirnas.txt'
    download_and_unzip(url, minotar_human_dir())
    url = 'http://www.flyrnai.org/supplement/HumanConservedSeedsConservedTargs.zip'
    download_and_unzip(url, minotar_human_dir())
    url = 'http://www.flyrnai.org/supplement/DrosUniqueMirnas.txt'
    download_and_unzip(url, minotar_fly_dir())
    url = 'http://www.flyrnai.org/supplement/DrosConservedSeedsConservedTargs.zip'
    download_and_unzip(url, minotar_fly_dir())


def concatenate_minotar_fly():
    '''
    Concatenate the dozens of files downloaded from minotar, each of which
    encodes the targetscan mir family in the filename, into a single csv file
    where the (unencoded) mir family is prepended as the first field of each
    line.
    '''
    targetfiles_dir = os.path.join(minotar_fly_dir(), 'ConservedSeedsConservedTargs')
    def family_parser(fn):
        return os.path.splitext(os.path.basename(fn).replace('_', '/'))[0]

    concatenate_minotar_sub(targetfiles_dir, family_parser, minotar_fly_csv_path())


def concatenate_minotar_human():
    '''
    Concatenate the dozens of files downloaded from minotar, each of which
    encodes the targetscan mir family in the filename, into a single csv file
    where the (unencoded) mir family is prepended as the first field of each
    line.
    '''
    targetfiles_dir = os.path.join(minotar_human_dir(), 'ConservedSeedsConservedTargs')
    def family_parser(fn):
        return os.path.splitext(os.path.basename(fn).replace(':', '/'))[0]

    concatenate_minotar_sub(targetfiles_dir, family_parser, minotar_human_csv_path())


def concatenate_minotar_sub(files_dir, family_parser, outfile):
    with open(outfile, 'w') as fh:
        for fn in os.listdir(files_dir):
            if not fn.endswith('.txt'):
                print 'Warning: file not ending in ".txt" dir={}, file={}'.format(files_dir, fn)
                continue
            mir_family = family_parser(fn)
            print mir_family
            with open(os.path.join(files_dir, fn)) as infh:
                for line in infh:
                    if not line.strip():
                        raise Exception('Unexpected blank line', fn, line)
                    if line.strip().startswith('#'):
                        raise Exception('Unexpected comment line', fn, line)
                    fh.write('{},{}'.format(mir_family, line))






