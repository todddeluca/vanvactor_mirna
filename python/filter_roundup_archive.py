
"""
For the paper, a working link to the roundup data is (reasonably) required. However the existing link to the
Roundup website is broken, since the website was decommissioned around 2014, and it is now 2020. The good news
is that the entire roundup dataset version 4 is available on o2.hms.harvard.edu. Unfortunately the orthologs file
for divergence=0.8 and p-value=1e-5 is 13GB. We only need orthologs for homo sapiens (9606) and drosophila
melanogaster (7227) and data repositories like Zenodo have file size limits. (And storing an extra 13GB is wasteful.)

The solution is to extract only the orthologs for hse and dme into a file and archive that in zenodo.

Orthologs for Roundup v4 computed using a divergence theshold of 0.8 and evalue threshold of 1e-5 were downloaded
from O2, the new version of the orchestra compute cluster at Harvard medical school. Specifically:

    scp o2:/n/groups/public+cbi/sites/roundup/datasets/4/download/roundup-4-orthologs_0.8_1e-5.txt.gz ~/data/roundup/datasets/4/.

The file was unzipped and this script was run to create the file roundup-v4-orthologs_for_7227_9606_0.8_1e-5.txt,
which was then gzipped and uploaded to Zenodo.
"""

from pathlib import Path


def main():
    div = '0.8'
    evalue = '1e-5'
    dme_taxon = '7227'
    hse_taxon = '9606'
    # dme_taxon = '329726'
    # hse_taxon = '99883'
    # dme_taxon = '243243'
    # hse_taxon = '583346'
    full_archive_path = Path(f'~/data/roundup/datasets/4/roundup-4-orthologs_{div}_{evalue}.txt').expanduser()
    filtered_archive_path = Path(f'~/data/roundup/datasets/4/roundup-v4-orthologs_for_{dme_taxon}_{hse_taxon}_{div}_{evalue}.txt').expanduser()

    state = 'SEARCH'
    with open(full_archive_path) as infh, open(filtered_archive_path, 'w') as outfh:
        for line in infh:
            if line.startswith('PA'):
                pa, qt, st, d, e = line.strip().split()
                if qt == dme_taxon and st == hse_taxon and d == div and e == evalue:
                    outfh.write(line)
                    state = 'FOUND'
            elif state == 'FOUND':
                outfh.write(line)
                if line.startswith('//'):
                    break  # done
            else:
                continue


if __name__ == '__main__':
    main()
