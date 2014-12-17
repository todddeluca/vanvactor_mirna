
'''
Contains code for creating a semantic database and querying it to answer Van Vactor
miRNA collaboration questions.


The "short_" prefix indicates sparql queries based on the simplified edges that
map mirbase id to gene (mir target prediction), gene to gene (orthology), and
mirbase id to mirbase id (homology)

The "merged_" prefix indicates sparql queries that combine targets for
microcosm and targetscan and for ibanez five-prime and full-sequence homology
methods.
'''

import cStringIO
import collections
import csv
import datetime
import functools
import itertools
import json
import logging
import os
import re
import subprocess
import sys
import urllib

import rdflib
import requests.exceptions
import scipy.stats

import util
import config
import secrets
import virtuoso
# import stardog
import sparql
from mirna import mcneill, ibanez
from mirna.util import makedirs, map_flatten_set_list, call
from mirna.core import (mirbase_id_iri, ensembl_gene_iri, dro_taxon,
                        homo_taxon, mus_taxon, cel_taxon)
from mirna.uniprot import (uniprot_fly_rdf_path,
                           uniprot_flybase_annotation_rdf_path,
                           uniprot_flybase_gene_rdf_path,
                           uniprot_flybase_transcript_rdf_path,
                           uniprot_homo_ensembl_gene_rdf_path,
                           uniprot_homo_ensembl_transcript_rdf_path,
                           uniprot_homo_rdf_path, uniprot_homo_refseq_rdf_path,
                           download_uniprot_fly_rdf,
                           download_uniprot_flybase_annotation_rdf,
                           download_uniprot_flybase_gene_rdf,
                           download_uniprot_flybase_transcript_rdf,
                           download_uniprot_homo_ensembl_gene_rdf,
                           download_uniprot_homo_ensembl_transcript_rdf,
                           download_uniprot_homo_rdf,
                           download_uniprot_homo_refseq_rdf,
                           set_uniprot_version)
from mirna.microcosm import (download_microcosm_targets, microcosm_rdf_path,
                             write_microcosm_rdf, microcosm_fly_mirs)
from mirna.mcneill import (write_mcneill_cns_screen_rdf,
                           write_mcneill_muscle_screen_rdf,
                           mcneill_cns_screen_rdf_path,
                           mcneill_muscle_screen_rdf_path)
from mirna.mirbase import (download_mirbase_aliases, mirbase_rdf_path,
                           write_mirbase_rdf)
from mirna.targetscan import (download_targetscan, targetscan_rdf_path,
                              write_targetscan_rdf)
from mirna.affymetrix import (affymetrix_fly_annotations_rdf_path,
                              write_affymetrix_fly_annotations_rdf)
from mirna.flybase import (flybase_rdf_path, write_flybase_rdf)
from mirna.ibanez import (ibanez_five_prime_mir_homologs_rdf_path,
                          ibanez_fullseq_mir_homologs_rdf_path,
                          write_ibanez_five_prime_mir_homologs_rdf,
                          write_ibanez_fullseq_mir_homologs_rdf)
from mirna.roundup import (download_roundup_from_sparql, roundup_rdf_path,
                           roundup_sparql_rdf_path,
                           write_roundup_orthologs_rdf)
from mirna.synaptomedb import (synaptomedb_v1_rdf_path, write_synaptomedb_rdf)
from mirna.nmjrnai import (get_nmj_rnai_genes,
                           get_nmj_rnai_gain_of_function_genes)



VIRT7_DB = 'virt7'
STARDOG_DB = 'stardog'
# DB_TYPE = VIRT7_DB
DB_TYPE = STARDOG_DB
virt7 = virtuoso.Virtuoso7(secrets.virtuoso_dba_password,
                           config.virtuoso_load_dir)
virt_sparq = sparql.Sparql(virt7.sparql_endpoint()) # used for querying


# Structure of a Dmel Transcript Annotation ID
# http://flybase.org/static_pages/docs/nomenclature/nomenclature3.html#2.4.
# Annotation IDs are represented in a common way: a species-specific 2
# letter prefix followed by a four or five digit integer. For historical
# reasons, there are two 2-letter prefixes for D. melanogaster: CG for
# protein-coding genes and CR for non-protein-coding-genes.
# And from personal communication with Josh Goodman, "The '-RX' and '-PX'
# suffixes (X being one or more letters) are given to transcripts and
# polypeptides respectively."
#

# species names used for db tables and microcosm files
mus = 'mus_musculus'
dro = 'drosophila_melanogaster'
cel = 'caenorhabditis_elegans'
homo = 'homo_sapiens'

# mirna target prediction databases
MICROCOSM = 'microcosm'
TARGETSCAN = 'targetscan'
targets_dbs = [MICROCOSM, TARGETSCAN]
MERGED = 'merged' # both prediction databases

# The seven mirs that were screened in muscle and CNS tissue to examine the
# differential gene expression effects.
orginal_screened_mirs = ['dme-miR-34', 'dme-miR-92b', 'dme-miR-137',
                         'dme-miR-190', 'dme-miR-219', 'dme-miR-276a',
                         'dme-miR-277']
# The current mirbase ids for the seven mirs screened by McNeill and Van Vactor.
screened_mirs = ['dme-miR-34-5p', 'dme-miR-92b-3p', 'dme-miR-137-3p',
                 'dme-miR-190-5p', 'dme-miR-219-5p', 'dme-miR-276a-3p',
                 'dme-miR-277-3p']

muscle_tissue = mcneill.MUSCLE
cns_tissue = mcneill.CNS
TISSUES = [muscle_tissue, cns_tissue]

# The twenty-six (formerly twenty-seven before the non-existent dme-miR-953 was
# dropped) mirs that were functionally validated by McNeill # and Van Vactor,
# by looking for morphological changes in the muscle phenotype.
original_validated_mirs = [
    'dme-miR-8', 'dme-miR-13a', 'dme-miR-14', 'dme-miR-34', 'dme-miR-92a',
    'dme-miR-92b', 'dme-miR-190', 'dme-miR-137', 'dme-miR-219', 'dme-miR-276a',
    'dme-miR-277', 'dme-miR-279', 'dme-miR-287', 'dme-miR-304', 'dme-miR-308',
    'dme-miR-313', 'dme-miR-314', 'dme-miR-316', 'dme-miR-932', 'dme-miR-969',
    'dme-miR-970', 'dme-miR-978', 'dme-miR-979', 'dme-miR-982', 'dme-miR-999',
    'dme-miR-1014']
# Twenty-six mirs have current miRBase ids.
validated_mirs = [
    u'dme-miR-8-3p', u'dme-miR-13a-3p', u'dme-miR-14-3p', u'dme-miR-34-5p',
    u'dme-miR-92a-3p', u'dme-miR-92b-3p', u'dme-miR-190-5p', u'dme-miR-137-3p',
    u'dme-miR-219-5p', u'dme-miR-276a-3p', u'dme-miR-277-3p',
    u'dme-miR-279-3p', u'dme-miR-287-3p', u'dme-miR-304-5p', u'dme-miR-308-3p',
    u'dme-miR-313-3p', u'dme-miR-314-3p', u'dme-miR-316-5p', u'dme-miR-932-5p',
    u'dme-miR-969-5p', u'dme-miR-970-3p', u'dme-miR-978-3p', u'dme-miR-979-3p',
    u'dme-miR-982-5p', u'dme-miR-999-3p', u'dme-miR-1014-3p']

# Copied from a list of 147? mirs emailed from Elizabeth McNeill (most of
# mirbase drosophila mirnas from a few years before 2013?)
original_147_mirs = [
    'bantam', 'let-7', 'mir-1', 'mir-2a', 'mir-2b', 'mir-2c', 'mir-3', 'mir-4',
    'mir-5', 'mir-6', 'mir-7', 'mir-8', 'mir-9a', 'mir-9b', 'mir-9c', 'mir-10',
    'mir-11', 'mir-12', 'mir-13a', 'mir-13b', 'mir-14', 'mir-31a', 'mir-31b',
    'mir-33', 'mir-34', 'mir-79', 'mir-87', 'mir-92a', 'mir-92b', 'mir-100',
    'mir-124', 'mir-125', 'mir-133', 'mir-137', 'mir-184', 'mir-190',
    'mir-193', 'mir-210', 'mir-219', 'mir-252', 'mir-263a', 'mir-263b',
    'mir-274', 'mir-275', 'mir-276', 'mir-276a', 'mir-276b', 'mir-277',
    'mir-278', 'mir-279', 'mir-281', 'mir-2811', 'mir-2812', 'mir-282',
    'mir-283', 'mir-284', 'mir-285', 'mir-286', 'mir-287', 'mir-288',
    'mir-289', 'mir-303', 'mir-304', 'mir-305', 'mir-306', 'mir-307',
    'mir-308', 'mir-309', 'mir-310', 'mir-311', 'mir-312', 'mir-313',
    'mir-314', 'mir-315', 'mir-316', 'mir-317', 'mir-318', 'mir-375',
    'mir-927', 'mir-929', 'mir-932', 'mir-954', 'mir-955', 'mir-956',
    'mir-957', 'mir-958', 'mir-959', 'mir-960', 'mir-961', 'mir-962',
    'mir-963', 'mir-964', 'mir-965', 'mir-966', 'mir-967', 'mir-968',
    'mir-969', 'mir-970', 'mir-971', 'mir-972', 'mir-973', 'mir-974',
    'mir-975', 'mir-976', 'mir-977', 'mir-978', 'mir-979', 'mir-980',
    'mir-981', 'mir-982', 'mir-983', 'mir-984', 'mir-985', 'mir-986',
    'mir-987', 'mir-988', 'mir-989', 'mir-990', 'mir-991', 'mir-992',
    'mir-993', 'mir-994', 'mir-995', 'mir-996', 'mir-997', 'mir-998',
    'mir-999', 'mir-1000', 'mir-1001', 'mir-1002', 'mir-1003', 'mir-1004',
    'mir-1005', 'mir-1006', 'mir-1007', 'mir-1008', 'mir-1009', 'mir-1010',
    'mir-1011', 'mir-1012', 'mir-1013', 'mir-1014', 'mir-1015', 'mir-1016',
    'mir-1017', 'mir-iab-4-3p', 'mir-iab-4-5p', ]

# The only original mirbase ids that could not be mapped to a current mirbase
# id were mir-276, mir-2811 and mir-2812.  mir-2811 is a dead id, mir-2812 is
# an id for bro-mir-2812 (i.e. wrong species), and mir-276 legitimately looked
# like it should have been searched for under dme-miR-276a and dme-miR-276b.
# However, since both of those mirs were independently searched for, no loss of
# information occurred by excluding mir-276.
# There 144 out of 147 that were mapped to a current mirbase id.
original_147_updated_data = [
    ['bantam', u'dme-bantam', u'dme-bantam-3p'],
    ['let-7', u'dme-let-7', u'dme-let-7-5p'],
    ['mir-1', u'dme-miR-1', u'dme-miR-1-3p'],
    ['mir-2a', u'dme-miR-2a', u'dme-miR-2a-3p'],
    ['mir-2b', u'dme-miR-2b', u'dme-miR-2b-3p'],
    ['mir-2c', u'dme-miR-2c', u'dme-miR-2c-3p'],
    ['mir-3', u'dme-miR-3', u'dme-miR-3-3p'],
    ['mir-4', u'dme-miR-4', u'dme-miR-4-3p'],
    ['mir-5', u'dme-miR-5', u'dme-miR-5-5p'],
    ['mir-6', u'dme-miR-6', u'dme-miR-6-3p'],
    ['mir-7', u'dme-miR-7', u'dme-miR-7-5p'],
    ['mir-8', u'dme-miR-8', u'dme-miR-8-3p'],
    ['mir-9a', u'dme-miR-9a', u'dme-miR-9a-5p'],
    ['mir-9b', u'dme-miR-9b', u'dme-miR-9b-5p'],
    ['mir-9c', u'dme-miR-9c', u'dme-miR-9c-5p'],
    ['mir-10', u'dme-miR-10', u'dme-miR-10-5p'],
    ['mir-11', u'dme-miR-11', u'dme-miR-11-3p'],
    ['mir-12', u'dme-miR-12', u'dme-miR-12-5p'],
    ['mir-13a', u'dme-miR-13a', u'dme-miR-13a-3p'],
    ['mir-13b', u'dme-miR-13b', u'dme-miR-13b-3p'],
    ['mir-14', u'dme-miR-14', u'dme-miR-14-3p'],
    ['mir-31a', u'dme-miR-31a', u'dme-miR-31a-5p'],
    ['mir-31b', u'dme-miR-31b', u'dme-miR-31b-5p'],
    ['mir-33', u'dme-miR-33', u'dme-miR-33-5p'],
    ['mir-34', u'dme-miR-34', u'dme-miR-34-5p'],
    ['mir-79', u'dme-miR-79', u'dme-miR-79-3p'],
    ['mir-87', u'dme-miR-87', u'dme-miR-87-3p'],
    ['mir-92a', u'dme-miR-92a', u'dme-miR-92a-3p'],
    ['mir-92b', u'dme-miR-92b', u'dme-miR-92b-3p'],
    ['mir-100', u'dme-miR-100', u'dme-miR-100-5p'],
    ['mir-124', u'dme-miR-124', u'dme-miR-124-3p'],
    ['mir-125', u'dme-miR-125', u'dme-miR-125-5p'],
    ['mir-133', u'dme-miR-133', u'dme-miR-133-3p'],
    ['mir-137', u'dme-miR-137', u'dme-miR-137-3p'],
    ['mir-184', u'dme-miR-184', u'dme-miR-184-3p'],
    ['mir-190', u'dme-miR-190', u'dme-miR-190-5p'],
    ['mir-193', u'dme-miR-193', u'dme-miR-193-3p'],
    ['mir-210', u'dme-miR-210', u'dme-miR-210-3p'],
    ['mir-219', u'dme-miR-219', u'dme-miR-219-5p'],
    ['mir-252', u'dme-miR-252', u'dme-miR-252-5p'],
    ['mir-263a', u'dme-miR-263a', u'dme-miR-263a-5p'],
    ['mir-263b', u'dme-miR-263b', u'dme-miR-263b-5p'],
    ['mir-274', u'dme-miR-274', u'dme-miR-274-5p'],
    ['mir-275', u'dme-miR-275', u'dme-miR-275-3p'],
    ['mir-276', u'dme-miR-276', None],
    ['mir-276a', u'dme-miR-276a', u'dme-miR-276a-3p'],
    ['mir-276b', u'dme-miR-276b', u'dme-miR-276b-3p'],
    ['mir-277', u'dme-miR-277', u'dme-miR-277-3p'],
    ['mir-278', u'dme-miR-278', u'dme-miR-278-3p'],
    ['mir-279', u'dme-miR-279', u'dme-miR-279-3p'],
    ['mir-281', u'dme-miR-281', u'dme-miR-281-3p'],
    ['mir-2811', u'dme-miR-2811', None],
    ['mir-2812', u'dme-miR-2812', None],
    ['mir-282', u'dme-miR-282', u'dme-miR-282-5p'],
    ['mir-283', u'dme-miR-283', u'dme-miR-283-5p'],
    ['mir-284', u'dme-miR-284', u'dme-miR-284-3p'],
    ['mir-285', u'dme-miR-285', u'dme-miR-285-3p'],
    ['mir-286', u'dme-miR-286', u'dme-miR-286-3p'],
    ['mir-287', u'dme-miR-287', u'dme-miR-287-3p'],
    ['mir-288', u'dme-miR-288', u'dme-miR-288-3p'],
    ['mir-289', u'dme-miR-289', u'dme-miR-289-5p'],
    ['mir-303', u'dme-miR-303', u'dme-miR-303-5p'],
    ['mir-304', u'dme-miR-304', u'dme-miR-304-5p'],
    ['mir-305', u'dme-miR-305', u'dme-miR-305-5p'],
    ['mir-306', u'dme-miR-306', u'dme-miR-306-5p'],
    ['mir-307', u'dme-miR-307', u'dme-miR-307a-3p'],
    ['mir-308', u'dme-miR-308', u'dme-miR-308-3p'],
    ['mir-309', u'dme-miR-309', u'dme-miR-309-3p'],
    ['mir-310', u'dme-miR-310', u'dme-miR-310-3p'],
    ['mir-311', u'dme-miR-311', u'dme-miR-311-3p'],
    ['mir-312', u'dme-miR-312', u'dme-miR-312-3p'],
    ['mir-313', u'dme-miR-313', u'dme-miR-313-3p'],
    ['mir-314', u'dme-miR-314', u'dme-miR-314-3p'],
    ['mir-315', u'dme-miR-315', u'dme-miR-315-5p'],
    ['mir-316', u'dme-miR-316', u'dme-miR-316-5p'],
    ['mir-317', u'dme-miR-317', u'dme-miR-317-3p'],
    ['mir-318', u'dme-miR-318', u'dme-miR-318-3p'],
    ['mir-375', u'dme-miR-375', u'dme-miR-375-3p'],
    ['mir-927', u'dme-miR-927', u'dme-miR-927-5p'],
    ['mir-929', u'dme-miR-929', u'dme-miR-929-3p'],
    ['mir-932', u'dme-miR-932', u'dme-miR-932-5p'],
    ['mir-954', u'dme-miR-954', u'dme-miR-954-5p'],
    ['mir-955', u'dme-miR-955', u'dme-miR-955-5p'],
    ['mir-956', u'dme-miR-956', u'dme-miR-956-3p'],
    ['mir-957', u'dme-miR-957', u'dme-miR-957-3p'],
    ['mir-958', u'dme-miR-958', u'dme-miR-958-3p'],
    ['mir-959', u'dme-miR-959', u'dme-miR-959-3p'],
    ['mir-960', u'dme-miR-960', u'dme-miR-960-5p'],
    ['mir-961', u'dme-miR-961', u'dme-miR-961-5p'],
    ['mir-962', u'dme-miR-962', u'dme-miR-962-5p'],
    ['mir-963', u'dme-miR-963', u'dme-miR-963-5p'],
    ['mir-964', u'dme-miR-964', u'dme-miR-964-5p'],
    ['mir-965', u'dme-miR-965', u'dme-miR-965-3p'],
    ['mir-966', u'dme-miR-966', u'dme-miR-966-5p'],
    ['mir-967', u'dme-miR-967', u'dme-miR-967-5p'],
    ['mir-968', u'dme-miR-968', u'dme-miR-968-5p'],
    ['mir-969', u'dme-miR-969', u'dme-miR-969-5p'],
    ['mir-970', u'dme-miR-970', u'dme-miR-970-3p'],
    ['mir-971', u'dme-miR-971', u'dme-miR-971-3p'],
    ['mir-972', u'dme-miR-972', u'dme-miR-972-3p'],
    ['mir-973', u'dme-miR-973', u'dme-miR-973-5p'],
    ['mir-974', u'dme-miR-974', u'dme-miR-974-5p'],
    ['mir-975', u'dme-miR-975', u'dme-miR-975-5p'],
    ['mir-976', u'dme-miR-976', u'dme-miR-976-3p'],
    ['mir-977', u'dme-miR-977', u'dme-miR-977-3p'],
    ['mir-978', u'dme-miR-978', u'dme-miR-978-3p'],
    ['mir-979', u'dme-miR-979', u'dme-miR-979-3p'],
    ['mir-980', u'dme-miR-980', u'dme-miR-980-3p'],
    ['mir-981', u'dme-miR-981', u'dme-miR-981-3p'],
    ['mir-982', u'dme-miR-982', u'dme-miR-982-5p'],
    ['mir-983', u'dme-miR-983', u'dme-miR-983-5p'],
    ['mir-984', u'dme-miR-984', u'dme-miR-984-5p'],
    ['mir-985', u'dme-miR-985', u'dme-miR-985-3p'],
    ['mir-986', u'dme-miR-986', u'dme-miR-986-5p'],
    ['mir-987', u'dme-miR-987', u'dme-miR-987-5p'],
    ['mir-988', u'dme-miR-988', u'dme-miR-988-3p'],
    ['mir-989', u'dme-miR-989', u'dme-miR-989-3p'],
    ['mir-990', u'dme-miR-990', u'dme-miR-990-5p'],
    ['mir-991', u'dme-miR-991', u'dme-miR-991-3p'],
    ['mir-992', u'dme-miR-992', u'dme-miR-992-3p'],
    ['mir-993', u'dme-miR-993', u'dme-miR-993-3p'],
    ['mir-994', u'dme-miR-994', u'dme-miR-994-5p'],
    ['mir-995', u'dme-miR-995', u'dme-miR-995-3p'],
    ['mir-996', u'dme-miR-996', u'dme-miR-996-3p'],
    ['mir-997', u'dme-miR-997', u'dme-miR-997-5p'],
    ['mir-998', u'dme-miR-998', u'dme-miR-998-3p'],
    ['mir-999', u'dme-miR-999', u'dme-miR-999-3p'],
    ['mir-1000', u'dme-miR-1000', u'dme-miR-1000-5p'],
    ['mir-1001', u'dme-miR-1001', u'dme-miR-1001-5p'],
    ['mir-1002', u'dme-miR-1002', u'dme-miR-1002-5p'],
    ['mir-1003', u'dme-miR-1003', u'dme-miR-1003-3p'],
    ['mir-1004', u'dme-miR-1004', u'dme-miR-1004-3p'],
    ['mir-1005', u'dme-miR-1005', u'dme-miR-1005-3p'],
    ['mir-1006', u'dme-miR-1006', u'dme-miR-1006-3p'],
    ['mir-1007', u'dme-miR-1007', u'dme-miR-1007-3p'],
    ['mir-1008', u'dme-miR-1008', u'dme-miR-1008-3p'],
    ['mir-1009', u'dme-miR-1009', u'dme-miR-1009-3p'],
    ['mir-1010', u'dme-miR-1010', u'dme-miR-1010-3p'],
    ['mir-1011', u'dme-miR-1011', u'dme-miR-1011-3p'],
    ['mir-1012', u'dme-miR-1012', u'dme-miR-1012-3p'],
    ['mir-1013', u'dme-miR-1013', u'dme-miR-1013-3p'],
    ['mir-1014', u'dme-miR-1014', u'dme-miR-1014-3p'],
    ['mir-1015', u'dme-miR-1015', u'dme-miR-1015-3p'],
    ['mir-1016', u'dme-miR-1016', u'dme-miR-1016-3p'],
    ['mir-1017', u'dme-miR-1017', u'dme-miR-1017-3p'],
    ['mir-iab-4-3p', u'dme-miR-iab-4-3p', u'dme-miR-iab-4-3p'],
    ['mir-iab-4-5p', u'dme-miR-iab-4-5p', u'dme-miR-iab-4-5p'],
]

# The updated mirbase id is in the third column.
# the_147_mirs = [d[2] for d in original_147_updated_data if d[2] is not None]
# print the_147_mirs
# print len(the_147_mirs)
# 144
# the_147_mirs should really be the_144_mirs.
the_147_mirs = [
    u'dme-bantam-3p', u'dme-let-7-5p', u'dme-miR-1-3p', u'dme-miR-2a-3p',
    u'dme-miR-2b-3p', u'dme-miR-2c-3p', u'dme-miR-3-3p', u'dme-miR-4-3p',
    u'dme-miR-5-5p', u'dme-miR-6-3p', u'dme-miR-7-5p', u'dme-miR-8-3p',
    u'dme-miR-9a-5p', u'dme-miR-9b-5p', u'dme-miR-9c-5p', u'dme-miR-10-5p',
    u'dme-miR-11-3p', u'dme-miR-12-5p', u'dme-miR-13a-3p', u'dme-miR-13b-3p',
    u'dme-miR-14-3p', u'dme-miR-31a-5p', u'dme-miR-31b-5p', u'dme-miR-33-5p',
    u'dme-miR-34-5p', u'dme-miR-79-3p', u'dme-miR-87-3p', u'dme-miR-92a-3p',
    u'dme-miR-92b-3p', u'dme-miR-100-5p', u'dme-miR-124-3p', u'dme-miR-125-5p',
    u'dme-miR-133-3p', u'dme-miR-137-3p', u'dme-miR-184-3p', u'dme-miR-190-5p',
    u'dme-miR-193-3p', u'dme-miR-210-3p', u'dme-miR-219-5p', u'dme-miR-252-5p',
    u'dme-miR-263a-5p', u'dme-miR-263b-5p', u'dme-miR-274-5p',
    u'dme-miR-275-3p', u'dme-miR-276a-3p', u'dme-miR-276b-3p',
    u'dme-miR-277-3p', u'dme-miR-278-3p', u'dme-miR-279-3p', u'dme-miR-281-3p',
    u'dme-miR-282-5p', u'dme-miR-283-5p', u'dme-miR-284-3p', u'dme-miR-285-3p',
    u'dme-miR-286-3p', u'dme-miR-287-3p', u'dme-miR-288-3p', u'dme-miR-289-5p',
    u'dme-miR-303-5p', u'dme-miR-304-5p', u'dme-miR-305-5p', u'dme-miR-306-5p',
    u'dme-miR-307a-3p', u'dme-miR-308-3p', u'dme-miR-309-3p',
    u'dme-miR-310-3p', u'dme-miR-311-3p', u'dme-miR-312-3p', u'dme-miR-313-3p',
    u'dme-miR-314-3p', u'dme-miR-315-5p', u'dme-miR-316-5p', u'dme-miR-317-3p',
    u'dme-miR-318-3p', u'dme-miR-375-3p', u'dme-miR-927-5p', u'dme-miR-929-3p',
    u'dme-miR-932-5p', u'dme-miR-954-5p', u'dme-miR-955-5p', u'dme-miR-956-3p',
    u'dme-miR-957-3p', u'dme-miR-958-3p', u'dme-miR-959-3p', u'dme-miR-960-5p',
    u'dme-miR-961-5p', u'dme-miR-962-5p', u'dme-miR-963-5p', u'dme-miR-964-5p',
    u'dme-miR-965-3p', u'dme-miR-966-5p', u'dme-miR-967-5p', u'dme-miR-968-5p',
    u'dme-miR-969-5p', u'dme-miR-970-3p', u'dme-miR-971-3p', u'dme-miR-972-3p',
    u'dme-miR-973-5p', u'dme-miR-974-5p', u'dme-miR-975-5p', u'dme-miR-976-3p',
    u'dme-miR-977-3p', u'dme-miR-978-3p', u'dme-miR-979-3p', u'dme-miR-980-3p',
    u'dme-miR-981-3p', u'dme-miR-982-5p', u'dme-miR-983-5p', u'dme-miR-984-5p',
    u'dme-miR-985-3p', u'dme-miR-986-5p', u'dme-miR-987-5p', u'dme-miR-988-3p',
    u'dme-miR-989-3p', u'dme-miR-990-5p', u'dme-miR-991-3p', u'dme-miR-992-3p',
    u'dme-miR-993-3p', u'dme-miR-994-5p', u'dme-miR-995-3p', u'dme-miR-996-3p',
    u'dme-miR-997-5p', u'dme-miR-998-3p', u'dme-miR-999-3p',
    u'dme-miR-1000-5p', u'dme-miR-1001-5p', u'dme-miR-1002-5p',
    u'dme-miR-1003-3p', u'dme-miR-1004-3p', u'dme-miR-1005-3p',
    u'dme-miR-1006-3p', u'dme-miR-1007-3p', u'dme-miR-1008-3p',
    u'dme-miR-1009-3p', u'dme-miR-1010-3p', u'dme-miR-1011-3p',
    u'dme-miR-1012-3p', u'dme-miR-1013-3p', u'dme-miR-1014-3p',
    u'dme-miR-1015-3p', u'dme-miR-1016-3p', u'dme-miR-1017-3p',
    u'dme-miR-iab-4-3p', u'dme-miR-iab-4-5p']

# May 2014 Davie grant miRs
# "for the grant this week, I am mainly interested in the specific targets for five miRs"
# These are the updated mirbase ids for the 5 (actually 7) mirs davie sent
# in the email.  This list was generated by copying and pasting the
# corresponding updated id from `the_147_mirs`.
may_mirs = [u'dme-let-7-5p', u'dme-miR-13a-3p', u'dme-miR-13b-3p',
            u'dme-miR-34-5p', u'dme-miR-92a-3p', u'dme-miR-92b-3p',
            u'dme-miR-137-3p']


# miRNAs with a NMJ phenotype (i.e. in the validated_mirs) AND conserved miR in
# human and fly AND whose targets rank highly in overlap with conserved
# synaptic genes.
# This list comes from file Elizabeth emailed on 2014/2/14, now located at
# deploy/vanvactor_mirna/data/20140214_Hit_target_ranking.docx
original_validated_conserved_and_ranked_mirs = ['miR-8-3p', 'miR-13a-3p',
                                                'miR-14-3p', 'miR-34-5p',
                                                'miR-92a-3p', 'miR-92b-3p',
                                                'miR-137-3p', 'miR-190-5p',
                                                'miR-219-5p', 'miR-276a-3p',
                                                'miR-277-3p', 'miR-279-3p',
                                                'miR-304-5p', 'miR-308-3p',
                                                'miR-314-3p', 'miR-316-5p',
                                                'miR-1014-3p']
validated_conserved_and_ranked_mirs = ['dme-' + m for m in original_validated_conserved_and_ranked_mirs]
# sanity check: these mirs are a subset of the validated mirs.
assert not set(validated_conserved_and_ranked_mirs) - set(validated_mirs)


# Fake Named Graph URIs used to:
# - select which tissue or ibanez homology method to use in queries
# - organize loading and reloading of data into the database.
affymetrix_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_affymetrix')
flybase_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_flybase')
microcosm_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_microcosm')
mirbase_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_mirbase')
roundup_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_roundup')
synaptomedb_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_synaptomedb')
targetscan_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_targetscan')
uniprot_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_uniprot')
mcneill_cns_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_mcneill_cns')
mcneill_muscle_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_mcneill_muscle')
ibanez_five_prime_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_ibanez_five_prime')
ibanez_full_sequence_ng = rdflib.URIRef('http://purl.example.com/graph/mirna_ibanez_full_sequence')


# taxons and pairs used for roundup orthologs
taxmap = {dro: dro_taxon, homo: homo_taxon, mus: mus_taxon, cel: cel_taxon}
TAXON_TO_NAME = {cel_taxon: cel, dro_taxon: dro, homo_taxon: homo, mus_taxon: mus}
# Orthologs between human and the other 3 species.
ROUNDUP_PAIRS = [(dro_taxon, homo_taxon), (mus_taxon, homo_taxon), (cel_taxon, homo_taxon)]

IBANEZ_METHOD_TO_NG = {ibanez.FIVE_PRIME: ibanez_five_prime_ng,
                       ibanez.FULL_SEQ: ibanez_full_sequence_ng}
ibanez_methods = [ibanez.FIVE_PRIME, ibanez.FULL_SEQ]
TARGETS_DB_TO_NG = {MICROCOSM: microcosm_ng, TARGETSCAN: targetscan_ng}


##############################
# RDF GRAPH DATABASE FUNCTIONS


def stop_stardog():
    user = secrets.stardog_admin_user
    password = secrets.stardog_admin_password
    print 'To stop stardog, run the following command:'
    print 'stardog-admin server stop --username {} --passwd {}'.format(user, password)


def start_stardog():
    # In order to use Stardog in federated queries, the authentication needs to
    # be disabled for the default anonymous user.
    # http://www.stardog.com/docs/faq/.
    # Disable authentication for the `anonymous` user
    user = secrets.stardog_admin_user
    password = secrets.stardog_admin_password
    print 'To start stardog, run the following command:'
    print 'stardog-admin --disable-security server start --username {} --passwd {}'.format(user, password)


def stardog_sparql():
    return sparql.Sparql(
        'http://localhost:5822/{db}/query'.format(db=stardog_database_name()),
        auth=(secrets.stardog_user, secrets.stardog_password))


def stardog_construct_query(query, deduplicate=True):
    '''
    Run a construct query and return (deduplicated) triples as N-triples.

    :param query: a CONSTRUCT query string.
    :param deduplicate: Stardog does not necessarily return deduplicated
    triples from a construct query.  If deduplicate is True (the default),
    duplicate triples will be removed.
    '''
    print 'stardog_construct_query'
    print query
    out = stardog_sparql().query(query, accept='text/plain')
    if deduplicate:
        uniques = set(line for line in cStringIO.StringIO(out) if line.strip())
        out = ''.join(uniques)
    return out


def stardog_json_query(query):
    print 'stardog_json_query'
    print query
    try:
        return stardog_sparql().query(query,
                                      accept='application/sparql-results+json')
    except requests.exceptions.HTTPError as e:
        print e
        print dir(e)
        print vars(e)
        print vars(e.response)
        print e.response.content
        raise


def virtuoso_json_query(query):
    print 'virtuoso_json_query'
    print query
    return virt_sparq.query(query, accept='application/sparql-results+json')


def sparql_json_query(query):
    if DB_TYPE == VIRT7_DB:
        return virtuoso_json_query(query)
    elif DB_TYPE == STARDOG_DB:
        return stardog_json_query(query)
    else:
        raise Exception('Unrecognized DB_TYPE', DB_TYPE)


def query_for_id_tuples(query, bindings):
    '''
    Make a sparql query, parse out the bindings (with no type conversion) from
    the results for the given bindings.  The binding values should be IRIs of
    the form "http://example.com/path/to/id".  Return a list of tuples of ids
    (not IRIs), where the ids are ordered according to `bindings`.
    '''
    result = sparql_json_query(query)
    # IRIs are like u'http://purl.targetscan.org/mir_family/miR-33'
    # or u'http://purl.targetscan.org/mir_family/miR-279%2F286%2F996'
    id_tuples = []
    for row in result['results']['bindings']:
        iris = [row[binding]['value'] for binding in bindings]
        ids = tuple(iri_to_id(iri) for iri in iris)
        id_tuples.append(ids)
    return id_tuples


def query_for_ids(query, binding):
    '''
    Wrap a sparql query into a call to stardog, and parse out the id value of 
    the results for the given binding.  The binding values should be IRIs of
    the form "http://example.com/path/to/id".  Return a list of ids.
    '''
    result = sparql_json_query(query)
    # IRIs are like u'http://purl.targetscan.org/mir_family/miR-33'
    # or u'http://purl.targetscan.org/mir_family/miR-279%2F286%2F996'
    iris = [b[binding]['value'] for b in result['results']['bindings']]
    return iris_to_ids(iris)



def iris_to_ids(iris, token='/'):
    '''
    Convert a list of URIs/IRIs into plain old identifiers.  The ids are also
    URL unquoted, which should be safe, considering that, as IRIs, they should
    be quoted, right?  Or is that just a URL thing?
    iris: a list of IRIs in the form: 'http://example.com/path/to/identifier'
    token: a character (or string) to split the identifier from the rest of
    the IRI.  In RDF, it is typically a slash ('/') or a hash ('#').
    '''
    return [iri_to_id(iri, token) for iri in iris]

def iri_to_id(iri, token='/'):
    '''
    Convert a URIs/IRIs into a plain old identifier.  The id are also
    URL unquoted, which should be safe, considering that, as IRIs, they should
    be quoted, right?  Or is that just a URL thing?
    iri: str.  An IRI in the form: 'http://example.com/path/to/identifier'
    token: a character (or string) to split the identifier from the rest of
    the IRI.  In RDF, it is typically a slash ('/') or a hash ('#').
    '''
    return urllib.unquote(iri.rsplit(token, 1)[-1])


def prefixes():
    '''
    Return a string containing several common PREFIX clauses used in the
    mirna rdf database.
    '''
    return '''
    PREFIX syndb:<http://purl.synaptomedb.org/owl/>
    PREFIX mb:<http://purl.mirbase.org/owl/>
    PREFIX ts:<http://purl.targetscan.org/owl/>
    PREFIX ex:<http://purl.example.com/owl/>
    PREFIX obo:<http://purl.org/obo/owl/obo#>
    PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
    PREFIX up:<http://purl.uniprot.org/core/>
    PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
    PREFIX db:<http://purl.example.com/database/>
    '''


def stardog_database_name():
    '''
    The name of the database for the vanvactor mirna collaboration
    '''
    return 'mirna'


def drop_stardog_database():
    cmd = 'time stardog-admin db drop --username {} --passwd {} {}'.format(
        secrets.stardog_admin_user, secrets.stardog_admin_password,
        stardog_database_name())
    call(cmd)


def stardog_drop_graph(graph):
    cmd = 'time stardog data remove --username {} --passwd {}'.format(
        secrets.stardog_admin_user, secrets.stardog_admin_password)
    cmd += ' --named-graph {} {}'.format(graph, stardog_database_name())
    call(cmd)


def stardog_load_graph(filename, graph=None):
    '''
    :param graph: A named graph URI.  If None, use the default graph.
    '''
    cmd = 'time stardog data add --username {} --passwd {}'.format(
        secrets.stardog_admin_user, secrets.stardog_admin_password)
    if graph:
        cmd += ' --named-graph {}'.format(graph)
    cmd += ' {} {}'.format(stardog_database_name(), filename)
    call(cmd)


def load_stardog_database():
    '''
    Data is loaded when the database is created b/c loading is much faster
    that way, according to the stardog docs.
    '''
    cmd = 'time stardog-admin db create --username {} --passwd {}'.format(
        secrets.stardog_admin_user, secrets.stardog_admin_password)
    cmd += ' --name {}'.format(stardog_database_name())
    # Quickly bulk load data at database creation time.
    # Bulk load quads if different graphs are required.
    # cmd += ' ' + ' '.join(all_rdf_paths())
    call(cmd)

    # Load triples into different named graphs.
    for filename, graph in all_rdf_files_and_graphs():
        stardog_load_graph(filename, graph)


def stardog_default_graph_rdf_files_and_graphs():
    return [
        (affymetrix_fly_annotations_rdf_path(), affymetrix_ng),
        (flybase_rdf_path(), flybase_ng),
        (microcosm_rdf_path(), microcosm_ng),
        (mirbase_rdf_path(), mirbase_ng),
        (synaptomedb_v1_rdf_path(), synaptomedb_ng),
        (targetscan_rdf_path(), targetscan_ng),
        # UniProt from the idmapping.dat file
        # (uniprot_rdf_path(), uniprot_ng),
        # UniProt from beta.sparql.uniprot.org
        (uniprot_homo_ensembl_gene_rdf_path(), uniprot_ng),
        (uniprot_homo_ensembl_transcript_rdf_path(), uniprot_ng),
        (uniprot_homo_refseq_rdf_path(), uniprot_ng),
        (uniprot_flybase_transcript_rdf_path(), uniprot_ng),
        (uniprot_flybase_gene_rdf_path(), uniprot_ng),
        (uniprot_flybase_annotation_rdf_path(), uniprot_ng),
        (uniprot_homo_rdf_path(), uniprot_ng),
        (uniprot_fly_rdf_path(), uniprot_ng),
        # Roundup from the downloaded txt files
        (roundup_rdf_path(), roundup_ng),
        # Roundup from sparql.roundup.hms.harvard.edu
        (roundup_sparql_rdf_path(), roundup_ng),
    ]


def stardog_named_graph_rdf_files_and_graphs():
    return [
        (ibanez_five_prime_mir_homologs_rdf_path(), ibanez_five_prime_ng),
        (ibanez_fullseq_mir_homologs_rdf_path(), ibanez_full_sequence_ng),
        (mcneill_cns_screen_rdf_path(), mcneill_cns_ng),
        (mcneill_muscle_screen_rdf_path(), mcneill_muscle_ng),
    ]


def all_rdf_files_and_graphs():
    return [
        (affymetrix_fly_annotations_rdf_path(), affymetrix_ng),
        (flybase_rdf_path(), flybase_ng),
        (ibanez_five_prime_mir_homologs_rdf_path(), ibanez_five_prime_ng),
        (ibanez_fullseq_mir_homologs_rdf_path(), ibanez_full_sequence_ng),
        (mcneill_cns_screen_rdf_path(), mcneill_cns_ng),
        (mcneill_muscle_screen_rdf_path(), mcneill_muscle_ng),
        (microcosm_rdf_path(), microcosm_ng),
        (mirbase_rdf_path(), mirbase_ng),
        (synaptomedb_v1_rdf_path(), synaptomedb_ng),
        (targetscan_rdf_path(), targetscan_ng),
        # UniProt from the idmapping.dat file
        # (uniprot_rdf_path(), uniprot_ng),
        # UniProt from beta.sparql.uniprot.org
        (uniprot_homo_ensembl_gene_rdf_path(), uniprot_ng),
        (uniprot_homo_ensembl_transcript_rdf_path(), uniprot_ng),
        (uniprot_homo_refseq_rdf_path(), uniprot_ng),
        (uniprot_flybase_transcript_rdf_path(), uniprot_ng),
        (uniprot_flybase_gene_rdf_path(), uniprot_ng),
        (uniprot_flybase_annotation_rdf_path(), uniprot_ng),
        (uniprot_homo_rdf_path(), uniprot_ng),
        (uniprot_fly_rdf_path(), uniprot_ng),
        # Roundup from the downloaded txt files
        (roundup_rdf_path(), roundup_ng),
        # Roundup from sparql.roundup.hms.harvard.edu
        (roundup_sparql_rdf_path(), roundup_ng),
    ]


def all_rdf_paths():
    return [filename for filename, graph in all_rdf_files_and_graphs()]


def load_virtuoso_database():
    for filename, graph in all_rdf_files_and_graphs():
        # drop any existing graph
        virt7.drop_graph(graph, silent=True)
        # load files into graph, one at a time
        print 'loading graph {}'.format(graph)
        virt7.load_graph(graph, filename)


def load_rdf_database():
    if DB_TYPE == VIRT7_DB:
        load_virtuoso_database()
    elif DB_TYPE == STARDOG_DB:
        drop_stardog_database()
        load_stardog_database()
    else:
        raise Exception('Unrecognized DB_TYPE', DB_TYPE)


###################
# CONSTRUCT QUERIES
#
# Construct, save, and load edges that simplify the relationships in the
# data and speed queries.


def constructed_path(basename):
    return os.path.join(config.datadir, 'constructed', basename)


def count_lines(text):
    count = 0
    for line in cStringIO.StringIO(text):
        count += 1
    return count


def construct_and_save(query, basename, approx_num_triples):
    filename = constructed_path(basename)
    nt = stardog_construct_query(query)
    count = count_lines(nt)
    print '{count} lines in CONSTRUCT result.'.format(count=count)
    assert count >= approx_num_triples
    with open(filename, 'w') as fh:
        fh.write(nt)


def construct_save_and_load(query, basename, graph, approx_num_triples):
    filename = constructed_path(basename)
    construct_and_save(query, filename)
    # stardog_load_graph(filename, graph)


def constructed_microcosm_fly_path():
    fn = 'fly_microcosm-v5_mirbase_id_to_flybase_gene.nt'
    return constructed_path(fn)


def constructed_microcosm_human_path():
    fn = 'human_microcosm-v5_mirbase_id_to_ensembl_gene.nt'
    return constructed_path(fn)


def constructed_targetscan_fly_path():
    fn = 'fly_targetscan-6.2_mirbase_id_to_flybase_gene.nt'
    return constructed_path(fn)


def constructed_targetscan_human_path():
    fn = 'human_targetscan-6.2_mirbase_id_to_ensembl_gene.nt'
    return constructed_path(fn)


def constructed_roundup_human_fly_orthologs_path():
    fn = 'roundup-4_human_ensembl_gene_to_flybase_gene.nt'
    return constructed_path(fn)


def constructed_ibanez_five_prime_path():
    fn = 'ibanez-2008-five_prime-current-mirbase-homologs.nt'
    return constructed_path(fn)


def constructed_ibanez_fullseq_path():
    fn = 'ibanez-2008-fullseq-current-mirbase-homologs.nt'
    return constructed_path(fn)


def construct_targetscan_fly_mirna_id_to_gene_edges():
    '''
    Construct edges from current mirbase id to flybase gene
    by going from current mirbase id to mirbase accession to
    targetscan mir family to flybase gene.
    '''
    query = prefixes() + '''
    CONSTRUCT { ?cmid mb:targets ?fg . }
    WHERE {
    # fly flybase genes
    ?fg up:organism taxon:7227 .
    ?fg up:database db:flybase_gene .
    # targetscan mir family targeting fly genes
    ?tmf mb:targets ?fg .
    ?tmf up:database db:targetscan_mir_family .
    # convert targetscan mir family to mature mirbase accs
    ?tmf rdfs:seeAlso ?macc .
    ?macc up:database db:mirbase_acc .
    ?macc a mb:mature_mirna .
    ?macc up:organism taxon:7227 .
    # convert mirbase accs to current mirbase ids
    ?macc mb:current_id ?cmid .
    ?cmid up:database db:mirbase_id .
    }
    '''
    fn = constructed_targetscan_fly_path()
    approx_num_triples = 15000
    construct_and_save(query, fn, approx_num_triples)


def construct_targetscan_human_mirna_id_to_gene_edges():
    '''
    Construct edges from current mirbase id to ensembl gene
    by going from current mirbase id to mirbase accession to
    targetscan mir family to refseq transcript to uniprot to 
    ensembl gene.
    '''
    query = prefixes() + '''
    CONSTRUCT { ?cmid mb:targets ?heg . }
    WHERE {
    # human ensembl gene
    ?heg up:organism taxon:9606 .
    ?heg up:database db:ensembl_gene .
    # human uniprot
    ?hu rdfs:seeAlso ?heg .
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    # human refseq transcript
    ?hu rdfs:seeAlso ?ht .
    ?ht up:database db:ncbi_refseq .
    ?ht up:organism taxon:9606 .
    # targetscan mir family targeting human refseq transcript
    ?hmf mb:targets ?ht .
    ?hmf up:database db:targetscan_mir_family .
    # targetscan mir family to mature mirbase accs
    ?hmf rdfs:seeAlso ?macc .
    ?macc up:organism taxon:9606 .
    ?macc up:database db:mirbase_acc .
    ?macc a mb:mature_mirna .
    # mirbase acc to current mirbase id
    ?macc mb:current_id ?cmid .
    ?cmid up:database db:mirbase_id .
    }
    '''
    fn = constructed_targetscan_human_path()
    approx_num_triples = 170000
    construct_and_save(query, fn, approx_num_triples)


def construct_microcosm_fly_mirna_id_to_gene_edges():
    '''
    Construct edges from current mirbase id to flybase gene
    by going from current mirbase id to mirbase accession to
    mirbase id to flybase annotation id to flybase transcript
    to uniprot to flybase gene
    '''
    query = prefixes() + '''
    CONSTRUCT { ?cmid mb:targets ?fbg . }
    WHERE {
    # fly mirbase id
    ?mid up:database db:mirbase_id .
    ?mid up:organism taxon:7227 .
    # mirbase id to mature mirbase accession
    ?macc mb:has_id ?mid .
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .
    # mirbase acc to current mirbase id
    ?macc mb:current_id ?cmid .
    ?cmid up:database db:mirbase_id .
    ?cmid up:organism taxon:7227 .
    # mirbase id to targeted flybase annotation ids
    ?mid mb:targets ?fba .
    ?fba up:database db:flybase_annotation .
    # flybase annotation id to flybase transcript
    ?fbt rdfs:seeAlso ?fba .
    ?fbt up:database db:flybase_transcript .
    # flybase transcript to uniprot
    ?u rdfs:seeAlso ?fbt .
    ?u up:database db:uniprot .
    # uniprot to flybase gene .
    ?u rdfs:seeAlso ?fbg .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    }
    '''
    fn = constructed_microcosm_fly_path()
    approx_num_triples = 15000
    construct_and_save(query, fn, approx_num_triples)


def construct_microcosm_human_mirna_id_to_gene_edges():
    '''
    Construct edges from current mirbase id to ensembl gene
    by going from current mirbase id to mirbase accession to
    mirbase id to ensembl transcript to uniprot to ensembl gene.
    '''
    query = prefixes() + '''
    CONSTRUCT { ?cmid mb:targets ?heg . }
    WHERE {
    # current human mirbase id to mirbase acc
    ?macc mb:current_id ?cmid .
    ?cmid up:database db:mirbase_id .
    ?cmid up:organism taxon:9606 .
    # mirbase id to mature mirbase accession
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .
    # mirbase id
    ?macc mb:has_id ?mid .
    ?mid up:database db:mirbase_id .
    # mirbase id to targeted ensembl transcripts
    ?mid mb:targets ?enst .
    ?enst up:database db:ensembl_transcript .
    # ensembl transcript to uniprot
    ?u rdfs:seeAlso ?enst .
    ?u up:database db:uniprot .
    # uniprot to human ensembl gene .
    ?u rdfs:seeAlso ?heg .
    ?heg up:database db:ensembl_gene .
    ?heg up:organism taxon:9606 .
    }
    '''
    fn = constructed_microcosm_human_path()
    approx_num_triples = 28000
    construct_and_save(query, fn, approx_num_triples)



def construct_roundup_human_gene_to_fly_gene_edges():
    '''
    Construct edges from human ensembl genes to flybase genes by going
    from ensembl gene to human uniprot to fly uniprot to flybase gene.
    '''
    query = prefixes() + '''
    CONSTRUCT {
    ?heg obo:orthologous_to ?fbg .
    ?fbg obo:orthologous_to ?heg .
    }
    WHERE {
    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .
    # uniprot to human ensembl gene .
    ?hu rdfs:seeAlso ?heg .
    ?heg up:database db:ensembl_gene .
    ?heg up:organism taxon:9606 .
    # uniprot to flybase gene .
    ?du rdfs:seeAlso ?fbg .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    }
    '''
    fn = constructed_roundup_human_fly_orthologs_path()
    approx_num_triples = 10000
    construct_and_save(query, fn, approx_num_triples)


def construct_ibanez_fullseq_edges():
    query = prefixes() + '''
    CONSTRUCT {
    ?dcmid obo:homologous_to ?hcmid .
    }
    WHERE {
    # homologous mirbase ids
    GRAPH <''' + ibanez_full_sequence_ng + '''> {
        ?dmid obo:homologous_to ?hmid .
    }
    # current mirbase id for fly
    ?dmid up:database db:mirbase_id .
    ?dmacc mb:has_id ?dmid .
    ?dmacc a mb:mature_mirna .
    ?dmacc up:database db:mirbase_acc .
    ?dmacc mb:current_id ?dcmid .
    ?dcmid up:database db:mirbase_id .
    ?dcmid up:organism taxon:7227 .
    # current mirbase id for human
    ?hmid up:database db:mirbase_id .
    ?hmacc mb:has_id ?hmid .
    ?hmacc a mb:mature_mirna .
    ?hmacc up:database db:mirbase_acc .
    ?hmacc mb:current_id ?hcmid .
    ?hcmid up:database db:mirbase_id .
    ?hcmid up:organism taxon:9606 .
    }
    '''
    fn = constructed_ibanez_fullseq_path()
    approx_num_triples = 50
    construct_and_save(query, fn, approx_num_triples)


def construct_ibanez_five_prime_edges():
    query = prefixes() + '''
    CONSTRUCT {
    ?dcmid obo:homologous_to ?hcmid .
    }
    WHERE {
    # homologous mirbase ids
    GRAPH <''' + ibanez_five_prime_ng + '''> {
        ?dmid obo:homologous_to ?hmid .
    }
    # current mirbase id for fly
    ?dmid up:database db:mirbase_id .
    ?dmacc mb:has_id ?dmid .
    ?dmacc a mb:mature_mirna .
    ?dmacc up:database db:mirbase_acc .
    ?dmacc mb:current_id ?dcmid .
    ?dcmid up:database db:mirbase_id .
    ?dcmid up:organism taxon:7227 .
    # current mirbase id for human
    ?hmid up:database db:mirbase_id .
    ?hmacc mb:has_id ?hmid .
    ?hmacc a mb:mature_mirna .
    ?hmacc up:database db:mirbase_acc .
    ?hmacc mb:current_id ?hcmid .
    ?hcmid up:database db:mirbase_id .
    ?hcmid up:organism taxon:9606 .
    }
    '''
    fn = constructed_ibanez_five_prime_path()
    approx_num_triples = 50
    construct_and_save(query, fn, approx_num_triples)


def load_constructed_edges():
    pairs = [(constructed_microcosm_fly_path(), microcosm_ng),
             (constructed_microcosm_human_path(), microcosm_ng),
             (constructed_targetscan_fly_path(), targetscan_ng),
             (constructed_targetscan_human_path(), targetscan_ng),
             (constructed_roundup_human_fly_orthologs_path(), roundup_ng),
             (constructed_ibanez_five_prime_path(), ibanez_five_prime_ng),
             (constructed_ibanez_fullseq_path(), ibanez_full_sequence_ng),
            ]
    for filename, graph in pairs:
        stardog_load_graph(filename, graph)


################
# SPARQL QUERIES

def update_mirbase_ids(mirbase_ids):
    '''
    Take a list of mirbase ids for MATURE mirbase accessions, and return a dict
    that maps each id to the current mirbase id for the accession the mirbase
    id maps to.  Possible data quirks.  If the mirbase id does not map, it will
    not be in the returned dict.  If the mirbase id maps to multiple mature
    mirbase accs (I do not *think* this happens), the returned dict will only
    map the original to one of the current mirbase ids for those mirbase accs.
    Finally, the original mirbase id could be the same as the current
    mirbase id if it is up-to-date.
    '''
    # Generate URIs for the mirbase ids
    # Empirical tests showed that the FULL Head HTTP error occurs for 25+
    # mirbase ids, but not for 20 or fewer.  So run queries in groups of 20.

    old2new = {}
    all_mids = [mirbase_id_iri(mirbase_id) for mirbase_id in mirbase_ids]
    for mids in util.groupsOfN(all_mids, 20):
        # Query for a mapping from original id to current id
        query = prefixes() + '''
        SELECT DISTINCT ?dm_old ?dm
        WHERE {
        VALUES ?dm_old { ''' + ' '.join(['<{}>'.format(mid) for mid in mids]) + ''' }
        # convert old mirbase ids to mirbase accs
        ?dma mb:has_id ?dm_old .
        # FILTER ?dm_old IN ( ''' + ', '.join('<{}>'.format(mid) for mid in mids) + ''' )
        ?dm_old up:database db:mirbase_id .
        # convert mirbase accs to current mirbase ids
        ?dma up:database db:mirbase_acc .
        ?dma a mb:mature_mirna .
        ?dma mb:current_id ?dm .
        }
        '''
        print 'query:', query
        result = sparql_json_query(query)
        # iris are like u'http://purl.targetscan.org/mir_family/miR-33'
        # or u'http://purl.targetscan.org/mir_family/miR-279%2F286%2F996'
        lookup = dict((os.path.basename(b['dm_old']['value']), 
                    os.path.basename(b['dm']['value'])) for b in
                    result['results']['bindings'])
        old2new.update(lookup)
    return old2new # map old id to new id


def original_147_mirs_searchname(mir):
    '''
    Change strings like mir-3 into dme-miR-3 when searching for the most recent
    version.
    '''
    if not mir.startswith('dme-'):
        mir = u'dme-' + mir
    if mir.find('-mir-') != -1:
        mir = re.sub('-mir-', '-miR-', mir)
    return mir


def update_original_147_mirs():
    orig2search = {mir: original_147_mirs_searchname(mir) for mir in original_147_mirs}
    search_mirs = orig2search.values()
    print 'The ids and search terms:'
    for mir in original_147_mirs:
        print mir, '===>', orig2search[mir]
    print 'Updating mirbase ids'
    search2latest = update_mirbase_ids(search_mirs) # [:20])  # 20 <=x 25 < 30 < 50 < 80
    for mir in original_147_mirs:
        search = orig2search[mir]
        latest = search2latest.get(search)
        print repr([mir, search, latest])


def update_validated_mirs_to_current_mirbase_ids():
    '''
    Used to update a list of mirs from Elizabeth to their current mirbase ids,
    so they can be more easily used in the database queries
    '''
    lookup = update_mirbase_ids(original_validated_mirs)
    for mir in original_validated_mirs:
        print mir, lookup.get(mir)

    current_mirs = [lookup[mir] for mir in original_validated_mirs if mir in lookup]
    print current_mirs
    print len(current_mirs)
    return current_mirs


def merged_fly_mir_targets_query_sub(mirbase_id, where):
    current_mirbase_id = update_mirbase_ids([mirbase_id])[mirbase_id]
    mir = mirbase_id_iri(current_mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?fbg
    WHERE {
    ''' + where(mir) + '''
    }
    '''
    return query_for_ids(query, 'fbg')


def merged_fly_mir_targets_where(mir):
    return '''
    # flybase genes targeted by this mirbase id
    ?dmid mb:targets ?fbg .
    FILTER ( ?dmid = <''' + mir + '''> )
    ?dmid up:database db:mirbase_id .
    ?dmid up:organism taxon:7227 .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    '''


def merged_conserved_fly_mir_targets_where(mir):
    return merged_fly_mir_targets_where(mir) + '''
    # conserved in human
    ?heg obo:orthologous_to ?fbg .
    ?heg up:database db:ensembl_gene .
    ?heg up:organism taxon:9606 .
    '''


def merged_conserved_synaptic_fly_mir_targets_where(mir):
    return merged_conserved_fly_mir_targets_where(mir) + '''
    # synaptic gene
    ?heg up:classifiedWith syndb:synaptic .
    '''


def merged_conserved_synaptic_fly_targets(mirbase_id):
    '''
    Return the flybase genes that are conserved in human and targeted by
    `mirbase_id` in the predicted targets databases.
    '''
    return merged_fly_mir_targets_query_sub(
        mirbase_id, merged_conserved_synaptic_fly_mir_targets_where)


def merged_conserved_fly_targets(mirbase_id):
    '''
    Return the flybase genes that are conserved in human and targeted by
    `mirbase_id` in the predicted targets databases.
    '''
    return merged_fly_mir_targets_query_sub(
        mirbase_id, merged_conserved_fly_mir_targets_where)


def merged_fly_targets(mirbase_id):
    '''
    Return the flybase genes targeted by `mirbase_id` in the predicted targets
    databases.
    '''
    return merged_fly_mir_targets_query_sub(
        mirbase_id, merged_fly_mir_targets_where)


def conserved_genes_where():
    return '''
    # conserved genes
    ?heg up:database db:ensembl_gene .
    ?heg up:organism taxon:9606 .
    ?heg obo:orthologous_to ?fbg .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    '''


def conserved_synaptic_genes_where():
    return conserved_genes_where() + '''
    # synaptic genes
    ?heg up:classifiedWith syndb:synaptic .
    '''


def short_human_mirs_targeting_conserved_synapse_genes(targets_db):
    '''
    Return a list of human mirs that are predicted to target genes that are
    synaptic genes (according the synapsedb) and orthologous to fly genes.
    '''
    targets_ng = TARGETS_DB_TO_NG[targets_db]
    query = prefixes() + '''
    SELECT DISTINCT ?cmid
    WHERE {
    ''' + conserved_synaptic_genes_where() + '''
    # current human mirbase ids targeting ensembl genes
    GRAPH <''' + targets_ng + '''>
    {
    ?cmid mb:targets ?heg .
    }
    ?cmid up:database db:mirbase_id .
    ?cmid up:organism taxon:9606 .
    ?macc mb:current_id ?cmid .
    }
    '''
    return query_for_ids(query, 'cmid')



def short_fly_mirs_targeting_conserved_synapse_genes(targets_db):
    '''
    Return a list of fly mirs that are predicted to target genes that are
    orthologous to human synaptic genes (according the synapsedb).
    '''
    targets_ng = TARGETS_DB_TO_NG[targets_db]
    query = prefixes() + '''
    SELECT DISTINCT ?cmid
    WHERE {
    ''' + conserved_synaptic_genes_where() + '''
    # current fly mirbase ids targeting flybase genes
    GRAPH <''' + targets_ng + '''>
    {
    ?cmid mb:targets ?fbg .
    }
    ?cmid up:database db:mirbase_id .
    ?cmid up:organism taxon:7227 .
    ?macc mb:current_id ?cmid .
    }
    '''
    return query_for_ids(query, 'cmid')


def short_mirs_targeting_conserved_synapse_genes(targets_db, species):
    '''
    Return a list of mirs that are predicted to target genes that are
    conserved synaptic genes (according the synaptomedb and roundup).
    '''
    taxon = taxmap[species]
    query = prefixes() + '''
    SELECT DISTINCT ?cmid
    WHERE {
    ''' + conserved_synaptic_genes_where() + '''
    ''' + targets_in_graph_where(targets_db, species) + '''
    ?cmid up:database db:mirbase_id .
    ?macc mb:current_id ?cmid .
    ?cmid up:organism taxon:''' + taxon + ''' .
    }
    '''
    return query_for_ids(query, 'cmid')


def targets_in_graph_where(targets_db, species):
    '''
    Return a graph pattern that matches `?cmid mb:targets ?heg .` or 
    `?cmid mb:tagets ?fbg .`, possibly wrapped in a named graph stanza.
    If targets_db is MICROCOSM, use the microcosm graph, if it is TARGETSCAN,
    use the targetscan graph, and if it is MERGED, do not wrap clause in
    a graph.
    e.g.
    GRAPH <http://example.com/my_named_graph>
    {
    ?cmid mb:targets ?fbg .
    }

    species: dro or homo.
    targets_db: MICROCOSM, TARGETSCAN, or MERGED.
    '''
    if species == dro:
        gene_var = '?fbg'
    elif species == homo:
        gene_var = '?heg'
    else:
        raise Exception('Unknown species.', species)
    cmid_clause = '''
    # current mirbase ids targeting genes
    ?cmid mb:targets ''' + gene_var + ''' .
    '''

    if targets_db == MERGED:
        return cmid_clause
    elif targets_db in targets_dbs:
        targets_ng = TARGETS_DB_TO_NG[targets_db]
        where = 'GRAPH <' + targets_ng + '> {\n'
        where += cmid_clause
        where += '}\n'
        return where
    else:
        raise Exception('Unknown targets_db.', targets_db)


def merged_human_mirs_targeting_conserved_synapse_genes():
    '''
    Return human mirbase ids predicted to target a conserved synapse gene.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?hmid
    WHERE {
    # human mirbase ids targeting human ensembl genes
    ?hmid mb:targets ?heg .
    ?hmid up:database db:mirbase_id .
    ?hmid up:organism taxon:9606 .
    ''' + conserved_synaptic_genes_where() + '''
    }
    '''
    return query_for_ids(query, 'hmid')


def merged_fly_mirs_targeting_conserved_synapse_genes():
    '''
    Return fly mirbase ids predicted to target a conserved synapse gene.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?dmid
    WHERE {
    # fly mirbase ids targeting flybase genes
    ?dmid mb:targets ?fbg .
    ?dmid up:database db:mirbase_id .
    ?dmid up:organism taxon:7227 .
    ''' + conserved_synaptic_genes_where() + '''
    }
    '''
    return query_for_ids(query, 'dmid')


def short_human_to_fly_conserved_synapse_genes_table():
    '''
    Return a list of tuples of human ensembl gene and flybase gene, where
    the human gene is a synaptic gene (i.e. in synaptome db) and the fly
    gene is orthologous to the human gene.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?heg ?fbg
    WHERE {
    # conserved synaptic genes
    ?heg up:classifiedWith syndb:synaptic .
    ?heg up:organism taxon:9606 .
    ?heg up:database db:ensembl_gene .
    ?heg obo:orthologous_to ?fbg .
    ?fbg up:organism taxon:7227 .
    ?fbg up:database db:flybase_gene .
    }
    '''
    return query_for_id_tuples(query, ['heg', 'fbg'])


def fly_conserved_synapse_genes_human_orthologs_map():
    '''
    Return a mapping from fly genes that are conserved synapse genes to their
    human orthologs (that are synaptic genes).
    '''
    human_to_fly_table = short_human_to_fly_conserved_synapse_genes_table()
    fly_to_human_genes = collections.defaultdict(set)
    for human_gene, fly_gene in human_to_fly_table:
        fly_to_human_genes[fly_gene].add(human_gene)

    return fly_to_human_genes


def short_human_to_fly_conserved_genes_table():
    '''
    Return a list of tuples of human ensembl gene and flybase gene that are
    orthologous to each other.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?heg ?fbg
    WHERE {
    # conserved genes
    ?heg up:organism taxon:9606 .
    ?heg up:database db:ensembl_gene .
    ?heg obo:orthologous_to ?fbg .
    ?fbg up:organism taxon:7227 .
    ?fbg up:database db:flybase_gene .
    }
    '''
    return query_for_id_tuples(query, ['heg', 'fbg'])


def fly_conserved_genes():
    human, fly = zip(*short_human_to_fly_conserved_genes_table())
    return list(set(fly))


def fly_genes():
    '''
    return a list of flybase genes.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?fbg
    WHERE {
    # drosophila melanogaster flybase genes
    ?fbg up:organism taxon:7227 .
    ?fbg up:database db:flybase_gene .
    }
    '''
    return query_for_ids(query, 'fbg')


def human_to_fly_conserved_synapse_genes_table():
    '''
    Return a list of tuples of human ensembl gene and flybase gene.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?heg ?fbg
    WHERE {
    # human synaptic ensembl genes
    ?heg up:classifiedWith syndb:synaptic .
    ?heg up:organism taxon:9606 .
    ?heg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?heg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # flybase genes
    ?du rdfs:seeAlso ?fbg .
    ?fbg up:database db:flybase_gene .
    }
    '''
    return query_for_id_tuples(query, ['heg', 'fbg'])


def human_conserved_synapse_genes():
    '''
    Synapse genes are human genes from SynapseDB.
    Conserved synapse genes are the subset that have fly orthologs.
    Return a list of human ensembl synaptic genes with fly orthologs.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?heg
    WHERE {
    # human synaptic ensembl genes
    ?heg up:classifiedWith syndb:synaptic .
    ?heg up:organism taxon:9606 .
    ?heg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?heg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # flybase genes
    ?du rdfs:seeAlso ?dg .
    ?dg up:database db:flybase_gene .
    }
    '''
    return query_for_ids(query, 'heg')


def fly_conserved_synapse_genes():
    '''
    Synapse genes are human genes from SynapseDB.
    Conserved synapse genes are the subset that have fly orthologs.
    Return a list of fly genes that are orthologous to human synaptic genes.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?dg
    WHERE {
    # human synaptic ensembl genes
    ?heg up:classifiedWith syndb:synaptic .
    ?heg up:organism taxon:9606 .
    ?heg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?heg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # flybase genes
    ?du rdfs:seeAlso ?dg .
    ?dg up:database db:flybase_gene .
    }
    '''
    return query_for_ids(query, 'dg')


def targetscan_human_mirs_targeting_conserved_synapse_genes():
    '''
    Synapse genes are human genes from SynapseDB.
    Conserved synapse genes are the subset that have fly orthologs.
    Return the human mirs predicted by targetscan to target one or more of
    the conserved synapse genes.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?hm
    WHERE {
    # human synaptic genes
    ?hg up:classifiedWith syndb:synaptic .
    ?hg up:organism taxon:9606 .
    ?hu rdfs:seeAlso ?hg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # human refseq transcripts
    ?hu rdfs:seeAlso ?ht .
    ?ht up:database db:ncbi_refseq .
    ?ht up:organism taxon:9606 .

    # targetscan mirs targeting human refseq transcripts
    ?hmf mb:targets ?ht .
    ?hmf up:database db:targetscan_mir_family .

    # convert targetscan mir family into current mirbase accs
    ?hmf rdfs:seeAlso ?hma .
    ?hma up:organism taxon:9606 .

    # convert mirbase accs into current mirbase ids
    ?hma up:database db:mirbase_acc .
    ?hma a mb:mature_mirna .
    ?hma mb:current_id ?hm .
    }
    '''
    return query_for_ids(query, 'hm')


def targetscan_fly_mirs_targeting_conserved_synapse_genes():
    '''
    Synapse genes are human genes from SynapseDB.
    Conserved synapse genes are the subset that have fly orthologs.
    Return the fly mirs that map to targetscan mir families predicted by
    targetscan to target one or more of fly orthologs of the conserved synapse
    genes.

    fly targetscan family targeting flybase gene linked to uniprot orthologous
    to human uniprot linked to human ensembl gene classified with synaptic.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?dm
    WHERE {
    # human synaptic ensembl genes
    ?hg up:classifiedWith syndb:synaptic .
    ?hg up:organism taxon:9606 .
    ?hg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?hg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # fly genes
    ?du rdfs:seeAlso ?dg .
    ?dg up:database db:flybase_gene .
    # targetscan mir family targeting fly genes
    ?dmf mb:targets ?dg .
    ?dmf up:database db:targetscan_mir_family .
    # convert targetscan mir family to mirbase accs
    ?dmf rdfs:seeAlso ?dma .
    ?dma up:organism taxon:7227 .
    # convert mirbase accs to current mirbase ids
    ?dma up:database db:mirbase_acc .
    ?dma a mb:mature_mirna .
    ?dma mb:current_id ?dm .
    }
    '''
    return query_for_ids(query, 'dm')


def microcosm_human_mirs_targeting_conserved_synapse_genes():
    '''
    Return a list of human mirs that are predicted to target genes that are
    synaptic genes (according the synapsedb) and orthologous to fly genes.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?hm
    WHERE {

    # synaptic human ensembl genes
    ?hg up:classifiedWith syndb:synaptic .
    ?hg up:organism taxon:9606 .
    ?hg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?hg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # human ensembl transcripts
    ?hu rdfs:seeAlso ?ht .
    ?ht up:database db:ensembl_transcript .
    ?ht up:organism taxon:9606 .

    # microcosm mirbase ids targeting transcripts
    ?hm_old mb:targets ?ht .

    # convert microcosm mirbase ids into mirbase accs
    ?hma mb:has_id ?hm_old .
    ?hma up:organism taxon:9606 .

    # convert mirbase accs into current mirbase ids
    ?hma up:database db:mirbase_acc .
    ?hma a mb:mature_mirna .
    ?hma mb:current_id ?hm .
    }
    '''
    return query_for_ids(query, 'hm')


def microcosm_fly_mirs_targeting_conserved_synapse_genes():
    '''
    fly targetscan family targeting flybase gene linked to uniprot orthologous
    to human uniprot linked to human ensembl gene classified with synaptic.
    '''
    query = prefixes() + '''
    SELECT DISTINCT ?dm
    WHERE {
    # human synaptic ensembl genes
    ?hg up:classifiedWith syndb:synaptic .
    ?hg up:organism taxon:9606 .
    ?hg up:database db:ensembl_gene .
    ?hu rdfs:seeAlso ?hg .

    # human - fly orthologs
    ?hu up:organism taxon:9606 .
    ?hu up:database db:uniprot .
    ?du obo:orthologous_to ?hu .
    ?du up:organism taxon:7227 .
    ?du up:database db:uniprot .

    # fly transcrips
    ?du rdfs:seeAlso ?dt .
    ?dt up:database db:flybase_transcript .

    # fly annotations
    ?dt rdfs:seeAlso ?da .
    ?da up:database db:flybase_annotation .

    # microcosm mirbase ids targeting fly annotation ids
    ?dm_old mb:targets ?da .

    # convert old mirbase ids to mirbase accs
    ?dma mb:has_id ?dm_old .
    ?dma up:database db:mirbase_acc .
    ?dma a mb:mature_mirna .

    # convert mirbase accs to current mirbase ids
    ?dma mb:current_id ?dm .
    ?dm up:database db:mirbase_id .
    }
    '''
    return query_for_ids(query, 'dm')


def microcosm_predicted_fly_mir_targets(mirbase_id):
    '''
    Given a mirbase id, go to a mature mirbase accession and back to mirbase
    ids (since we are not sure if the input mirbase_id is the same as one of
    the mirbase_ids used in microcosm), then to flybase accession targets and eventually to flybase genes.
    '''
    mir = mirbase_id_iri(mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?gene
    WHERE {{
    # original mirbase id to mature mirbase accession
    ?macc mb:has_id <{mid_orig}> .
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .

    # mirbase acc to all mirbase ids
    ?macc mb:has_id ?mid .
    ?mid up:database db:mirbase_id .

    # mirbase id to targeted flybase annotation ids
    ?mid mb:targets ?fba .
    ?fba up:database db:flybase_annotation .
    ?fba up:organism taxon:7227 .

    # flybase annotation id to flybase transcript
    ?fbt rdfs:seeAlso ?fba .
    ?fbt up:database db:flybase_transcript .

    # flybase transcript to uniprot
    ?u rdfs:seeAlso ?fbt .
    ?u up:database db:uniprot .

    # uniprot to flybase gene .
    ?u rdfs:seeAlso ?gene .
    ?gene up:database db:flybase_gene .
    ?gene up:organism taxon:7227 .
    }}
    '''.format(mid_orig=mir)
    return query_for_ids(query, 'gene')


def targetscan_predicted_fly_mir_targets(mirbase_id):
    '''
    Given a mirbase id, go to a mature mirbase accession, to a targetscan mir family,
    to predicted targeted flybase genes.
    '''
    mir = mirbase_id_iri(mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?gene
    WHERE {{
    # original mirbase id to mature mirbase accession
    ?macc mb:has_id <{mid_orig}> .
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .

    # mirbase acc to targetscan mir family
    ?tmf rdfs:seeAlso ?macc .
    ?tmf up:database db:targetscan_mir_family .

    # targetscan mir family to targeted flybase genes
    ?tmf mb:targets ?gene .
    ?gene up:database db:flybase_gene .
    ?gene up:organism taxon:7227 .
    }}
    '''.format(mid_orig=mir)
    return query_for_ids(query, 'gene')



def short_conserved_fly_targets_of_mirs(mirbase_ids, targets_db):
    return map_flatten_set_list(short_conserved_fly_targets, mirbase_ids, 
                                [targets_db] * len(mirbase_ids))



def short_conserved_fly_targets(mirbase_id, targets_db):
    '''
    Return the flybase genes that are conserved in human and targeted by
    `mirbase_id` in the predicted targets database `targets_db`.
    '''
    current_mirbase_id = update_mirbase_ids([mirbase_id])[mirbase_id]
    mir = mirbase_id_iri(current_mirbase_id)
    targets_ng = TARGETS_DB_TO_NG[targets_db]
    query = prefixes() + '''
    SELECT DISTINCT ?fbg
    WHERE {
    # flybase genes targeted by this mirbase id
    GRAPH <''' + targets_ng + '''> {
        ?dmid mb:targets ?fbg .
        }
    FILTER ( ?dmid = <''' + mir + '''> )
    ?dmid up:database db:mirbase_id .
    ?dmid up:organism taxon:7227 .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    # conserved in human
    ?heg obo:orthologous_to ?fbg .
    ?heg up:database db:ensembl_gene .
    ?heg up:organism taxon:9606 .
    }
    '''
    return query_for_ids(query, 'fbg')


def microcosm_fly_targets_table():
    '''
    Return a list of tuples, where each tuple contains an original mirbase id
    from microcosm, the current mirbase id, and a flybase gene targeted by
    the mirbase id.
    '''
    # original_mirs = microcosm_fly_mirs()
    # to_current_mir = update_mirbase_ids(original_mirs)
    # # to_original_mir = dict(value, key for to_current_mir.values())

    # table = []
    # for mir in original_mirs[:5]:
        # current_mir = to_current_mir.get(mir)
        # targets = sorted(short_fly_targets(mir, MICROCOSM))
        # if not targets:
            # table.append((mir, current_mir, None))
        # else:
            # for target in targets:
                # table.append((mir, current_mir, target))

    # return table

    query = prefixes() + '''
    SELECT DISTINCT ?dmid ?fbg
    WHERE {
    # flybase genes targeted by mirs in microcosm
    GRAPH <''' + TARGETS_DB_TO_NG[MICROCOSM] + '''> {
        ?dmid mb:targets ?fbg .
        }
    ?dmid up:database db:mirbase_id .
    ?dmid up:organism taxon:7227 .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    }
    '''
    table = query_for_id_tuples(query, ['dmid', 'fbg'])
    return sorted(set(table))




def fly_targets_set(mirbase_id, targets_db):
    '''
    mirbase_id: a fly mirbase id, e.g. dme-miR-8-3p
    Return the set of flybase genes that are targeted by the mirbase id, using
    the specified targets database, e.g. MICROCOSM, TARGETSCAN, or both (i.e.
    MERGED).
    '''
    if targets_db == MERGED:
        return set(short_fly_targets(mirbase_id, MICROCOSM)) | set(short_fly_targets(mirbase_id, TARGETSCAN))
    else:
        return set(short_fly_targets(mirbase_id, targets_db))



def short_fly_targets(mirbase_id, targets_db):
    '''
    Return the flybase genes targeted by `mirbase_id` in the predicted targets
    database `targets_db`.
    '''
    current_mirbase_id = update_mirbase_ids([mirbase_id])[mirbase_id]
    mir = mirbase_id_iri(current_mirbase_id)
    targets_ng = TARGETS_DB_TO_NG[targets_db]
    query = prefixes() + '''
    SELECT DISTINCT ?fbg
    WHERE {
    # flybase genes targeted by this mirbase id
    GRAPH <''' + targets_ng + '''> {
        ?dmid mb:targets ?fbg .
        }
    FILTER ( ?dmid = <''' + mir + '''> )
    ?dmid up:database db:mirbase_id .
    ?dmid up:organism taxon:7227 .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    }
    '''
    return query_for_ids(query, 'fbg')


def merged_fly_orthologs_of_targets_of_human_mirs_homologous_to_fly_mir(
    mirbase_id):
    '''
    Return flybase genes that are orthologous to human ensembl genes that
    are targeted by human mirbase ids that are homologous to fly `mirbase_id`.
    '''
    current_mirbase_id = update_mirbase_ids([mirbase_id])[mirbase_id]
    mir = mirbase_id_iri(current_mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?fbg
    WHERE {
    # get human mirbase ids homologous to this fly mirbase id
    ?dmid obo:homologous_to ?hmid .
    FILTER ( ?dmid = <''' + mir + '''> )
    ?dmid up:database db:mirbase_id .
    ?dmid up:organism taxon:7227 .
    ?hmid up:database db:mirbase_id .
    ?hmid up:organism taxon:9606 .
    # human mirbase id to ensembl gene
    ?hmid mb:targets ?heg .
    ?heg up:database db:ensembl_gene .
    ?heg up:organism taxon:9606 .
    # flybase gene orthologous to ensembl gene
    ?heg obo:orthologous_to ?fbg .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    }
    '''
    return query_for_ids(query, 'fbg')


def short_fly_orthologs_of_targets_of_human_mirs_homologous_to_fly_mir(
    mirbase_id, ibanez_method, targets_db):
    '''
    Return flybase genes that are orthologous to human ensembl genes that
    are targeted (according to `targets_db` by human mirbase ids that are
    homologous (according to `ibanez_method`) to fly `mirbase_id`.
    '''
    current_mirbase_id = update_mirbase_ids([mirbase_id])[mirbase_id]
    mir = mirbase_id_iri(current_mirbase_id)
    method_ng = IBANEZ_METHOD_TO_NG[ibanez_method]
    targets_ng = TARGETS_DB_TO_NG[targets_db]
    query = prefixes() + '''
    SELECT DISTINCT ?fbg
    WHERE {
    # get human mirbase ids homologous to this fly mirbase id
    GRAPH <''' + method_ng + '''> {
        ?dmid obo:homologous_to ?hmid .
    }
    FILTER ( ?dmid = <''' + mir + '''> )
    ?dmid up:database db:mirbase_id .
    ?dmid up:organism taxon:7227 .
    ?hmid up:database db:mirbase_id .
    ?hmid up:organism taxon:9606 .
    # human mirbase id to ensembl gene
    GRAPH <''' + targets_ng + '''> {
        ?hmid mb:targets ?heg .
        }
    ?heg up:database db:ensembl_gene .
    ?heg up:organism taxon:9606 .
    # flybase gene orthologous to ensembl gene
    ?heg obo:orthologous_to ?fbg .
    ?fbg up:database db:flybase_gene .
    ?fbg up:organism taxon:7227 .
    }
    '''
    return query_for_ids(query, 'fbg')


def affymetrix_to_flybase_table():
    query = prefixes() + '''
    SELECT DISTINCT ?affy ?fbg
    WHERE {
        ?affy up:database db:affymetrix_probeset .
        ?affy up:organism taxon:7227 .
        ?affy rdfs:seeAlso ?fbg .
        ?fbg up:database db:flybase_gene .
        ?fbg up:organism taxon:7227 .
    }
    '''
    return query_for_id_tuples(query, ['affy', 'fbg'])


def screened_fly_mir_targets_in_tissue(mirbase_id, tissue):

    # mcneill screen edges
    # graph.add((mir, regulates_pred, probeset))
    # graph.add((mir, db_pred, mirbase_id_db))
    # graph.add((mir, organism_pred, dro_iri))
    # graph.add((probeset, db_pred, affymetrix_probeset_db))
    # graph.add((probeset, organism_pred, dro_iri))

    # affy edges
    # graph.add((probeset, see_also_pred, gene))
    # graph.add((probeset, db_pred, affymetrix_probeset_db))
    # graph.add((probeset, organism_pred, dro_iri))
    # graph.add((gene, db_pred, flybase_gene_db))
    # graph.add((gene, organism_pred, dro_iri))

    if tissue == muscle_tissue:
        tissue_ng = mcneill_muscle_ng
    elif tissue == cns_tissue:
        tissue_ng = mcneill_cns_ng
    else:
        raise Exception('Unknown tissue {}. mirbase_id {}.'.format(
            tissue, mirbase_id))

    mir = mirbase_id_iri(mirbase_id)
    query = prefixes() + '''
    SELECT DISTINCT ?gene
    WHERE {{
    # original mirbase id to mature mirbase accession
    ?macc mb:has_id <{mid_orig}> .
    ?macc a mb:mature_mirna .
    ?macc up:database db:mirbase_acc .

    # mirbase acc to all mirbase ids
    ?macc mb:has_id ?mid .
    ?mid up:database db:mirbase_id .

    # mirbase id regulates an affymetrix probe
    GRAPH <{tissue_ng}> {{
        ?mid ex:regulates ?affy .
        ?affy up:database db:affymetrix_probeset .
    }}

    # affy probe to flybase gene
    ?affy rdfs:seeAlso ?gene .
    ?gene up:database db:flybase_gene .
    ?gene up:organism taxon:7227 .
    }}
    '''.format(tissue_ng=tissue_ng, mid_orig=mir)
    return query_for_ids(query, 'gene')


####################
# WRITE OUTPUT FILES



def results_dir():
    dn = datetime.datetime.now().strftime('%Y%m%d') + '_results'
    return makedirs(os.path.join(config.datadir, dn))


def write_set_overlap_file(set1, set2, name1, name2, filename):
    '''
    Given two sets, write out a csv file containing the set differences and 
    set intersection, placed in 3 columns, in order to be excel friendly.
    The first row is column headers based on name1 and name2.  The columns 
    are <name1>_not_<name2>, <name1>_and_<name2>, <name2>_not_<name1>.

    set1: a set.  Should have csv friendly values, since no quoting of any kind
    is done.
    set2: a set.
    name1: A name for set 1 used to build the column headers.  keep it simple.
    Avoid commas and special characters.
    name2: Like name 1, but for set 2.
    filename: where to write the csv file.
    '''
    set1 = set(set1)
    set2 = set(set2)
    one_not_two = sorted(set1 - set2)
    two_not_one = sorted(set2 - set1)
    one_and_two = sorted(set1 & set2)
    makedirs(os.path.dirname(filename))
    with open(filename, 'w') as fh:
        fh.write('{one}_not_{two},{one}_and_{two},{two}_not_{one}\n'.format(
            one=name1, two=name2))
        for row in itertools.izip_longest(one_not_two, one_and_two, two_not_one, fillvalue=''):
            fh.write(','.join([str(i) for i in row]) + '\n')


def write_list_file(items, filename):
    makedirs(os.path.dirname(filename))
    with open(filename, 'w') as fh:
        fh.write(''.join([item + '\n' for item in sorted(items)]))


def write_csv_file(rows, filename, headers=None):
    '''
    :param rows: list of rows, where each row is a list of value.
    :param headers: list of strings, used as a header row.  If None, no header
    row will be output.
    '''
    with open(filename, 'w') as fh:
        writer = csv.writer(fh)
        if headers:
            writer.writerow(headers)
        for row in rows:
            writer.writerow(row)


# PHASE 0
# Write useful list of genes and lookup tables.

def write_conserved_synaptic_genes():
    dn = os.path.join(results_dir(), 'conserved_synapse_genes')
    fn = os.path.join(dn, 'fly_conserved_synapse_genes.txt')
    genes = fly_conserved_synapse_genes()
    write_list_file(genes, fn)
    fn = os.path.join(dn, 'human_conserved_synapse_genes.txt')
    genes = human_conserved_synapse_genes()
    write_list_file(genes, fn)


def write_affymetrix_to_flybase_table():
    fn = os.path.join(results_dir(), 'affymetrix_to_flybase_genes.csv')
    write_csv_file(affymetrix_to_flybase_table(), fn,
                   headers=('affymetrix_probeset_id', 'flybase_gene'))


def write_human_to_fly_conserved_synapse_genes_table():
    dn = os.path.join(results_dir(), 'conserved_synapse_genes')
    fn = os.path.join(dn, 'human_to_fly_conserved_synapse_genes.csv')
    write_csv_file(human_to_fly_conserved_synapse_genes_table(), fn,
                   headers=('ensembl_gene', 'flybase_gene'))


# PHASE I
# Phase one focuses on miRs that target conserved synaptic genes.
# What are the conserved synaptic genes in fly and human?
# What are the microRNA targeting those genes in fly and human, for microcosm,
# targetscan, and merged?
# What are the overlaps of miRs for targetscan and microcosm?
# What is the overlap of the functionally validated fly mirs and the
# (microcosm, targetscan, merged) fly mirs targeting conserved synapse genes?


def phase_1_mirs_fn(method, species, dirname):
    basename = '{method}_predicted_{species}_mirs_'
    basename += 'targeting_conserved_synaptic_genes.txt'
    return os.path.join(dirname, basename.format(
        method=method, species=species))


def write_phase_1():
    dn = os.path.join(results_dir(), 'phase1')

    targ_hm = short_mirs_targeting_conserved_synapse_genes(TARGETSCAN, homo)
    targ_fm = short_mirs_targeting_conserved_synapse_genes(TARGETSCAN, dro)
    micr_hm = short_mirs_targeting_conserved_synapse_genes(MICROCOSM, homo)
    micr_fm = short_mirs_targeting_conserved_synapse_genes(MICROCOSM, dro)
    merg_hm = short_mirs_targeting_conserved_synapse_genes(MERGED, homo)
    merg_fm = short_mirs_targeting_conserved_synapse_genes(MERGED, dro)
    # functionally validated
    vali_fm = set(validated_mirs)

    assert set(merg_hm) == (set(targ_hm) | set(micr_hm))
    assert set(merg_fm) == (set(targ_fm) | set(micr_fm))

    write_list_file(targ_hm, phase_1_mirs_fn('targetscan', 'human', dn))
    write_list_file(targ_fm, phase_1_mirs_fn('targetscan', 'fly', dn))
    write_list_file(micr_hm, phase_1_mirs_fn('microcosm', 'human', dn))
    write_list_file(micr_fm, phase_1_mirs_fn('microcosm', 'fly', dn))
    write_list_file(merg_hm, phase_1_mirs_fn('merged', 'human', dn))
    write_list_file(merg_fm, phase_1_mirs_fn('merged', 'fly', dn))
    fn = os.path.join(dn, 'validated_fly_mirs.txt')
    write_list_file(vali_fm, fn)

    fn = os.path.join(dn, 'overlap_between_microcosm_and_targetscan_human' +
                      '_mirs_targeting_conserved_synaptic_genes.csv')
    write_set_overlap_file(micr_hm, targ_hm, 'microcosm', 'targetscan', fn)
    fn = os.path.join(dn, 'overlap_between_microcosm_and_targetscan_fly' +
                      '_mirs_targeting_conserved_synaptic_genes.csv')
    write_set_overlap_file(micr_fm, targ_fm, 'microcosm', 'targetscan', fn)
    fn = os.path.join(dn, 'overlap_between_validated_and_targetscan_fly' +
                      '_mirs_targeting_conserved_synaptic_genes.csv')
    write_set_overlap_file(vali_fm, targ_fm, 'validated', 'targetscan', fn)
    fn = os.path.join(dn, 'overlap_between_validated_and_microcosm_fly' +
                      '_mirs_targeting_conserved_synaptic_genes.csv')
    write_set_overlap_file(vali_fm, micr_fm, 'validated', 'microcosm', fn)
    fn = os.path.join(dn, 'overlap_between_validated_and_merged_fly' +
                      '_mirs_targeting_conserved_synaptic_genes.csv')
    write_set_overlap_file(vali_fm, merg_fm, 'validated', 'merged', fn)


def write_targetscan_human_mirs_targeting_conserved_synapse_genes():
    dn = os.path.join(results_dir(), 'phase1')
    mirs = sorted(targetscan_human_mirs_targeting_conserved_synapse_genes())
    fn = os.path.join(dn, 'predicted_targetscan_human_mirs_targeting_conserved_synapse_genes.txt')
    write_list_file(mirs, fn)


def write_microcosm_human_mirs_targeting_conserved_synapse_genes():
    dn = os.path.join(results_dir(), 'phase1')
    mirs = sorted(microcosm_human_mirs_targeting_conserved_synapse_genes())
    fn = os.path.join(dn, 'predicted_microcosm_human_mirs_targeting_conserved_synapse_genes.txt')
    write_list_file(mirs, fn)


def write_overlap_between_microcosm_human_mirs_and_targetscan_human_mirs():
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'overlap_between_microcosm_and_targetscan_human_mirs.csv')
    micro = set(microcosm_human_mirs_targeting_conserved_synapse_genes())
    targ = set(targetscan_human_mirs_targeting_conserved_synapse_genes())
    write_set_overlap_file(micro, targ, 'microcosm', 'targetscan', fn)


def write_targetscan_fly_mirs_targeting_conserved_synapse_genes():
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'predicted_targetscan_fly_mirs_targeting_conserved_synapse_genes.txt')
    mirs = sorted(targetscan_fly_mirs_targeting_conserved_synapse_genes())
    write_list_file(mirs, fn)


def write_microcosm_fly_mirs_targeting_conserved_synapse_genes():
    '''
    Write a list of fly mirs that are predicted to target genes that are
    orthologous to human synaptic genes (according the synapsedb).
    '''
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'predicted_microcosm_fly_mirs_targeting_conserved_synapse_genes.txt')
    mirs = sorted(microcosm_fly_mirs_targeting_conserved_synapse_genes())
    write_list_file(mirs, fn)


def write_overlap_between_microcosm_fly_mirs_and_targetscan_fly_mirs():
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'overlap_between_microcosm_and_targetscan_fly_mirs.csv')
    micro = set(microcosm_fly_mirs_targeting_conserved_synapse_genes())
    targ = set(targetscan_fly_mirs_targeting_conserved_synapse_genes())
    write_set_overlap_file(micro, targ, 'microcosm', 'targetscan', fn)


def write_overlap_between_validated_fly_mirs_and_microcosm_predicted_fly_mirs():
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'overlap_between_validated_and_predicted_microcosm_fly_mirs.csv')
    predicted = set(microcosm_fly_mirs_targeting_conserved_synapse_genes())
    validated = set(validated_mirs)
    write_set_overlap_file(validated, predicted, 'validated', 'predicted', fn)


def write_overlap_between_validated_fly_mirs_and_targetscan_predicted_fly_mirs():
    dn = os.path.join(results_dir(), 'phase1')
    fn = os.path.join(dn, 'overlap_between_validated_and_predicted_targetscan_fly_mirs.csv')
    predicted = set(targetscan_fly_mirs_targeting_conserved_synapse_genes())
    validated = set(validated_mirs)
    write_set_overlap_file(validated, predicted, 'validated', 'predicted', fn)


# PHASE II
# SETS and OVERLAPS for
# For EACH SCREENED MIR
#     MICROCOSM TARGETED FLY GENES (FBgn)
#     TARGETSCAN TARGETED FLY GENES (FBgn)
#     MERGED TARGETED FLY GENES (FBgn)
#     MICROCOSM TARGETED HUMAN GENES (ENSG)
#     TARGETSCAN TARGETED HUMAN GENES (ENSG)
#     MERGED TARGETED HUMAN GENES (ENSG)
#     For EACH TISSUE TYPE
#         SCREENED REGULATED FLY GENES (FBgn)
# - TARGETS OF 7 SCREENED MIRS in 2 TISSUES VS PREDICTED TARGETS of MICROCOSM and TARGETSCAN
#
# FOR EACH VALIDATED FLY MIR:
#     PREDICTED TARGETS OF THE FLY MIR COMPARED TO (UNION OF PREDICTED TARGETS OF THE HUMAN HOMOLOGS OF THE FLY MIR)


def write_overlap_between_screened_and_microcosm_and_targetscan_fly_mir_targets():
    dn = os.path.join(results_dir(), 'phase2', 'overlap_of_screened_and_predicted_target_genes')
    for mir in screened_mirs:
        microcosm_targets = set(microcosm_predicted_fly_mir_targets(mir))
        targetscan_targets = set(targetscan_predicted_fly_mir_targets(mir))
        mt_fn = os.path.join(dn, '{}_overlap_between_microcosm_and_targetscan_targets.csv'.format(mir))
        write_set_overlap_file(microcosm_targets, targetscan_targets, 'microcosm', 'targetscan', mt_fn)
        merg_fg = merged_fly_targets(mir)
        for tissue in TISSUES:
            screened_targets = set(screened_fly_mir_targets_in_tissue(mir, tissue))
            sm_fn = os.path.join(dn, '{}_{}_overlap_between_screened_and_microcosm_targets.csv'.format(mir, tissue))
            write_set_overlap_file(screened_targets, microcosm_targets, 'screened', 'microcosm', sm_fn)
            st_fn = os.path.join(dn, '{}_{}_overlap_between_screened_and_targetscan_targets.csv'.format(mir, tissue))
            write_set_overlap_file(screened_targets, targetscan_targets, 'screened', 'targetscan', st_fn)
            merg_fn = os.path.join(dn, '{}_{}_overlap_between_screened_and_merged_targets.csv'.format(mir, tissue))
            write_set_overlap_file(screened_targets, merg_fg, 'screened', 'merged', merg_fn)


def write_overlap_of_screened_fly_mir_targets_and_nmj_rnai_genes():
    dn = os.path.join(results_dir(), 'phase2', 'overlap_of_nmj_rnai_genes_and_screened_fly_targets')
    nmj_rnai = set(get_nmj_rnai_genes())
    for mir in screened_mirs:
        for tissue in TISSUES:
            screened_targets = set(screened_fly_mir_targets_in_tissue(mir, tissue))
            fn = os.path.join(dn, '{}_{}_overlap_between_nmj_rnai_genes_and_screened_fly_targets.csv'.format(mir, tissue))
            write_set_overlap_file(nmj_rnai, screened_targets, 'nmj_rnai_genes', 'screened', fn)


def write_overlap_of_validated_fly_mir_targets_and_nmj_rnai_genes():
    dn = os.path.join(results_dir(), 'phase2', 'overlap_of_nmj_rnai_genes_and_predicted_fly_targets')
    nmj_rnai = set(get_nmj_rnai_genes())
    for mir in validated_mirs:
        microcosm_targets = set(microcosm_predicted_fly_mir_targets(mir))
        fn = os.path.join(dn, '{}_overlap_of_nmj_rnai_genes_and_microcosm_fly_targets.csv'.format(mir))
        write_set_overlap_file(nmj_rnai, microcosm_targets, 'nmj_rnai_genes', 'microcosm', fn)
        targetscan_targets = set(targetscan_predicted_fly_mir_targets(mir))
        fn = os.path.join(dn, '{}_overlap_of_nmj_rnai_genes_and_targetscan_fly_targets.csv'.format(mir))
        write_set_overlap_file(nmj_rnai, targetscan_targets, 'nmj_rnai_genes', 'targetscan', fn)
        merg_fg = merged_fly_targets(mir)
        fn = os.path.join(dn, '{}_overlap_of_nmj_rnai_genes_and_merged_fly_targets.csv'.format(mir))
        write_set_overlap_file(nmj_rnai, merg_fg, 'nmj_rnai_genes', 'merged', fn)


def write_overlap_of_validated_mir_targets_and_conserved_targets_of_human_mir_homologs():
    dn = os.path.join(results_dir(), 'phase2', 'overlap_of_validated_mir_targets_and_conserved_targets_of_human_mir_homologs')
    # filename template
    fnt = 'overlap_of_{}_{}_targets_and_conserved_human_mir_{}_homolog_targets.csv'
    mfnt = 'overlap_of_{}_merged_targets_and_conserved_human_mir_homolog_targets.csv'
    for mir in validated_mirs:
        # merged fly targets
        merg_fg = merged_conserved_fly_targets(mir)
        # merged fly ortholog targets
        merg_fog = merged_fly_orthologs_of_targets_of_human_mirs_homologous_to_fly_mir(mir)
        fn = os.path.join(dn, mfnt.format(mir))
        write_set_overlap_file(merg_fg, merg_fog, 'fly_targets',
                               'human_target_orthologs', fn)
        for targets_db, ibanez_method in itertools.product(targets_dbs,
                                                           ibanez_methods):
            fn = os.path.join(dn, fnt.format(mir, targets_db, ibanez_method))
            print targets_db, ibanez_method, mir, fn
            fly_targets = short_conserved_fly_targets(mir, targets_db)
            fly_ortholog_targets = short_fly_orthologs_of_targets_of_human_mirs_homologous_to_fly_mir(
                mir, ibanez_method, targets_db)
            write_set_overlap_file(fly_targets, fly_ortholog_targets,
                                    'fly_targets', 'human_target_orthologs',
                                    fn)


# PHASE III
# Rank the 27 functionally validated miRs by percentage of predicted targets that are also in the NMJ RNAi genes.
# Rank the 27 functionally validated miRs by percentage of predicted targets that are also conserved synaptic genes.
# Also rank them by hypergeometric distribution.  
#   - sets: predicted targets, conserved synaptic genes, background: all conserved genes (or all flybase genes?)
#   - sets: predicted targets, NMJ RNAi genes. background: all flybase genes
#   - examine the sizes of the sets and the overlaps, so we can confirm that the statistics look reasonable.

def hypergeometric_pvalue(background_genes, target_genes, goal_genes):
    '''
    Compute the p-value that the sample target_genes, sampled from
    background_genes, contains as many objects from goal_genes as it does, when
    sampling at random.
    background_genes must contain every gene in target_genes and goal_genes.
    '''
    # The cdf function returns the probability that X <= x
    # The p-value is the probability of X >= x.
    # Since the hypergeometric distribution is discrete, P(X >= x) = 1 - P(X <= (x-1)).
    assert not target_genes - background_genes
    if goal_genes - background_genes:
        print 'len(background_genes)'
        print len(background_genes)
        print 'goal_genes'
        print goal_genes
        print 'goal_genes - background_genes'
        print goal_genes - background_genes
        print 'len(goal_genes)'
        print len(goal_genes)
        print 'len(goal_genes - background_genes)'
        print len(goal_genes - background_genes)
    assert not goal_genes - background_genes

    # total number of balls
    M = len(background_genes)
    # number of black balls
    n = len(goal_genes)
    # number of balls sampled
    N = len(target_genes)
    # number of black balls in sample
    x = len(target_genes & goal_genes)
    pvalue = 1.0 - scipy.stats.hypergeom.cdf(x - 1, M, n, N, loc=0)
    return pvalue


def write_ranked_147_mir_targets():
    print 'write_ranked_147_mir_targets'
    # write_ranked_mir_targets(the_147_mirs[:2], '147', print_validated=True)
    write_ranked_mir_targets(the_147_mirs, '147', print_validated=True)


def write_ranked_validated_mir_targets():
    print 'write_ranked_validated_mir_targets'
    # # If you want to shorten validated_mirs for testing
    # global validated_mirs
    # validated_mirs = validated_mirs[:2]
    write_ranked_mir_targets(validated_mirs, 'validated')


def write_ranked_mir_targets(mirs, mirs_name, print_validated=False):
    '''
    mirs: a list of current mirbase ids
    mirs_name: a string used to construct file names and table column names.
    print_validated: If True, print a column saying whether or not each miR
    is part of the validated set of mirs.
    '''
    print 'write_ranked_mir_targets'
    dn = makedirs(os.path.join(results_dir(), 'phase3',
                               'ranked_{name}_mirs'.format(name=mirs_name)))
    fnt = 'ranked_{name}_mirs_using_{targets}_and_{goal}_genes.csv'

    # background gene sets:  Used as the background for hypergeometric test
    conserved_genes = set(fly_conserved_genes())
    all_genes = set(fly_genes()) # all flybase Drosophila melanogaster genes in the datastore.

    # "other" gene sets, for lack of a better word.  The more of these genes in
    # the targets, the higher the ranking of the miR.
    # Note: a few NMJ RNAi FBgn ids are not in 'all_genes', perhaps because
    # they come from different versions of the flybase database.
    # Specifically 10 of the 468 RNAi genes are not in `all_genes`
    # Also 23 of the 201 RNAi gain of function genes are not in all_genes.
    # Most of these are from non-Dmel species, some were created after we
    # downloaded our data (UniProt 2013_04), and the remaining ones are 
    # inexplicably missing.
    rnai = set(get_nmj_rnai_genes())
    print 'rnai genes not in all genes:'
    print rnai - all_genes
    print len(rnai - all_genes)
    rnai = rnai & all_genes

    rnai_gof = set(get_nmj_rnai_gain_of_function_genes())
    print 'len(rnai_gof):', len(rnai_gof)
    print 'rnai gof genes not in all genes:'
    print rnai_gof - all_genes
    print len(rnai_gof - all_genes)
    rnai_gof = rnai_gof & all_genes

    conserved_synapse = set(fly_conserved_synapse_genes())

    other_pairs = [(rnai, 'nmj_rnai'), 
                   (conserved_synapse, 'conserved_synapse'),
                   (rnai_gof, 'nmj_rnai_gof')]
    for other_genes, other_name in other_pairs:
        # lookup tables of genes targeted by each mir for each method (microcosm, targetscan, both)
        microcosm_targets = {mir: set(short_fly_targets(mir, MICROCOSM)) for mir in mirs}
        targetscan_targets = {mir: set(short_fly_targets(mir, TARGETSCAN)) for mir in mirs}
        merged_targets = {mir: microcosm_targets[mir] | targetscan_targets[mir] for mir in mirs}
        for mir_targets, targets_name in ((microcosm_targets, MICROCOSM),
                                          (targetscan_targets, TARGETSCAN),
                                          (merged_targets, 'merged')):
            fn = os.path.join(dn, fnt.format(name=mirs_name, targets=targets_name, goal=other_name))
            print 'writing', fn
            write_ranked_mir_targets_sub(fn, mirs, mirs_name, mir_targets,
                                         conserved_genes, all_genes,
                                         other_genes, other_name,
                                         print_validated=print_validated)


def test_merged_targets_equals_microcosm_union_targetscan():
    '''
    The merged targets should equal the union of the microcosm targets and
    targetscan targets for a mir.  Test that this is true using the 147
    mirs.
    An exception is raised if it is not true for a mir.
    '''
    for mir in the_147_mirs:
        print 'testing', mir
        merged = merged_fly_targets(mir)
        microcosm = short_fly_targets(mir, MICROCOSM)
        targetscan = short_fly_targets(mir, TARGETSCAN)
        if set(merged) != set(microcosm) | set(targetscan):
            print 'merged != microcosm + targetscan'
            raise Exception()


def write_ranked_mir_targets_sub(filename, mirs, mirs_name, mir_targets,
                                 conserved_genes, all_genes,
                                 other_genes, other_name,
                                 print_validated=True):
    headers = ['{name}_group_mir'.format(name=mirs_name)]
    if print_validated:
        headers += ['validated']
    headers += ['num_target_genes', 'num_{other}_genes', 'num_in_intersection',
                'percent_target_genes_in_{other}_genes',
                'percent_{other}_genes_in_target_genes',
                'num_conserved_target_genes', 'num_conserved_{other}_genes',
                'num_in_conserved_intersection', 'all_hypergeometric_pvalue',
                'conserved_hypergeometric_pvalue',
                'num_all_genes', 'num_conserved_genes']
    template = ','.join('{}' for h in headers) + '\n'
    others = set(other_genes)
    num_others = len(others)
    conserved_others = others & conserved_genes
    num_conserved_others = len(conserved_others)
    with open(filename, 'w') as fh:
        fh.write(template.format(*headers).format(other=other_name))
        for mir in mirs:
            targets = set(mir_targets[mir])
            num_targets = len(targets)
            conserved_targets = targets & conserved_genes
            num_conserved_targets = len(conserved_targets)
            num_intersection = len(targets & others)
            num_conserved_intersection = len(conserved_targets & conserved_others)
            if num_targets > 0:
                percent_targets_in_others = num_intersection / float(num_targets)
            else:
                percent_targets_in_others = 'NaN'
            percent_others_in_targets = num_intersection / float(num_others)
            all_genes_hypergeom = hypergeometric_pvalue(
                background_genes=all_genes, target_genes=targets,
                goal_genes=others)
            conserved_genes_hypergeom = hypergeometric_pvalue(
                background_genes=conserved_genes, target_genes=conserved_targets,
                goal_genes=conserved_others)
            fields = [mir]
            if print_validated:
                fields += [mir in validated_mirs]
            fields += [num_targets, num_others,num_intersection,
                       percent_targets_in_others, percent_others_in_targets,
                       num_conserved_targets, num_conserved_others,
                       num_conserved_intersection, all_genes_hypergeom,
                       conserved_genes_hypergeom,
                       len(all_genes), len(conserved_genes)]
            fh.write(template.format(*fields))

# PHASE IV

def write_conserved_mir_targets_tables():
    '''
    For the validated, conserved and ranked mirs, write a table containing, 
    for each mir, the intersection of the predicted targets of the mir (using
    targetscan, microcosm, or both) with the fly orthologs of synaptomedb genes
    (i.e. conserved synaptic genes).
    '''
    mirs = validated_conserved_and_ranked_mirs
    conserved_synapse = set(fly_conserved_synapse_genes())
    human_synaptic_orthologs = fly_conserved_synapse_genes_human_orthologs_map()
    targets_dbs = [MERGED, MICROCOSM, TARGETSCAN]

    # map (mir, target_db) to the mir's targets in that targeting database
    to_targets = dict(
        ((mir, targets_db), fly_targets_set(mir, targets_db)) for mir in
        mirs for targets_db in targets_dbs)

    # filter to just conserved synaptic targets
    to_conserved_synaptic_targets = dict(
        (key, targets & conserved_synapse) for key, targets in to_targets.items())

    # map set of fly genes to set of fly-gene, human-gene pairs, for fly genes
    # with orthologous human genes or (fly-gene, None) if the fly gene has no
    # ortholog.
    def fly_to_pairs(fly_gene, orthologs):
        return set((fly_gene, ortholog) for ortholog in orthologs.get(fly_gene, [None]))

    # Fun fact: this can be done on one sickeningly complicated dict comprehension at the expense of any readability
    to_targets_and_ologs = {}
    for key, targets in to_conserved_synaptic_targets.items():
        # for each fly target, get the set of fly target, human ortholog pairs.
        # There will usually only be one ortholog, so this is usually a set
        # with one pair.  However there can be no ortholog or multiple
        # orthologs.  No ortholog results in a pair of (fly_gene, None).
        # Multiple orthologs result in multiple pairs.
        lists_of_pairs = [fly_to_pairs(fly_gene, human_synaptic_orthologs) for fly_gene in targets]
        # flatten the lists of pairs into a single set of pairs
        fly_human_pairs = set(itertools.chain.from_iterable(lists_of_pairs))
        to_targets_and_ologs[key] = fly_human_pairs

    # Transform to an output table with 2 rows per mir, one with the fly gene
    # targets and the other with the human orthologs.  Maintain the alignment
    # of the human-fly gene pairs, so we can tell what is orthologous.
    targets_db_tables = collections.defaultdict(list)
    for (mir, targets_db), pairs in to_targets_and_ologs.items():
        if pairs:
            fly_genes, human_genes = zip(*pairs)
        else:
            # Some targets database (e.g. microcosm) have no targets for some mirs.
            fly_genes = []
            human_genes = []
        fly_row = [mir] + list(fly_genes)
        human_row = [mir] + list(human_genes)
        targets_db_tables[targets_db].append(fly_row)
        targets_db_tables[targets_db].append(human_row)

    # Transpose the tables, so that each column starts with a mir and has genes
    # below it.  Alternating columns will have fly or human genes. And the
    # fly genes and human genes (for a given pair of columsn) that are adjacent
    # are orthologs.  Note that rows in the input columns can be different 
    # lengths for different mirs.
    targets_db_column_tables = {targets_db: itertools.izip_longest(*rows, fillvalue=None) for targets_db, rows in targets_db_tables.items()}

    dn = makedirs(os.path.join(results_dir(), 'mirs_and_target_genes'))
    fnt = 'mir_and_conserved_synaptic_fly_and_human_{targets}_targets.csv'
    for targets_db, rows in targets_db_column_tables.items():
        fn = os.path.join(dn, fnt.format(targets=targets_db))
        write_csv_file(rows=rows, filename=fn)



# May 9 Davie Grant: let 7, miR-13a, miR-13b miR-92a, miR-92b, miR-34, miR-137

def write_may_grant_mir_targets_tables():
    '''
    For each of the 7 mirs that Davie is interested in using in the 2014 May
    grant proposal, find the predicted targets (using microcosm, targetscan
    or both) and the overlap with conserved synaptic genes (synaptomedb genes
    with fly orthologs) and with NMJ RNAi GoF/LoF genes.
    '''
    pass
    mirs = may_mirs
    print mirs
    conserved_synapse = set(fly_conserved_synapse_genes())
    rnai = set(get_nmj_rnai_genes())
    rnai_gof = set(get_nmj_rnai_gain_of_function_genes())
    print len(conserved_synapse), len(rnai), len(rnai_gof)
    filter_to_genes = {
        'conserved_synaptic': conserved_synapse,
        'NMJ_RNAi': rnai,
        'NMJ_RNAi_GoF': rnai_gof
    }
    filters = sorted(filter_to_genes.keys())
    targets_dbs = [MERGED, MICROCOSM, TARGETSCAN]

    # for each targets database, map each mir to its gene targets
    # map each targets database to a map of mir targets.
    to_targets = dict(
        (targets_db, 
         dict((mir, fly_targets_set(mir, targets_db)) for mir in mirs))
        for targets_db in targets_dbs)

    # for each combination of filter and targets_db, create a map from mir to filtered mir targets
    to_filtered_targets = {}
    for filter_, targets_db in itertools.product(filters, targets_dbs):
        filtered_targets = {}
        for mir in mirs:
            filtered_targets[mir] = to_targets[targets_db][mir] & filter_to_genes[filter_]
        to_filtered_targets[(filter_, targets_db)] = filtered_targets

    # Create a table, one row for each mir and its filtered target genes
    # for each combination of filter and targets_db
    to_tables = {}
    for key in to_filtered_targets:
        table = [[mir] + sorted(to_filtered_targets[key][mir]) for mir in mirs]
        to_tables[key] = table

    # Transpose the tables, so each column starts with a mir and has genes.
    # Note that columns can be different lengths, since different mirs will
    # have different numbers of (filtered) targets
    to_column_tables = {key: itertools.izip_longest(*rows, fillvalue=None) for key, rows in to_tables.items()}

    dn = makedirs(os.path.join(results_dir(), 'mirs_and_target_genes'))
    # filename template
    fnt = 'mir_and_{filter}_fly_{targets}_targets.csv'
    for (filter_, targets_db), rows in to_column_tables.items():
        fn = os.path.join(dn, fnt.format(filter=filter_, targets=targets_db))
        write_csv_file(rows=rows, filename=fn)




# RANDOM AND ONE-OFF TASKS

def write_microcosm_mir_to_gene_table():
    fn = os.path.join(results_dir(), 'microcosm_mirbase_id_to_flybase_gene_id.csv')
    headers = ['mirbase_id', 'flybase_gene_id']
    table = microcosm_fly_targets_table()
    write_csv_file(rows=table, filename=fn, headers=headers)
    # for mirbase_id, flybase_gene in table:
        # print mirbase_id, flybase_gene
    # # for mir, current_mir, gene in table:
        # # print mir, current_mir, gene
    # print len(table)


###############
# MAIN WORKFLOW


def workflow():

    # Clear out a graph loaded into the wrong graph URI
    # Slow for big graphs.
    # virt.clear_graph('http://purl.roundup.hms.harvard.edu/graph/uniprot_taxonomy')

    if False: # Done tasks
        pass

        # Download initial data files
        download_microcosm_targets()
        download_mirbase_aliases()
        download_targetscan()
        # Not Implemented b/c not available online or not available anymore for
        # the version we want (uniprot 2013_04) or just not done yet (Roundup,
        # Flybase)
        # download_affymetrix_fly_annotations_file()
        # download_ibanez_fly_human_homologs()
        # download_mcneill_screen()
        # download_roundup_orthologs()
        # download_synaptome_v1()
        # download_flybase_transcript_reporting_xml()

        # Convert databases to RDF
        write_affymetrix_fly_annotations_rdf()
        write_flybase_rdf()
        write_ibanez_five_prime_mir_homologs_rdf()
        write_ibanez_fullseq_mir_homologs_rdf()
        write_mcneill_cns_screen_rdf()
        write_mcneill_muscle_screen_rdf()
        write_microcosm_rdf()
        write_mirbase_rdf()
        write_synaptomedb_rdf()
        write_targetscan_rdf()

        # Construct Roundup Orthologs from sparql.roundup.hms.harvard.edu
        download_roundup_from_sparql()
        # Or get them from the roundup text files from a raw data download
        write_roundup_orthologs_rdf()

        # Get UniProt ID Mapping RDF from beta.sparql.uniprot.org.
        set_uniprot_version()
        download_uniprot_fly_rdf()
        download_uniprot_homo_rdf()
        download_uniprot_flybase_annotation_rdf()
        download_uniprot_flybase_gene_rdf()
        download_uniprot_flybase_transcript_rdf()
        download_uniprot_homo_refseq_rdf()
        download_uniprot_homo_ensembl_transcript_rdf()
        download_uniprot_homo_ensembl_gene_rdf()
        # Or get UniProt ID Mapping RDF from idmapping.dat (UniProt ftp file.)
        # write_uniprot_rdf()

        load_rdf_database()

        # Simplify and speed queries by constructing edges between microRNAs
        # and genes
        construct_microcosm_human_mirna_id_to_gene_edges()
        construct_microcosm_fly_mirna_id_to_gene_edges()
        construct_targetscan_human_mirna_id_to_gene_edges()
        construct_targetscan_fly_mirna_id_to_gene_edges()
        construct_roundup_human_gene_to_fly_gene_edges()
        construct_ibanez_fullseq_edges()
        construct_ibanez_five_prime_edges()
        load_constructed_edges()

        # PHASE ZERO
        write_conserved_synaptic_genes()
        write_affymetrix_to_flybase_table()
        write_human_to_fly_conserved_synapse_genes_table()

        # This step requires that the mirbase data be loaded, so ids can be updated.
        # The step is part of a manual process that resulted in the_147_mirs being
        # assigned a literal array.  To fully automate, one could save the
        # array to a file when creating it and then read it from the file when
        # using it.
        update_original_147_mirs() # A step in the manual process to update ids

        # PHASE I
        write_phase_1()
        # write_targetscan_human_mirs_targeting_conserved_synapse_genes()
        # write_microcosm_human_mirs_targeting_conserved_synapse_genes()
        # write_overlap_between_microcosm_human_mirs_and_targetscan_human_mirs()
        # write_targetscan_fly_mirs_targeting_conserved_synapse_genes()
        # write_microcosm_fly_mirs_targeting_conserved_synapse_genes()
        # write_overlap_between_validated_fly_mirs_and_microcosm_predicted_fly_mirs()
        # write_overlap_between_validated_fly_mirs_and_targetscan_predicted_fly_mirs()
        # write_overlap_between_microcosm_fly_mirs_and_targetscan_fly_mirs()

        # PHASE II
        write_overlap_between_screened_and_microcosm_and_targetscan_fly_mir_targets()
        write_overlap_of_screened_fly_mir_targets_and_nmj_rnai_genes()
        write_overlap_of_validated_fly_mir_targets_and_nmj_rnai_genes()
        write_overlap_of_validated_mir_targets_and_conserved_targets_of_human_mir_homologs()

        # PHASE III
        write_ranked_validated_mir_targets()
        write_ranked_147_mir_targets()
        write_microcosm_mir_to_gene_table()

        # PHASE IV
        write_conserved_mir_targets_tables()

        # 2014 MAY GRANT
        write_may_grant_mir_targets_tables()

    else:
        start_stardog()
        # load_rdf_database()
        # load_constructed_edges()

        # write_ranked_validated_mir_targets()
        # write_ranked_147_mir_targets()
        # write_microcosm_mir_to_gene_table()
        write_may_grant_mir_targets_tables()

        # debug_pairs(short_human_to_fly_conserved_genes_table(), 'conserved genes')
        # debug_pairs(human_to_fly_conserved_synapse_genes_table(), 'conserved synapse genes')
        # debug_pairs(short_human_to_fly_conserved_synapse_genes_table(), 'short conserved synapse genes')

        # fg = fly_genes()
        # print len(fg)
        # print len(set(fg))
        # print len(fly_conserved_genes())

        # test_merged_targets_equals_microcosm_union_targetscan()

        # vs = set(validated_mirs)
        # t147s = set(the_147_mirs)
        # print 'Validated mirs should all be in the 147 mirs, I think:'
        # print 'Number of validated mirs not in the 147 mirs:', len(vs - t147s)
        # print 'Number of 147 mirs in validated mirs:', len(t147s & vs)
        # print 'Size of the 147 mirs set:', len(t147s)
        # print 'Size of the validated mirs set:', len(vs)
        pass


def debug_pairs(tbl, desc=None):
    '''
    Some functions return a list of pairs of human and fly genes.  This
    prints out how many pairs there are, how many unique pairs, and
    how many unique human and fly genes.
    '''
    print desc
    # print tbl
    print 'all pairs', len(tbl)
    print 'unique pairs', len(set(tbl))
    h, f = zip(*tbl)
    h = set(h)
    f = set(f)
    print 'unique human ids', len(h)
    print 'unique fly ids', len(f)


def main():
    workflow()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception('')
        raise


