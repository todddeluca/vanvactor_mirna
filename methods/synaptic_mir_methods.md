ROUGH DRAFT

# Synaptic miR Analysis Methods


## What data sets were downloaded? Where did they come from? What version?

In order to analyze predicted targets of miRs in relation to conserved synaptic genes, data was downloaded from the sources described below, transformed into RDF triples where necessary, and loaded into a StarDog 1.2 RDF database for querying and further processing.

### Microcosm

Predicted microRNA targets were downloaded from the microcosm database (version 5) for human and fly:

- ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt.mus_musculus.zip
- ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt.homo_sapiens.zip

### TargetScan

Predicted targets were downloaded from the TargetScan database for human (version 61) and fly (version 12):

- http://www.targetscan.org/vert_61/vert_61_data_download/miR_Family_Info.txt.zip
- http://www.targetscan.org/vert_61/vert_61_data_download/Predicted_Targets_Info.txt.zip
- http://www.targetscan.org/fly_12/fly_12_data_download/flybase2symbol_fly12_targetscan.txt
- http://www.targetscan.org/fly_12/fly_12_data_download/Conserved_Family_Conserved_Targets_Info.txt.zip
- http://www.targetscan.org/fly_12/fly_12_data_download/miR_Family_Info.txt.zip

In human, data was downloaded to map from mir "families" to NCBI RefSeq Ids and from mir "family" to miRBase Accession.  In fly, data was downloaded to map from mir "families" to the fly gene symbols they are predicted to target, from mir "family" to miRBase Accession, and from fly gene symbol to flybase gene id.

### miRBase

Aliases between miRBase Accessions and miRBase Ids were downloaded from miRBase (version 19):

- ftp://mirbase.org/pub/mirbase/19/aliases.txt.gz

### Roundup

Orthologs between human and fly proteins were downloaded from the Roundup orthology database (version 4).  This was done by going to http://roundup.hms.harvard.edu/download and downloading Homo sapiens (9606) and Drosophila melanogaster (7227) orthologs with divergence and evalue thresholds of 0.8 and 1e-5 respectively.

### Affymetrix

From Affymetrix, annotations were downloaded for release 33 of the Drosophila Genome 2.0 array.  This was done by going to http://www.affymetrix.com/support/technical/byproduct.affx?product=fly-20 and downloading Drosophila_2 Annotations, CSV format, Release 33 (8.4 MB, 10/30/12).  From this we parsed the mapping between affymetrix probeset ids and flybase gene ids.

### Ibanez

Homologous miRs between fly and human were taken from the 2008 Ibanez-Ventoso paper.

### McNeill

Experimental results for Affymetrix probeset ids that were differentially expressed when one of 7 miRs was perturbed in fly were obtained from McNeill (unpublished data).

26 fly miRs were functionally validated by McNeill et al. by looking for morphological changes in the muscle phenotype.

### SynaptomeDB

Human genes that are considered synaptic were downloaded from SynaptomeDB (version 1.06) at http://psychiatry.igm.jhmi.edu/SynaptomeDB/index.php.
[Ensembl Gene Ids]

### FlyBase

The relationships between FlyBase Annotation Ids and FlyBase Transcript Ids were downloaded from FlyBase version FB2013_02 at ftp://ftp.flybase.net/releases/FB2013_02/reporting-xml/FBtr.xml.gz.
[Used to map from]

### UniProt

RDF triples were downloaded from UniProt (version 2013_07) via the SPARQL endpoint http://beta.sparql.uniprot.org/.  The following data was downloaded:

- Human and Fly UniProt accessions
- Mappings from UniProt accessions to FlyBase Accessions, FlyBase Gene Ids, FlyBase Transcript Ids, Human NCBI RefSeq Ids, Human Ensembl Transcript Ids, and Human Ensembl Gene Ids.

### NMJ RNAi Genes

[TODO] These came from Elizabeth.  I'm not sure where she got them.  There is a list of Gain of Function and Loss of Function genes, as well as a list of just Gain of Function genes.

## Relationships Between miRs and the Genes

SPARQL Construct queries were written to map from miRs to genes and load the resulting RDF triples into the database for further querying.  The specific mappings are detailed below.

### Microcosm.

The relationship between human genes and the miRs predicted to target them was established by mapping from a current miRBase Id to a miRBase Accession to a miRBase Id (possibly the current one) to an Ensembl Transcript Id to a UniProt Id to an Ensembl Gene Id.  This was done using data from miRBase, Microcosm, and UniProt.

The relationship between fly genes and miRs was established by mapping from a current miRBase Id to a miRBase Accession to a miRBase Id (possibly the current one) to a FlyBase Annotation Id to a FlyBase Transcript Id to a UniProt Id to a FlyBase Gene Id, using data from miRBase, Microcosm, FlyBase, and UniProt.

### TargetScan

The relationship between human genes and the miRs predicted to target them was established by mapping from a current miRBase Id to a miRBase Accession to a TargetScan miR Family to a RefSeq Transcript Id to a UniProt Id to an Ensembl Gene Id, using data from miRBase, TargetScan, and UniProt.

The relationship between fly genes and miRs was established by mapping from a current miRBase Id to a miRBase Accession to a TargetScan miR Family to a FlyBase Gene Id, using data from miRBase and TargetScan.


## Human-Fly Orthologs

Orthologous relationships between human and fly genes were established by mapping from human Ensembl Gene Ids to human UniProt Ids to fly UniProt Ids to FlyBase Gene Ids.  A SPARQL query was used to map the identifiers, creating new RDF triples to load into the database.

## Phase I

### Fly miRs that Target Conserved Synaptic Genes

Conserved synaptic genes are fly genes that are orthologous to human genes that are synaptic according to SynaptomeDB.

A SPARQL query was used to determine what fly miRBase Ids are predicted to target conserved synaptic genes.  This set of miRBase Ids was then compared to the set of 26 miRs functionally validated by McNiell.

## Phase II

For each of the seven miRs screened by McNeill, compare the differentially expressed genes in Muscle and Central Nervous System tissues from the screen with the predicted target genes (from Microcosm and TargetScan) as well as the NMJ RNAi GoF/LoF genes.

For each of the 26 functionally validated miRs from McNeill, compare the predicted target genes to the RNAi NMJ genes.  Also compare the predicted targets to the conserved fly genes of the human predicted targets of the human miR homologs of the fly miRs.

## Phase III

### Ranking 144 miRs

144 drosophila miRs were ranked according to:

- the percentage of their gene targets overlap with NMJ RNAi genes
- the percentage of their gene targets that overlap with the conserved synaptic fly genes
- a hypergeometric test p-value comparing the miR target genes to the NMJ RNAi genes, using all FlyBase Gene Ids as a background.
- a hypergeometric test p-value comparing the target genes of each miR to the conserved synaptic fly genes, using all conserved fly genes as a background.  

The hypergeometric test p-value was calculated using the cumulative density function of the hypergeometric distribution from the `stats` module of SciPy, an open source library of scientific tools written in Python.

### Rank the 27 functionally validated miRs by percentage of predicted targets that are also in the NMJ RNAi genes.
### Rank the 27 functionally validated miRs by percentage of predicted targets that are also conserved synaptic genes.
### Also rank them by hypergeometric distribution.  
- sets: predicted targets, conserved synaptic genes, background: all conserved genes (or all flybase genes?)
- sets: predicted targets, NMJ RNAi genes. background: all flybase genes
- examine the sizes of the sets and the overlaps, so we can confirm that the statistics look reasonable.


How were relationships between miRs and genes established in fly and in human?
What were the original edges and transitive edges that were generated?
How were the relationships between fly and human genes established, including original and transitive edges?
How were miRs ranked?  By percentage and by hypergeometric distribution?
How were the genes in the overlaps established?
