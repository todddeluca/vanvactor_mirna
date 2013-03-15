
# Configuration
datadir = '/Users/td23/deploy/vanvactor_mirna/data'

# get the latest version of bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")

# Venn Diagrams in R

# Can not get Vennerable to install
# biocLite("vennerable")
# install.packages("Vennerable", repos="http://R-Forge.R-project.org")
# install.packages(c("graph", "RBGL"), dependencies=TRUE)
# library(Vennerable)

biocLite('VennDiagram')
library(VennDiagram)

# Read in gene targets of 2b and 13b and get their intersection and set differences
targets.2b = read.table(paste(datadir, 'microcosm/v5/microcosm-v5-dme-miR-2b-targets-20130205.flybase.genes.txt', sep='/'), stringsAsFactors=F)[,1]
targets.2b
length(targets.2b)
targets.13b = read.table(paste(datadir, 'microcosm/v5/microcosm-v5-dme-miR-13b-targets-20130205.flybase.genes.txt', sep='/'), stringsAsFactors=F)[,1]
targets.13b
length(targets.13b)
both = intersect(targets.2b, targets.13b)
both
length(both)
just2b = setdiff(targets.2b, targets.13b)
just2b
length(just2b)
just13b = setdiff(targets.13b, targets.2b)
just13b
length(just13b)

venn.plot = draw.pairwise.venn(length(targets.2b), length(targets.13b), length(both), c('miR-13b target genes', 'miR-2b target genes'))
# Writing to file
tiff(filename = "~/tmp/microcosm-v5-dme-miR-2b-and-13b-flybase-gene-targets-20130205.venn.diagram.tiff", compression = "lzw");
grid.draw(venn.plot);
dev.off();


