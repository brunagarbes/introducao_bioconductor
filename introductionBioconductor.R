#What is bioconductor?

# Load the BiocInstaller package
library(BiocInstaller)

# Explicit syntax to check the Bioconductor version
BiocInstaller::biocVersion()

# When BiocInstaller is loaded use biocVersion alone
biocVersion()

# Load the BSgenome package
library(BSgenome)

# Check the version of the BSgenome package
packageVersion("BSgenome")

# Investigate about the a_genome using show()
show(a_genome)

# Investigate some other accesors

organism(a_genome)
provider(a_genome)
seqinfo(a_genome)

# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the head of seqnames and tail of seqlengths for yeastGenome
head(seqnames(yeastGenome))
tail(seqlengths(yeastGenome))

# Select chromosome M, alias chrM
getSeq(yeastGenome, "chrM")

# Count characters of the chrM sequence
nchar(yeastGenome$chrM)
# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the first 30 bases of each chromosome
getSeq(yeastGenome,start = 1, end = 30)