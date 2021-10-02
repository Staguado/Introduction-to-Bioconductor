####################################################################################################################

# By: Santiago Taguado Menza
# Introduction to Bioconductor: 
# April 21st, 2021
# Harvard University

####################################################################################################################

# Introduction to Bioconductor Section 2
# Bioconductor Basics with GRanges and Biostrings
# This week, we introduce core Bioconductor structures for representing genes and genetic sequences: 
# GRanges and Biostrings.

# Learn core Bioconductor data structures while exploring a ChIP-seq dataset.
# Use GRanges to organize genomic features such as genes and exons in R.
# Use Biostrings to organize genetic sequences in R.
# Discover how to efficiently manipulate GRanges and Biostrings objects using Bioconductor functions.

BiocManager::install(c("Homo.sapiens",
                       "GenomicFeatures",
                       "genomicsclass/ERBS",
                       "genomicsclass/ph525x"))
BiocManager::install("genomicsclass/ERBS")
library(ERBS)

####################################################################################################################

# Assessment : Setup and GRanges
data(HepG2)

# Question 1
# What is the class of HepG2?
class(HepG2)

# Question 2
# How many regions are represented?
length(HepG2)

####################################################################################################################

# Assessment: Introduction to Genomic Ranges

# In the video we used the values() method to extract meta-data on the regions. An alternative, and actually preferred 
# approach going forward, is mcols().

# Question 1
# What is the median of the signalValue column for the HepG2 data?
median(HepG2$signalValue)
dat_HepG2 = mcols(HepG2)

# Question 2
max(HepG2$signalValue)
which.max(HepG2$signalValue)
HepG2[120,]

# Question 3
# How many regions are from chromosome 16?
library(GenomicRanges)
x = seqnames(HepG2)
as.character(x)
table(x)[1:24]

# Question 4: Statistics on peak width
# Make a histogram of the widths of the regions from all chromosomes (not just chr16). 
# Note it has a heavy right tail.
# What is the median width?
hist(width(HepG2))
median(width(HepG2))

####################################################################################################################

# Assessment: iRanges

browseVignettes("IRanges")

# This is just a subset of all the possible operations, and remember, the rest are documented in the help pages mentioned 
# in the video and in the book page. We will first do a simple review of these operations, so that you get a sense 
# of using them in your R console. Then we will have a few questions which require more thought.

# Load the IRanges package. Define an integer range starting at 101 and ending at 200. 
# If we use the operation *2, this will zoom in, giving us a range with half the width. 

# Question 1
# What is the starting point of the resulting range?

(ir = IRanges(101, end = 200))

# Question 2
# Define an integer range starting at 101 and ending at 200. 
# If we use the operation narrow(x, start=20), what is the new starting point of the range?
  
narrow(ir,start=20)

# Question 3
# If we use the operation +25, what is the width of the resulting range?
ir + 25 

# Question 4
# Define an IRanges with starts at 1,11,21 and ends at 3,15,27. width() gives the widths for each range.
(ir_2 = IRanges(start = c(1,11,21),end = c(3,15,27)))
sum(width(ir_2))

# Question 5 
# Define an IRanges object, x, with the following set of ranges:
# Starts at 101,106,201,211,221,301,306,311,351,361,401,411,501
# Ends at 150,160,210,270,225,310,310,330,390,380,415,470,510

# Plot these ranges using the plotRanges() function in the ph525x package. You can install this library, 
# if you have not done so already, with the command install_github("genomicsclass/ph525x").
# What is the total width from 101 to 510 which is not covered by ranges in x?

BiocManager::install("genomicsclass/ph525x")
library(ph525x)
ir_3 = IRanges(start = c(101,106,201,211,221,301,306,311,351,361,401,411,501),
               end = c(150,160,210,270,225,310,310,330,390,380,415,470,510))
plotRanges(ir_3)
sum(width(gaps(ir_3)))

# Question 6
#By "disjoint ranges", we mean the following: for two ranges [1,10] and [6,15], there are three disjoint ranges 
# contained within: [1,5], [6,10], and [11,15].
# How many disjoint ranges are contained within the ranges in x from the previous question?
  
length(disjoin(ir_3))

# Question 7
# An intra-range function we didn't show in the video is resize().
# Set up a grid of 2 stacked plots:
par(mfrow=c(2,1))

# Now use plotRanges() to plot the x from last question, as well as resize(x,1). 
# You will have to set the xlim to make sure that the plots line up vertically. 
# For example, you can use plotRanges(x, xlim=c(0,600)).
# What is the best description for the operation resize(x,1)?  

(ir_4 = resize(ir_3,1))
plotRanges(ir_3,xlim = c(0,600))
plotRanges(ir_4, xlim = c(0,600))

# It gives you just the starting point of each range.
# From the man page: resize() resizes the ranges to the specified width where either the start, end, or center is 
# used as an anchor. The default is fix="start", so resize(x,1) gives you the starting integer of each range in x.

####################################################################################################################

# Assessment: GRanges

# Question 1
# Understanding strand orientation with resize
#In the first week, in the subsection "What We Measure and Why", we learned that DNA has two strands. 
# These two strands are often called plus, "+", and minus, "-".
# The GRanges object in the GenomicRanges package extends the concept of interval ranges in two major ways. 
# The ranges are now also identified by:
  
# 1. the chromosome we are referring to (in Bioconductor, this is called "seqnames")
# 2. the strand of the DNA we are referring to ("+" or "-"). No strand is labelled with a star, "*".

# Without these two pieces of information, a specification of a range of DNA would be ambiguous. 
# Let's make two ranges, with strand and chromosome information, and see how the range operations act based on strand.

x = GRanges("chr1", IRanges(c(1,101),c(50,150)), strand=c("+","-"))

# In the last assessment, we visualized IRanges with the plotRanges() function in the ph525x library. 
# We can get the internal IRanges from a GRanges object with the following code:

ranges(x)

# So let's define a new plotting function:
plotGRanges = function(x) plotRanges(ranges(x))

# Compare x and resize(x,1) using plotGRanges. The result of running resize(x,1) is two ranges of width 1 which start...
plotGRanges(resize(x,1))

# at the left-most point of the "+" strand ranges in x, and the right-most point of the "-" strand ranges in x

# Question 2
# Suppose we have two different sets of ranges, which overlap somewhat but not entirely. 
# This is the case for many genes, in which there are different versions of transcripts, also called isoforms. 
#The different transcripts consist of exons which end up in the final mRNA molecule, and a number of transcripts 
# can share exons or have exons which are overlapping but not identical ranges.
# We'll start with a toy example, and learn how to load real genes later:

x = GRanges("chr1", IRanges(c(101,201,401,501),c(150,250,450,550)), strand="+")
y = GRanges("chr1", IRanges(c(101,221,301,401,541),c(150,250,350,470,550)), strand="+")

# Plot these two sets of ranges using par(mfrow=c(2,1)) and two calls to plotGRanges.
# If we want to keep the information about which set the ranges belong to, we could combine the two GRanges into a 
# GRangesList:

par(mfrow = c(2,1))
plotGRanges(x)
plotGRanges(y)

u = GRangesList(x,y)
# Hint: use c(), disjoin() and %over%.

# Find the total width which is covered by ranges in both x and y.

disjoined = disjoin(c(x,y))
in.both = (disjoined %over% x & disjoined %over% y)
sum(width(disjoined[ in.both ]))

# Question 3 
# What is the total width which is in x or y but not in both?

disjoined = disjoin(c(x,y))
not.in.both = !(disjoined %over% x & disjoined %over% y)
sum(width(disjoined[ not.in.both ]))

# Question 4
# Define a new genomic range, z, which covers range(ranges(x)) but has the opposite strand.

z  = GRanges("chr1", IRanges(c(101),c(550)), strand="-")

# What is the number of ranges in x which overlap z according to the %over% command?
sum(x %over% z)

####################################################################################################################

# Assessment: Finding Overlaps

library(ERBS)
data(HepG2)
data(GM12878)

# Question 1
# Where does the 17th HepG2 region start?
HepG2[17]

# Question 2
# Use distanceToNearest() to find the closest region in GM12878 to the 17th region in HepG2. 
# What is the start site of this region?
  
distanceToNearest(HepG2[17],GM12878)
GM12878[945]

# Question 3: Measuring distance between closest regions
# What is the distance between the 17th region of HepG2 and its closest region in GM12878?
  
distance(HepG2[17],GM12878[945])

# Q4: Summarizing proximities of nearest regions in a pair of GRanges

dis = distanceToNearest(HepG2,GM12878)
dis_1 = mcols(dis)
mean(dis_1[1:303,] < 2000)

####################################################################################################################

# Assessment: Genes as GRanges

# Question 1
# What genome build was used here?
library(Homo.sapiens)
ghs = genes(Homo.sapiens)

# How many genes are represented in ghs?
seqnames = "chromosomes"
length(ghs)

# Question 2
# What is the chromosome with the most genes?
seqnames(ghs)
chr = seqnames(ghs)
as.character(chr)
max(table(chr)[1:24])

# Question 3
# Make a histogram of the widths of genes (use the width() on the GRanges object). 
# This width gives the number of basepairs from the start of the gene to the end, so including exons and introns.
# Which best describes the width of genes?
hist(width(ghs),breaks =1000)

# Skewed to the right

# Question 4
# What is the median of the width?
median(width(ghs))

####################################################################################################################

# Assessment: Finding and getting annotation for closest gene (Transcription start site)

# In this assessment (which deals with topics in several videos) we will find the closest genes to some of our 
# binding sites. We will use a consensus set of regions. In the video we did it like this:
  
library(ERBS)
data(HepG2)
data(GM12878)

res = findOverlaps(HepG2,GM12878)
erbs = HepG2[queryHits(res)]
erbs = granges(erbs)

# The following command is similar:
erbs2= intersect(HepG2,GM12878)

# Compare the results of these two commands.

# Which of the following is true of the consensus entities erbs and erbs2?
# Over 90% of these regions in these two objects are the same with the different regions being smaller in erbs2.

## first order them
(erbs3 = erbs[order(erbs),])
##confirm same chr
all( seqnames(erbs2)==seqnames(erbs3) )
mean( start(erbs2)==start(erbs3) & end(erbs2)==end(erbs3) )
##the intersection should be smaller
all( width(erbs2) <= width(erbs3) )

# Question 2
# What is the TSS (Transcription Start Site) of the gene with ID: 100113402?

library(Homo.sapiens)
ghs = genes(Homo.sapiens)
ghs["100113402"]

tss = start(ghs)
tssghs = resize(ghs,1)
tssghs["100113402"]

# Question 3 
# Now use the erbs regions defined in a previous question:
# What is the GENEID of the gene with TSS closest to the 4th region of erbs?
  
library(ERBS)
data(HepG2)
data(GM12878)
(res = findOverlaps(HepG2,GM12878))
(erbs = HepG2[queryHits(res)])
erbs = granges(erbs)
erbs[4]
nearest(erbs[4],tssghs)
tssghs[6316]

# Q4: Translating ENTREZ ID to symbol
# Use the select() function to determine the SYMBOL of this gene.
keys = as.character()
select(Homo.sapiens,)
keys = as.character(tssghs[1:23056]$GENEID)
(res = select(Homo.sapiens, keys = keys,
             columns = c("SYMBOL", "GENENAME"), keytype = "GENEID"))
res[6316,]

# Or

i = 6316
(gene = as.character(mcols(tssghs)$GENEID[i]))
select(Homo.sapiens,key=gene,column="SYMBOL",keytype="GENEID")

####################################################################################################################

# Assessment: Biostring

library(Biostrings)
eco <- DNAString("GGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGCTGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGA")

# Question 1
# How many bases are in eco?
  
length(eco)

# Question 2
# In DNA, the start codon is encoded by the trinucleotide ATG.
# How many potential start codons are in the eco sequence?

countPattern("ATG",eco)

# Question 3
# Find the locations of these ATG trinucleotides.
# What is the start location of the first ATG?

matchPattern("ATG",eco)

# Question 4
# Take the substring of eco that starts at the location found in the previous question. Translate this sequence 
# into amino acids.
# How long is the resulting AAString?

aa_full = translate(eco[17:181])
length(aa_full)

# Question 5
# In AAStrings, the * character represents a stop codon.
# What is the location of the stop codon in the AAString?

start(matchPattern("*", aa_full))

# Question 6
# Subset the AAString to end just before the stop codon; that is, do not include the stop codon or anything 
# following it in the AAString.

# How many amino acids are in the resulting AAString?
  
aa_full_1 = aa_full[1:52]
length(aa_full_1)

# What is the resulting amino acid sequence?
  
aa_full_1

# Qestion 7
# At pH 7, two amino acids are negatively charged (D, E) and three amino acids are positively charged (H, K, R). 
# Assume this peptide is at pH 7. The net charge of a peptide is the number of positively charged amino acids minus 
# the number of negatively charged amino acids.

# How many positively charged amino acids are in the sequence?
letterFrequency(aa_full_1,"H")
letterFrequency(aa_full_1,"K")
letterFrequency(aa_full_1,"R")
# How many negatively charged amino acids are in the sequence?
letterFrequency(aa_full_1,"D")
letterFrequency(aa_full_1,"E")
# What is the net charge of this peptide at pH 7?
letterFrequency(aa_full_1, "HKR") - letterFrequency(aa_full_1, "DE")

####################################################################################################################

# Sample Code on doing Motif Data Analysis

library(ERBS)
# Extracted from human liver cells
data(HepG2)

# load and inspect human reference genome
library(BSgenome.Hsapiens.UCSC.hg19)

# extract chromosome 17 sequence
c17 = Hsapiens$chr17
class(Hsapiens)
showMethods("getSeq")

# collection of DNA strings with ChIP-seq binding peaks
hepseq = getSeq(Hsapiens, HepG2)
length(HepG2)    # same number of sequences
width(HepG2)[1:5]    # widths match

# collection of shifted DNA strings with no relationship to binding sequences - essentially random
rhepseq = getSeq(Hsapiens, shift(HepG2, 2500))

# count occurrences of a motif in DNA sequences
mot = "TCAAGGTCA"
?vmatchPattern
vcountPattern(mot, hepseq)

# consider both forward matches and reverse complement matches 
sum(vcountPattern(mot, hepseq))    # forward pattern match
sum(vcountPattern(mot, reverseComplement(hepseq)))    # reverse pattern match

## compare motif occurrence in binding peak to random upstream sequences
# count of motifs in binding peaks
sum(vcountPattern(mot, hepseq)) +
  sum(vcountPattern(mot, reverseComplement(hepseq)))
# count of motifs in randomly selected regions of equal length
sum(vcountPattern(mot, rhepseq)) +
  sum(vcountPattern(mot, reverseComplement(rhepseq)))

# for real analysis, use MotifDb package, probabilistic binding packages like MEME and FIMO

####################################################################################################################

# Getting Sequence

library(ERBS)
library(GenomicRanges)
data(HepG2)
data(GM12878)

res = findOverlaps(HepG2,GM12878)
erbs = HepG2[queryHits(res)]
erbs = granges(erbs)
length(erbs)

# Question 1
# What genome build was used to create these regions?

genome(erbs)

# Question 2

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
HSapiens

# Use the getSeq() function to extract the sequence of each region in erbs. 
# Then compute the GC-content (the number of C's + the number of G's divided by the length of sequence) of each.

# What is the median GC-content?

x = getSeq(Hsapiens,erbs)
gc_cont = (vcountPattern("G",x) + vcountPattern("C",x))/width(x)
median(gc_cont)
hist(gc_cont)

y = getSeq(Hsapiens,shift(erbs,10000))
gc_cont_r = (vcountPattern("G",y) + vcountPattern("C",y))/width(y)
median(gc_cont_r)
hist(gc_cont_r)

# Why do I care for guanine and cytosine?

# A likely explanation for this is that we have higher GC-content in the promoter of genes 
# (upstream from the transcription start site). 
# And as we saw in the videos and earlier assessments, our binding sites tend to be close to the transcription 
# start sites.

200 - ( 135 + 1 + 0.35) 

####################################################################################################################

# Required Packages
library(GenomicRanges)
library(Biostrings)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ERBS)

# Verified Assessment: GRanges and Biostrings
# Question 1
# What genome build are these GRanges from?

g = genes(Homo.sapiens)
genome(g)[1]

# How many genes are in this build of the human genome?
length(g)

# How many unique seqlevels() are in the human genome?
length(unique(seqlevels(g)))

# Question 2
chr21_1 <- keepSeqlevels(g, "chr21", pruning.mode = "coarse")

# How many genes are on chromosome 21?
length(chr21_1)

# How many base pairs is the longest gene on chromosome 21?
max(width(chr21_1))

# What proportion of genes on chromosome 21 are on the + strand?
mean(strand(chr21_1) == "+")

# Question 3
# What is the median GC content of genes on chromosome 21?
z = getSeq(Hsapiens,chr21_1)
gc_cont = rowSums(letterFrequency(z,"GC"))/width(z)
median(gc_cont)

# Examine the sequence of the fifth gene on chr21. How many possible start codons (ATG) are in this sequence?
z[5]
vcountPattern("ATG",z[5])

# Suppose that the first ATG is the real start codon. What is the start location within this gene of the first ATG?
start(vmatchPattern("ATG",z[5]))

# Translate the sequence starting at the location identified in the previous question part. 
# Recall that * is the stop codon.
gene_5 = z[5]
translate(DNAStringSet(gene_5,start = 32))

# Question 5
# Let's use GRanges to investigate the SNCA gene. The SNCA gene encodes a protein called alpha-synuclein. 
# Abnormal accumulation of alpha-synuclein inside neurons is hypothesized to cause the clinical manifestations of 
# Parkinson's disease and other related neurodegenerative disorders known as synucleinopathies.

# The human SNCA gene is uniquely identified by the Gene ID, 6622. We will learn how to find the ID number of a desired 
# gene in a later section.
snca <- genes(Homo.sapiens, filter=list(GENEID="6622"))
width(snca)

# On what chromosome is the SNCA gene?
seqnames(snca)

# In what location is the TSS of SNCA?
resize(snca, 1)

# Get the sequence of SNCA. Suppose a protein binds to the DNA motif "ACTGTGAA". How many locations can this 
# protein bind to the SNCA gene?
l = getSeq(Hsapiens,snca)
motif = "ACTGTGAA"
vcountPattern(motif,l) + vcountPattern(motif,reverseComplement(l))

# In Questions 6 - 8 you will analyze the abundance of ESRRA binding sites in promoter regions relative to other 
# parts of the genome. Begin by defining the set of human genes and the set of binding sites in the GM12878 cell line.
data(GM12878)
class(g)

# Question 6
# Define the promoter region of a gene as the 2000 base pairs upstream of that gene. 
# Find the promoter regions for all genes in the genome.
# What is the range for the promoter of the 100th gene in g?
length(g)
promo_reg = flank(g, width = 2000)
range(promo_reg[100])

# Question 7
# Use findOverlaps() to find which GM12878 ESRRA binding sites overlap with promoter regions of genes.
# How many overlaps are there between GM12878 ESRRA binding sites and promoters?
findOverlaps(promo_reg,GM12878)

# How many unique GM12878 ESRRA binding sites overlap with promoters?
sum(GM12878 %over% promo_reg)
# How many unique promoters overlap with GM12878 ESRRA binding sites?
sum(promo_reg %over% GM12878)
# Or
length(unique(subjectHits(findOverlaps(GM12878, promo_reg))))

# Question 8 
# You can compare ESRRA binding in promoter regions to ESRRA binding in random non-promoter regions of a similar 
# length to determine whether ESRRA is more likely to bind promoters than control regions.
# Shift the promoters by 10000 bases and find the number of overlaps between GM12878 ESRRA binding sites and 
# the shifted regions.
random_shift = shift(promo_reg,10000)

# How many overlaps are there between GM12878 ESRRA binding sites and the shifted regions?
findOverlaps(random_shift,GM12878)

# How many unique GM12878 ESRRA binding sites overlap the shifted regions?
sum(GM12878 %over% random_shift)

# What is the ratio of the number of unique ESRRA binding sites overlapping promoters versus the number of unique 
# ESRRA binding sites overlapping shifted regions in GM12878 cells? That is, how much more likely is an ESRRA binding 
# site to be found in a promoter than a shifted region?
sum(GM12878 %over% promo_reg)/sum(GM12878 %over% random_shift)

