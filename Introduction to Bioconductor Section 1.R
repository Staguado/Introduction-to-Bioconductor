####################################################################################################################

# By: Santiago Taguado Menza
# Introduction to Bioconductor: 
# April 14th, 2021
# HarvardX

####################################################################################################################

# Getting Started with Bioconductor:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::version()

# install Bioconductor packages
# The code below
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install(c("genefilter","geneplotter"))

# load installed packages
library(BSgenome.Hsapiens.UCSC.hg19)
library(genefilter)
library(geneplotter)
library(Biobase)    
methods(class = ExpressionSet)

# vignettes teach you how to use various functions in a package
vignette(package="Biobase")
vignette("ExpressionSetIntroduction")
browseVignettes(package="Biobase")

# report key details about your R session
sessionInfo()

BiocManager::install(c("genefu",
                       "COPDSexualDimorphism",
                       "gwascat",
                       "hgu133a.db",
                       "genomicsclass/tissuesGeneExpression"))

# Mammaprint Signature Genes: 

# Try using your basic R skills to explore the genes used in the Mammaprint gene signature described 
# in the video. This diagnostic signature assesses the risk that a breast cancer will progress using 
# the gene expression levels of 70 genes.

# You can find this signature in the genefu package, a Bioconductor package that provides information and 
# functions relevant for gene expression analysis, with a particular emphasis on breast cancer. 
# (Like all packages used in this section, this was installed with BiocManager::install()  
# on the Section 1 Overview page. You can install missing packages at any time with BiocManager::install().)

# Information on the 70 gene signature used in the Mammaprint algorithm is in the sig.gene70 data.frame.  
# You can have a look at this:
  
BiocManager::install("genefu")
library(genefu)
data(sig.gene70)
dim(sig.gene70)
head(sig.gene70)[,1:6]

# You can see from this that there are 70 records in the data frame, and that there are diverse ways of 
# describing the "genes" in the signature.

# Question 1
# You can read about the development of the 70 gene signature used in MammaPrint on Wikipedia:
# https://en.wikipedia.org/wiki/MammaPrint#Development
# How was the MammaPrint 70 gene signature designed?

# The 70 genes most statistically associated with breast cancer progression in microarray studies were chosen.

# Question 2
# How many components of the signature have a missing value (NA) for the associated NCBI gene symbol?

sum(is.na(sig.gene70$NCBI.gene.symbol))

# Question 3 
# What is the NCBI gene symbol for the gene with the Description "cyclin E2"?
sig.gene70[which(sig.gene70$Description == c("cyclin E2")),]

# Question 4
# Kinases are important for cell-cell communications; see the Wikipedia entry on Kinase for some background.
# The grep() function takes a pattern and a vector and returns the indexes of the vector that 
# match the pattern. Remember that you can see documentation and learn how to use a function with ?grep.
# You can use grep() on the Description field of the sig.gene70 data frame to search for substrings of 
# long gene names.
# How many of the members of the 70-gene signature are genes coding for kinases?

sig.gene70[grep("kinase",sig.gene70$Description),]

####################################################################################################################

# Assessment Phenotypes

# In the video we talk about phenotypes. Here we show some examples of what we mean by phenotypes, 
# how they can be coded in R objects, and how we compute with them. 
# Later in the course we will perform analyses that statistically connect these phenotypes to measured 
# molecular outcomes. Here we explore the use of data frames to store phenoypes (columns) from several 
# individuals (rows).

# Load the COPDSexualDimorphism.data package. This package includes data that can be used to compare the incidence of 
# COPD and emphysema by gender and smoking status. Use the commands:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("COPDSexualDimorphism.data")
library(COPDSexualDimorphism.data)
browseVignettes("COPDSexualDimorphism.data")
data(lgrc.expr.meta)

# Chronic Obstructive Pulmonary Diseases are a group of diseases that cause airflow blockage and breathing-related
# problems. Within these, emphysema is a condition that involves damage to the walls of the alveoli of the lung

# Question 1
# What is the number of female participants in this study?:
# Try the table() function, which counts the number of elements with each value, on the gender column.

table(expr.meta$GENDER)

# Question 2
summary(expr.meta$pkyrs)

# What is the median of the distribution of pack years smoked in this cohort?
median(expr.meta$pkyrs)
  
# What is the mean number of pack years smoked in this cohort?
mean(expr.meta$pkyrs)

# What is the maximum number of pack years smoked in this cohort?
max(expr.meta$pkyrs)

# Question 3
# Quantile-quantile plots (QQ-plots) are a useful tool for evaluating the simularity of two distributions.
# One common use of QQ-plots is to check whether the distribution of a dataset is well-approximated by the
# normal distribution. You can review how to make normal QQ-plots in the textbook.

# True or False: The distribution of pack-years smoked is well-approximated by a Gaussian (normal) 
# probability distribution.

line = qnorm(seq(0.01,0.99,0.01))
qqplot(line,expr.meta$pkyrs)
qqline(expr.meta$pkyrs)

qqnorm(expr.meta$pkyrs)
qqline(expr.meta$pkyrs)

# Question 4
# The units with which a quantity is recorded are not always the most effective for statistical analysis.
# Use the command boxplot(pkyrs~gender, data=expr.meta) to examine the distributions of pack years by gender using the 
# boxplot statistics and outlier flagging.

# Which of the following is an aspect of the display that would suggest caution in using the t-test in 
# comparing males and females with respect to pack years smoked?

boxplot(pkyrs~gender, data=expr.meta)

# Specially, Distributions appear quite asymmetric, with long tails skewed towards high values.

# All of the propositions among the choices are true, but the validity of the t-test is compromised when applied
# to variables with skewed distributions. The test's validity does not depend on separation of median values
# or magnitude of quartile, so first and last choices are incorrect. 
# Difference in number of outliers flagged can be important, but the more systematic concern with 
# asymmetry of distributions is the target response for this problem.    

####################################################################################################################

# Assessment: Chromosomes & SNIPS

# Genome-wide association studies (GWAS) are a major tool of genetic epidemiologists. In a case-control 
# design, individuals with a specific disease (cases) are identified along with otherwise similar 
# individuals who are disease free (controls).  SNP chips or DNA sequencing is then used to obtain 
# individuals' genotypes for a large number of SNPs. The genotype distributions for all SNPs are compared 
# between cases and controls, and those SNPs exhibiting association with disease are investigated for 
# potential insight into disruption of gene regulation or gene function.

# The Bioconductor gwascat package includes information on a catalog of GWAS results. Load the 
# gwascat package and check the version of the GWAS catalog stored in GRCh37 (hg19) coordinates.

BiocManager::install("GenomeInfoDbData")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("Homo.sapiens")
BiocManager::install("gwascat")
library(gwascat)

data(ebicat37)

# Use table() on the CHR_ID column to count the number of GWAS hits on each chromosome. 
# Sort the chromosomes by decreasing number of hits so that the first element has the most GWAS hits.

(sort(table(ebicat37$CHR_ID),decreasing=TRUE))[1:15]

# Which chromosome has the most GWAS hits in the catalog? Use an integer.

sort(table(mcols(ebicat37)[,"DISEASE/TRAIT"]), decreasing=TRUE)[1:25]

####################################################################################################################

# Assessment: Molecular Biology Concepts

# In a typical somatic human cell, before DNA replication, how many DNA copies are there for most genes 
# (ignore X and Y chromosomes)?

# Humans are a diploid species, meaning the somatic cells typically contain two copies of the autosomal 
# (not X or Y) chromosomes, before S phase in which the chromosomes are duplicated.
# 2 Copies

# In a single cell, how many RNA copies of a given gene are there?

# It varies.
# The number of copies of RNA depends on how much the gene is being transcribed. “Housekeeping genes” such as 
# those used to make ribosomal or transfer RNA are transcribed at a high rate, while others, such as mRNA for some 
# transcription factors, are transcribed less frequently. Some genes are not transcribed at all in certain cell types.

# Which of the following is most likely to be inherited in all the cells of all your body across a 
# generation (e.g. from either of your parents to you)

# a SNP (single nucleotide polymorphism) or other variant in the DNA
# Genetic information, such as single nucleotide polymorphisms, have a chance of being transmitted across generations.
# mRNA transcripts and proteins such as transcription factors degrade over time and, most importantly, do not replicate 
# themselves. In other words, DNA (not proteins or RNA) is known as the main molecule of genetic inheritance. 
# (Side note: There are cases of mRNA and proteins being temporarily inherited, for example the mRNA which are 
# in the egg cell at the moment it is fertilized by the sperm cell.)

# As a yeast cell of species S. cerevisiae passes from cell-cycle phase G1 to phase S, a bud begins to form and enlarge.
# When we distinguish between yeast colonies in the pre-budding and budding conditions, we are comparing

# Phenotypic States 
# The existence or nonexistence of the bud is a low-dimensional characteristic of the organism.

####################################################################################################################

# Verified Assessment: Biology Background

# In the following set of questions, we will analyze a gene expression microarray dataset. 
# Along the way, we will use the dataset to answer questions about gene expression and phenotype.

# We will be using microarray data included in the tissuesGeneExpression package.
# Load the dataset:

library(tissuesGeneExpression)
data(tissuesGeneExpression)

# This code loads a matrix, e, containing gene expression measurements obtained using a microarray technology. 
# More specifically, the matrix contains (log-scale) intensity data for thousands of microarray probes 
# (rows of the matrix) and several samples (columns of the matrix). Using the probe IDs, we can map rows of 
# the matrix to genes.

head(e[,1:5])

# The code above also loads the tissue vector, which specifies the tissue type of each sample.

table(tissue)

# Later on, we will learn about Bioconductor structures that make it easy for us to work with these data (e and tissue) 
# as a single object. For now, let's see how we can work with these data in their current form to analyze 
# how gene expression relates to tissue phenotypes.

# Question 3
# Different data types can be used to answer different biological questions.
# Which of these studies can be performed with the experimental data in tissuesGeneExperiment?

# identify which known transcripts correlate with tissue phenotype

# Gene expression microarrays give expression profiles for a selection of known transcripts. 
# They do not measure SNPs or transcription factor binding, and unlike RNA-seq they cannot identify new RNA species 
# which are not included on the microarray.

# Question 4

# Look at the data for the feature with ID "209169_at". You can index the rows of e directly with this character 
# string. For example, mean(e["209169_at",]) is about 7.26. 
# The gene measured by this probe is differentially expressed by tissue type.

# Which two tissues have the highest mean expression of this gene?
# Condition the tissue vector on your tissue type to determine which columns to use.

mean((e["209169_at",tissue == "cerebellum"]))
mean((e["209169_at",tissue == "colon"]))
mean((e["209169_at",tissue == "endometrium"]))
mean((e["209169_at",tissue == "hippocampus"]))
mean((e["209169_at",tissue == "kidney"]))
mean((e["209169_at",tissue == "liver"]))
mean((e["209169_at",tissue == "placenta"]))

# Question 5
# Based on these data, what could you hypothesize about the most likely function of the gene measured by "209169_at"?
  
# The gene plays a role in the nervous system.

# Question 6
# The Affymetrix technology for mRNA abundance measurement is based on hybridization of highly processed biological 
# sample material to synthetic oligonucleotide "probes" that are on fixed positions of the microarray surface. 
# Bioconductor provides detailed information on the probe and array structure as published by Affymetrix.

# The probe IDs in rownames(e) are probe IDs from the Affymetrix Human GeneChip HGU133A microarray. 
# We can find the gene symbols corresponding to these probe IDs with this code 
# (we will learn more about how it works later in the course):
  
rownames(e)
library(hgu133a.db)# installed in section overview
symbol = mapIds(hgu133a.db, keys=rownames(e), column="SYMBOL", keytype="PROBEID")

# This creates a vector symbol with the gene symbols correponding to the probes in rownames(e).
# What gene symbol corresponds to the probe with ID "209169_at"?

symbol["209169_at"]
sum(symbol ==   "GPM6B" , na.rm = TRUE)

# Question 7

# We can use sum(symbol == "GAPDH", na.rm=TRUE) to count the number of array features that measure expression of 
# the gene GAPDH.
# How many features in this array measure expression of gene H2AFX?

sum(symbol == "GAPDH", na.rm=TRUE)
sum(symbol == "H2AX" , na.rm=TRUE)

# Question 8
# Which of these probes is not associated with the H2AFX gene?

symbol["205436_s_at"]
symbol["212524_x_at"]
symbol["212525_s_at"]
symbol["212526_s_at"]
symbol["213344_s_at"]

# Question 9
# Consider the following plot:
# Which of the following relationships are suggested by this plot?

boxplot(as.numeric(e["205436_s_at",])~tissue)

# Question 10 
# Below is a vector of 6 IDs which index features (rows) of the matrix 'e':
IDs = c("201884_at", "209169_at", "206269_at", "207437_at", "219832_s_at", "212827_at")
# Which of the following IDs appear to represent a gene specific to placenta?

library(rafalib)
mypar(2,3)  
myfun = function(x){
boxplot(as.numeric(e[x,])~tissue,
        main = x)
}  
sapply(IDs, myfun)
