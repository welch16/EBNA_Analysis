
## A brief summary

The purpose of this document is to outline the analysis made without
entering into specific details.

### Input files

For the analysis we used  as inputs the following datasets (everything
on the GM12878 cell line):

- peaks called for EBNA2, EBNA3A, EBNA3B, EBNA3C and RBPJ (JK92 and
  JK234) using mosaics 
  
- Peaks of transcription factors called by ENCODE's pipeline spp
  optimal (June 2012 series)

  * We also considered 
  

- Dnase Hypersensitive sites, for this there are three different
  samples, defined an overlap with a DHS if there is an overlap with a
  site in any of the 3 DHS files

- Histone reads files, for this we are considering the following histones:

  * H3k27me3
  * H3k4me3
  * H3k9ac
  * H3k9me3
  * H3k36me3

### First steps of processing

Using ENCODE's hg19 blacklist, we removed original peaks that
overlapped any of the blacklisted regions. For each of the remaining
peaks an indicator matrix was built.

Using this matrices, heatmaps were calculated where the clustering
algorithm used was the default one in R, i.e. euclidean distance and
complete agglomeration (the distance between two clusters is update
with the max distance)

Similarly heatmaps we generated using those columns and overlaps with
histone peaks and Dnase hypersensitive sites.
		

To summarize the columns (since there are several replicates for each
transcription factor), a list with the relationship between TF and its
replicates was manually curated. Then this analysis was repeated using
as columns the TF's considered in the curated list. A peak is said to
overlap a TF if if overlaps any of its replicates.

To summarize the rows, we used partition around medoids with 10
clusters. We followed two different strategies:

- For each cluster, get the top M TF's. Then built an incidence matrix
  such that each entry was 1 if any of those top M transcription
  factors were present in the cluster and 0 otherwise.

- For each cluster, we considered the local proportion of each
  TF. Then an incidence matrix was built where each entry is defined
  as 1 if the local proportion of a TF in the cluster in greater than
  the threshold and 0 otherwise.

This procedures were repeated several times considering the presence
of Dnase or histone columns in the incidence matrices.

### Proportion plots

After exploring the data using hierarchical clustering, partition
around medoids was used to further improve the analysis. The number of
clusters were selected using the silhouette index and for each cluster
the incidence proportion (number of overlaping peaks with TF peaks /
total number of peaks).


