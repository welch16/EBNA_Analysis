
## Building the incidence matrix

We considered a list of all the genomic intervals in the experiment as
the disjoint union of all the peaks in the EBNA2, EBNA3A, EBNA3B,
EBNA3C and RBPJ peak sets.

For each of the peaks in the list we built an incidence vector
representing which peak sets were used to construct it. That way if an
interval in the experiment is made as the union of an EBNA2 and an
EBNA3B peaks, then it's incidence vectors is going to be (1,0,1,0,0).

Additionally we built the ENCODE TF peaks overlaps indicators as 1
when any of the peaks used to construct the region overlaps the TF and
0 otherwise.



