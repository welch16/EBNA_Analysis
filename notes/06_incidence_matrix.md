
## Building the incidence matrix

We considered a two-step approach to construct a peak overlap
incidence matrix using the lists of peaks for EBNA2, EBNA3A, EBNA3B,
EBNA3C and RBPJ. Specifically,

1. We obtained the regions universe as the disjoint union of all the
   elements in the EBNA2, EBNA3A, EBNA3B, EBNA3C and RBPJ peak lists.

2. We built an incidence matrix by using the regions universe and the
   peak lists. Each entry of the incidence matrix was set to 1 when a
   peak in the list was used to construct the region and 0 otherwise.

Additionally, we constructed ENCODE transcription factors columns to
indicate whether or not the region overlapped a peak in the ENCODE
peak sets obtained from.
