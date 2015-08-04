
## EBNA annotation

Previously we annotated the curated EBV (EBNA2, EBNA3's and RBPJ)
peaklists to the following regions in the genome:

- gene - regions in the gene body of the genome

- 5p1 and 3p1 - regions such that the distance to the nearest genes to
  the 5' or 3' directions of the genes respectively are between 1 and 1,999 bp

- 5p2 and 3p2 - similarly as above but the distance is between 2k and 9,999 bp

- 5d and 3d - similarly as abpve but the distance is between 10k and 100k

- gd - The rest of the genome



Trying not to make a different annotation, we partitioned the genome
into the following regions:

- Promoters - which are the region that are 2k bp upstream from the
  reduced intervals in the gene body defined with the genes in [the
  last hg19 gencode
  annotation](http://www.gencodegenes.org/releases/19.html)

- Gene Body - this are regions in the genome that overlaps any of the
  genes in the last [hg19 gencode
  database](http://www.gencodegenes.org/releases/19.html), but don't
  overlap any of the regions previouly defined as promoters
 
- Integenic regions - Rest of the genome		



The whole genome frequencies were estimated as the ratio between the
number of bp in region (as defined above) and genome length, and the
region probabilities are defined as the number of labelled peaks
divided by total number of peaks following

- Promoters - 5p1
- Gene body - gene
- Intergenic - the remaining categories


<img src="../figures/EBV_annots/whole -1.png" title="plot of chunk whole " alt="plot of chunk whole " width="600" />

The estimated probabilites are given as:


|sample | Gene body| Intergenic| Promoter|
|:------|---------:|----------:|--------:|
|genome |    0.2697|     0.7147|   0.0156|
|EBNA2  |    0.4871|     0.4298|   0.0831|
|EBNA3A |    0.4107|     0.4028|   0.1866|
|EBNA3B |    0.4690|     0.4405|   0.0905|
|EBNA3C |    0.4298|     0.5464|   0.0238|
|JK92   |    0.4996|     0.4359|   0.0645|
|JK234  |    0.4765|     0.4884|   0.0352|
|RBPJ   |    0.4886|     0.4547|   0.0567|

We can see that all estimated probabilities are significantly
different from the whole genome probabilities:


```
## $EBNA2
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 5114.39, df = 2, p-value < 2.2e-16
## 
## 
## $EBNA3A
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 3430.394, df = 2, p-value < 2.2e-16
## 
## 
## $EBNA3B
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 1862.536, df = 2, p-value < 2.2e-16
## 
## 
## $EBNA3C
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 501.6196, df = 2, p-value < 2.2e-16
## 
## 
## $JK92
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 3737.446, df = 2, p-value < 2.2e-16
## 
## 
## $JK234
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 1092.665, df = 2, p-value < 2.2e-16
## 
## 
## $RBPJ
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 3807.25, df = 2, p-value < 2.2e-16
```

### Conditional on Dnase

Following the previous analysis, we repeated it but considering to
estimate the promoter and gene probabilities as the number of base
pairs in the union of the Dnase peaks that satisfy the criteria above
for promoters and genes respectively.






<img src="../figures/EBV_annots/whole_dnase -1.png" title="plot of chunk whole_dnase " alt="plot of chunk whole_dnase " width="600" />

And the estimated probabilities conditional on Dnase are:

The estimated probabilites are given as:


|sample | Gene body| Intergenic| Promoter|
|:------|---------:|----------:|--------:|
|genome |    0.3557|     0.5633|   0.0811|
|EBNA2  |    0.4901|     0.4190|   0.0909|
|EBNA3A |    0.4068|     0.3098|   0.2834|
|EBNA3B |    0.4812|     0.4164|   0.1024|
|EBNA3C |    0.4614|     0.4986|   0.0400|
|JK92   |    0.5062|     0.4243|   0.0695|
|JK234  |    0.4990|     0.4585|   0.0425|
|RBPJ   |    0.5029|     0.4318|   0.0653|


```
## $EBNA2
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 702.8648, df = 2, p-value < 2.2e-16
## 
## 
## $EBNA3A
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 665.3966, df = 2, p-value < 2.2e-16
## 
## 
## $EBNA3B
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 234.2236, df = 2, p-value < 2.2e-16
## 
## 
## $EBNA3C
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 107.4833, df = 2, p-value < 2.2e-16
## 
## 
## $JK92
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 733.9508, df = 2, p-value < 2.2e-16
## 
## 
## $JK234
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 323.594, df = 2, p-value < 2.2e-16
## 
## 
## $RBPJ
## 
## 	Chi-squared test for given probabilities
## 
## data:  counts
## X-squared = 791.8673, df = 2, p-value < 2.2e-16
```
