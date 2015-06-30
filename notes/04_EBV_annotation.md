
## EBNA annotation

Previously we annotated the curated EBV (EBNA2, EBNA3's and RBPJ)
peaklists to the following regions in the genome:

- gene - regions in the gene body of the genome

- 5p1 and 3p1 - regions such that the distance to the nearest genes to
  the 5' or 3' directions of the genes respectively are between 1 and 1,999 bp

- 5p2 and 3p2 - similarly as above but the distance is between 2k and 9,999 bp

- 5d and 3d - similar and the distance are between 10k and 100k

- gd - The rest of the genome



Trying not to make a different annotation, we partitioned the genome
into the following regions:

- Promoters - which are the region that are 2k bp upstream from the
  reduced intervals in the gene body defined with the genes in [the
  last hg19 gencode
  annotation](http://www.gencodegenes.org/releases/19.html). This
  correspond to 5p1 in the annotation above.

- Gene Body - this are regions in the genome that overlaps any of the
  genes in the last [hg19 gencode
  database](http://www.gencodegenes.org/releases/19.html), but don't
  overlap any of the regions previouly defined as promoters. This
  correspond to gene in the annotation above.
 
- Integenic regions - Rest of the genome. This corresponds to the rest



The whole genome frequencies were estimated as the ratio between the
number of bp in region (as defined above) and genome length, and the
region probabilities are defined as the number of labelled peaks
divided by total number of peaks following:

- Promoters - 5p1
- Gene body - gene
- Intergenic - rest


<img src="../figures/EBV_annots/whole -1.png" title="plot of chunk whole " alt="plot of chunk whole " width="600" />

The estimated probabilites are given as:


```
## Error in kable(dcast.data.table(sample ~ lab, data = whole, value.var = "probs"), : table format not implemented yet!
```

Clearly all estimated probabilities are significantly different from
the whole genome probabilities:


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


