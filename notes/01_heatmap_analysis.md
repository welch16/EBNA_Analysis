
### Heatmaps analysis

For this analysis we are considering the rows after already filtering
the peaks that overlap any of the blacklisted regions and the columns
are already grouped by transcription factors.



#### Filtering to peaks that overlap Dnase hypersensitive sites

We generate the matrices and proportion vectors:


```r
mats <- mcmapply(build_binary_matrix, names(ranges),ranges,
  MoreArgs = list(expr = "Dnase == 1",
             col_expr = cols_to_remove(pattern,'columns')),
			 SIMPLIFY=FALSE,mc.silent=TRUE,mc.cores = mc)

proportions <- lapply(mats,get_proportions)
```



![plot of chunk hm_dnaseFilter ](../figures/heatmaps/hm_dnaseFilter -1.png) ![plot of chunk hm_dnaseFilter ](../figures/heatmaps/hm_dnaseFilter -2.png) ![plot of chunk hm_dnaseFilter ](../figures/heatmaps/hm_dnaseFilter -3.png) ![plot of chunk hm_dnaseFilter ](../figures/heatmaps/hm_dnaseFilter -4.png) ![plot of chunk hm_dnaseFilter ](../figures/heatmaps/hm_dnaseFilter -5.png) 

If we calculate the proportions of each TF, and sort them we can get the following:

![plot of chunk proportions](../figures/heatmaps/proportions-1.png) ![plot of chunk proportions](../figures/heatmaps/proportions-2.png) ![plot of chunk proportions](../figures/heatmaps/proportions-3.png) ![plot of chunk proportions](../figures/heatmaps/proportions-4.png) ![plot of chunk proportions](../figures/heatmaps/proportions-5.png) 










