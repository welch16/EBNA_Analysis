
### Partition around medoids analysis

Following the heatmap analysis where we were able to asses that there
is in fact some structure in the data. We considered using the
partition around medoids clustering strategy, where the only parameter
that is needed to tune is the number of clusters. Therefore, we
clustered the data considering k between 2 and 15.








If we make boxplots of the silhouette index for each point and
separate them by number of clusters and dataset, we can see that the
biggest separation is occuring when there are two clusters, which
agrees with the heatmaps. When using a higher number of clusters we
can see that the range of the boxplots stabilizes:


![plot of chunk silhouette_idx ](../figures/pam_clusters/silhouette_idx -1.png) 

```
## Error in scale_color_brewer(palette = "Set1") + theme(legend.position = "none") + : non-numeric argument to binary operator
```

The median silhouette index for each k behaves as:


```
## Error in eval(expr, envir, enclos): could not find function "scale_x_continous"
```










