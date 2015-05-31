
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

The median silhouette index for each k behaves as:

![plot of chunk median_sil](../figures/pam_clusters/median_sil-1.png) 



In this case, since we see a decreasing tendency in the silhouette
index, and we want to consider a greater selection of clusters than
two, we consider K=6. So, we consider that amount of clusters and
calculate the number proportions of peaks overlapping a specific
transcription factor per cluster.





![plot of chunk propotions](../figures/pam_clusters/propotions-1.png) ![plot of chunk propotions](../figures/pam_clusters/propotions-2.png) ![plot of chunk propotions](../figures/pam_clusters/propotions-3.png) ![plot of chunk propotions](../figures/pam_clusters/propotions-4.png) ![plot of chunk propotions](../figures/pam_clusters/propotions-5.png) 










