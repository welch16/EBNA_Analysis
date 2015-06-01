
## Integrative cluster analysis



In the previous analysis we considered clustering the peak set
specific matrices filtering the regions that didn't overlap DHS and
then we calculated the frequency of overlapping peaks in a cluster
with a specific TF. That way, we considered two strategies to
integrate the clusters of the different peak sets.






That way we can build two different incidence matrices.

1. One where the rows represent one of the 10 clusters of each
   dataset (i.e. 50 rows) and the columns are the
   union of the top 10 TFs for each peaksets.
  
2. Again the rows represent one of the 10 clusters of each datasets
   (excepting the clusters where there are not TFs such that its
   overlap frequency is greater 75%). And the columns
   represent the union of all TFs such that its overlapping frequency
   respect the peaks in the cluster is greater than 75%.




   
### Using the union of the top 10 TF per cluster as columns

![plot of chunk hm_top ](../figures/cluster_integrative/hm_top -1.png) 

- This heatmap divides the drivers as:
  * EBF1, FOXM1, P300, BCL11A, PAX5, RUNX3, NFIC and ATF2
  * POL2, CTCF, PU1
  * TAF1, PML, MXI1, CHD2, ELF1, MAZ, YY1
  * IRF4, BATF, MTA3, MEF2A

Roughly we can see the rows are partitioned into four groups,
therefore we can see the cluster presence frequency for the columns of
the matrix:



![plot of chunk top_propotions](../figures/cluster_integrative/top_propotions-1.png) 



  

### Using the union of the TFs with clustered frequency greater than 75

  
   
![plot of chunk hm_cut ](../figures/cluster_integrative/hm_cut -1.png) 

- This one divides the drivers as:
  * EBF1, RUNX3, PAX5, BCL11A, ATF2, NFIC, BATF, P300
  * MAX, TAF1, ELF1, POL2, YY1, MAZ, MXI1
  * SP1, PML, NFATC1, TBLR1, MEF2A, BCL3, MTA3

In this heatmaps we can see that the rows are divided into 3 segments:
  


![plot of chunk cut_propotions](../figures/cluster_integrative/cut_propotions-1.png) 



