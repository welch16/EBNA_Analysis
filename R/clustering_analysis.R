
pam_silhouette_dt <- function(mat,k,...)
{
  clust = pam(mat,k,...)
  sil = silhouette(clust)
  sil_dt = data.table(names = row.names(sil),cluster =sil[,1],
    neighbor =sil[,2],width=sil[,3])
  return(sil_dt)
}
