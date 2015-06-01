
pam_silhouette_dt <- function(mat,k,...)
{
  clust = pam(mat,k,...)
  sil = silhouette(clust)
  sil_dt = data.table(names = row.names(sil),cluster =sil[,1],
    neighbor =sil[,2],width=sil[,3])
  return(sil_dt)
}



top_k <- function(V1,variable,k,decreasing = TRUE)
{
  names(V1) = as.character(variable)
  u = sort(V1,decreasing = decreasing)
  u = u[1:k]
  nm = names(u)
  names(u) = NULL
  return(list(u,nm))
}

  

top_strategy <- function(peak,props,top)
{
  top_TF <- props[set == peak, top_k(V1,variable,top),by=.(set,clustering)]
  TFs <- unique(top_TF[,(V2)])
  return(TFs)
}

cut_strategy <- function(peak,props,cutline)
{
  cut_TF <- props[set == peak & V1 >= cutline]
  TFs <- as.character(unique(cut_TF[,(variable)]))
  return(TFs)
}

build_integrative_mat <- function(prop)
{
  cols <- unique(prop[,(variable)])
  prop_loc <- copy(prop)
  prop_loc[,rows:= paste0(set,"-",clustering)]
  rows <- unique(prop_loc[,(rows)])
  mat <- lapply(cols,function(x,rows,prop_loc)
     ifelse(rows %in% prop_loc[variable == x,(rows)],1,0),rows,prop_loc)
  mat <- as.matrix(do.call(cbind,mat))
  colnames(mat) <- cols
  rownames(mat) <- rows
  return(mat)
}
