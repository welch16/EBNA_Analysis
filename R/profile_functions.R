
load_clean_profile <- function(file ,window_ext)
{
  load(file) ## profile
  profile[,match := paste(chr,match,sep = "_")]
  profile[,chr := NULL]
  profile <- profile[order(match,coord)]
  npeaks <- nrow(profile) / (2 * window_ext + 1)
  profile[,coord := rep(-window_ext:window_ext, npeaks)]
  return(profile)
}


grepV <- function(pattern,x,ignore.case = FALSE,perl = FALSE,value = FALSE,
  fixed = FALSE,useBytes = FALSE, invert = FALSE)
{
  return(x[grep(pattern,x,ignore.case = ignore.case,perl = perl,value = value,
  fixed = fixed,useBytes = useBytes, invert = invert)])
}

add_names <- function(set,data,mc = 24)
{
    out <- mcmapply(function(x,y){   
      v <- strsplit(y,"_")[[1]]
      aux <- copy(x)
      aux[,set := v[1]]
      z <- strsplit(v[2],"-")[[1]]
      aux[,histone := z[1]]
      aux[,rep := as.numeric(z[2])]
      return(aux)
    },data[[set]],names(data[[set]]),SIMPLIFY = FALSE,mc.cores = mc)
  out <- do.call(rbind,out)
  return(out)
}

add_match <- function(set,unified_lists,factor_overlaps)
{
  peaks <- cbind(unified_lists$peaks,unified_lists$overlaps,factor_overlaps)
  chr <- unique(peaks[,(seqnames)])  
  setkey(peaks,seqnames)
  idx <- which(peaks[[set]] == 1) 
  lengths <- peaks[idx][,length(start),by = seqnames]
  peaks <- peaks[idx]
  peaks[,match := ""]
  setkey(lengths,seqnames)
  matchs <- lengths[,paste(seqnames,1:V1,sep = "_"),by = seqnames][,(V1)]
  peaks[order(seqnames,start),match := matchs]
  return(peaks)
}
