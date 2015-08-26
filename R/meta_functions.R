
load_meta <- function(dir)
{
  files <- list.files(dir)
  meta <- files[grep("meta",files)]
  out <- data.table(read.table( file.path(dir,meta),sep = "\t",stringsAsFactors = FALSE,
                               header = TRUE))
  return(out)
}


load_files <- function(meta,dir,mc)
{  
  files <- paste0(meta[,(File.accession)],".bed.gz")
  bed_files <- mclapply(file.path(dir,files),read.table,
    stringsAsFactors = FALSE,header = FALSE,mc.cores = mc)
  gr <- lapply(bed_files,function(x){
    out <- GRanges(seqnames = x$V1,
                   ranges = IRanges(start = x$V2,end = x$V3),
                   strand = "*")
    return(out)})
  return(gr)
}
