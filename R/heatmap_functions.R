
build_binary_matrix <- 
function(peakSet,ranges,expr = NULL,col_expr = NULL,rem.Col.zero = TRUE  )
{
   # Returns a binary matrix built from the ranges objects
   columns = get_names(ranges[[peakSet]])
   columns = subset(columns,eval(parse(text = grepl_cmd("summit","columns"))))
   columns = subset(columns,eval(parse(text = grepl_cmd("region","columns"))))
   columns = subset(columns,eval(parse(text = grepl_cmd("allTF","columns"))))
   columns = subset(columns,eval(parse(text = grepl_cmd("GeneID","columns"))))
   columns = subset(columns,eval(parse(text = grepl_cmd("distance","columns"))))
   columns = subset(columns,eval(parse(text = grepl_cmd(peakSet,"columns"))))
   columns = subset(columns,eval(parse(text = grepl_cmd("JK","columns"))))
   if(!is.null(col_expr)){
     columns = subset(columns,eval(parse(text = col_expr)))
   }
   if(!is.null(expr)){
     sub = subset(ranges[[peakSet]],select = columns,
       subset = eval(parse(text = expr))) 	
   }else{
     sub = subset(ranges[[peakSet]],select = columns)
   }
   mat = as.matrix(sub@elementMetadata)
   colnames(mat) = columns
   if(rem.Col.zero) mat = mat[,colSums(mat)>0]
   rownames(mat) = paste(seqnames(sub),start(sub),sep = " - ")
   mat = mat[rowSums(mat)> 0,]
   return(mat)
}

get_proportions <-
function(mat)
{
  return(sort(colSums(mat) / dim(mat)[1],decreasing = TRUE))
}
