
build_binary_matrix <- 
function(peakSet,ranges,expr = NULL,col_expr = NULL,rem.Col.zero = TRUE  )
{
   # Returns a binary matrix built from the ranges objects
   columns = get_names(ranges)
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
     sub = subset(ranges,select = columns,
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

gg_heatmap <- function(mat,row.names = FALSE,title = NULL)
{

  hh = hclust(dist(mat))
  vv = hclust(dist(t(mat)))

  hh1 = as.dendrogram(hh)
  vv1 = as.dendrogram(vv)
  
  h_order = order.dendrogram(reorder(hh1,rowMeans(mat)))
  v_order = order.dendrogram(reorder(vv1,colMeans(mat)))

  row_names = rownames(mat)
  col_names = colnames(mat)

  df = melt(mat)
  names(df) = c("Var1","Var2","value")
  
  df$Var1 = factor(df$Var1,levels = row_names[h_order])
  df$Var2 = factor(df$Var2,levels = rev(col_names[v_order]))

  p1 <- ggplot(df,aes(Var2,Var1))+geom_tile(aes(fill= factor(value)))+
    scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0))+
      labs(x="",y="")+theme(legend.position = "bottom")+
      theme(axis.text.x = element_text(angle = -90,size= 6.2,hjust= 0),
            axis.text.y = element_text(size = 6.2),
            plot.margin = unit(c(0,1,.5,1),"cm"))+
      scale_fill_manual(name = "Overlap",values = c("lightblue","navyblue"),labels = c("no","yes"))+
      scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
  
  if(!row.names){
    p1 <- p1 + theme(axis.ticks = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
  }
  
  ddata_y <- dendro_data(vv1)

  ### Set up a blank theme
  theme_none <- theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_text(colour=NA),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank()     
      )
  

  # Dendrogram 2
  p3 <- ggplot(segment(ddata_y)) +
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
      theme_none+labs(x=NULL,y=NULL)+
      theme(plot.margin = unit(c(.5,0,0,0),"cm"))+
      scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))+
        theme(axis.ticks.y = element_blank())

  if(!is.null(title)){
    p3 <- p3 + ggtitle(title)
  }
  
  return(list(p1,p3))  
  
}

print_hm <- function(x)
{
  p1 <- x[[1]]
  p2 <- x[[2]]

  g <- ggplotGrob(p1)

  g <- gtable_add_rows(g, unit(1.2,"in"), 0)
  g <- gtable_add_grob(g, ggplotGrob(p2) ,
         t = 1, l=4, b=1, r=4)

  grid.newpage()
  grid.draw(g)
 
}
