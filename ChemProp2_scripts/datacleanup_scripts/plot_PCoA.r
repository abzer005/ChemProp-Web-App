plotPCoA <- function(ft, 
                     md, 
                     distmetric = c("bray","euclidean", "maximum", "manhattan","canberra", "binary", "minkowski"), 
                     title = 'Principal coordinates plot'){

    distm <- vegdist(ft, method = distmetric)
    d <- as.matrix(distm)
    
    #computing multi-dimensional scaling on distance matrix
    raw_pcoa <- cmdscale(distm,
                     k = 10, #k=10 gets first 10 principal coordinates
                     eig = T, 
                     add = T)  

    pcoa_pts <- as.data.frame(raw_pcoa$points) #getting the principal coordinates (or PCOs)
    var <- round(raw_pcoa$eig*100/sum(raw_pcoa$eig),1) # getting the variance explained by each PCo
    names(pcoa_pts)[1:10] <- paste0("PCoA",seq(1,10)) #naming the PCos
    
    identical(rownames(pcoa_pts),md$filename)
    
    display(data.frame(ID=c(1:ncol(md)), Column_Names_Metadata=colnames(md))) #prints the metadata column names
    
    #choose the metadata column to which PCoA should be displayed
    metadata_column <- as.numeric(readline("Enter the number corresponding to the metadata column you want to use for PCoA: "))
    
    resultPlot <- ggplot(pcoa_pts, aes(x = PCoA1, 
                                  y = PCoA2, 
                                  color = as.factor(md[,metadata_column]))) + 
        geom_point(size=3) +
        labs(color = colnames(md)[metadata_column]) + 
        ggtitle(title) + 
        xlab(paste('PCo1',var[1],'%', sep = ' ')) + 
        ylab(paste('PCo2',var[2],'%', sep = ' ')) 

  return(resultPlot)
}
