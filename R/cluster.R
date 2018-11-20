# Cluster observations into ~k clusters
SPADE.cluster <- function(tbl, k) {
	if (nrow(tbl) > 60000) {
		warning("Potentially too many observations for the clustering step",immediate=TRUE);
	}

	if (nrow(tbl) < k) {
		stop("Number of requested clusters exceeds number of events")
	}

	# Transpose table before call into row major order
	cluster = Rclusterpp.hclust(tbl);
	clust = list(assgn=cutree(cluster,k=k));

	# Invalid clusters have assgn == 0
	centers = c()
	is.na(clust$assgn) <- which(clust$assgn == 0)
	for (i in c(1:max(clust$assgn, na.rm=TRUE))) {  
		obs <- which(clust$assgn == i)
		if (length(obs) > 1) {
			centers <- rbind(centers,colMeans(tbl[obs,,drop=FALSE]))
			clust$assgn[obs] <- nrow(centers)
		} else {
			is.na(clust$assgn) <- obs
		}
	}
	return(list(centers=centers,assign=clust$assgn,hclust=cluster))
}

SPADE.clustersToMST <- function(centers, method="manhattan") {
	adjacency  <- dist(centers, method=method)
	full_graph <- graph.adjacency(as.matrix(adjacency),mode="undirected",weighted=TRUE)
	mst_graph  <- minimum.spanning.tree(full_graph)
	mst_graph
}

SPADE.writeGraph <- function(graph, outfilename) {
	 write.graph(graph, outfilename, format="gml")
}

SPADE.FCSToTreeCIPHE <- function (
	
	in.fcs.files, 
	
	cols = NULL, 
	
	k = 200, 
	
	arcsinh_cofactor = NULL, 
	
	transforms=flowCore::linearTransform(a=1),
	
	desired_samples = 50000, 
	
	comp = TRUE
) {
	
	
	out.data <- list(NULL,NULL,NULL, NULL)
	if (!is.null(arcsinh_cofactor)) 
	{
		warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
		transforms <- flowCore::arcsinhTransform(a = 0, b = 1/arcsinh_cofactor)
	}
	data = c()
	in.fcs.files.list <- in.fcs.files
	if(!is.list(in.fcs.files))
	{
		in.fcs.files.list <- list(in.fcs.files)
	}
	for (f in in.fcs.files.list) 
	{
		in_fcs <- f
		in_data <- exprs(in_fcs)
		params <- parameters(in_fcs)
		pd <- pData(params)
		if (is.null(cols)) 
		{
			cols <- as.vector(pd$name)
		}
		idxs <- match(cols, pd$name)
		if (any(is.na(idxs))) 
		{
			stop("Invalid column specifier")
		}
		data <- rbind(data, in_data[, idxs, drop = FALSE])
		colnames(data) <- pd$name[idxs]
	}
	if (nrow(data) > desired_samples) 
	{
		data <- data[sample(1:nrow(data), desired_samples), ]
	}
	else if (nrow(data) == 0) 
	{
		stop("Number of observations to cluster is zero. Did the density/downsampling warn about data similarity?")
	}
	clust <- SPADE.cluster(SPADE.transform.matrix(data, transforms),  k)
	
	temp.path <- tempfile(patter = "merge_order", fileext = "txt")
	write.table(clust$hclust$merge, file = temp.path, sep = "\t", 
				quote = F, row.names = F, col.names = F)
	out.data[[4]] <- read.table(temp.path,sep = "\t")
	
	
	temp.fcs.sp <- SPADE.build.flowFrame(subset(cbind(data, cluster = clust$assign), 
												  !is.na(clust$assign)))
	temp.fcs.sp.path <- tempfile(fileext = "fcs")
	write.FCS(temp.fcs.sp, temp.fcs.sp.path)
	out.data[[1]] <- read.FCS(temp.fcs.sp.path)
	
	
	
	graphfilename <- tempfile(pattern="mst",fileext = "gml")
	SPADE.writeGraph(SPADE.clustersToMST(clust$centers), graphfilename)
	out.data[[3]] <- read.gml(graphfilename)
	
	
	clusterfilename <- tempfile(pattern="clusters",fileext = "table")
	write.table(clust$centers, file = clusterfilename, row.names = FALSE, 
				col.names = colnames(data))
	out.data[[2]] <- read.table(clusterfilename, col.names = colnames(data))
	
	return(out.data)
}
