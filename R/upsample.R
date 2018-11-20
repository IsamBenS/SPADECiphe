# Transpose data before call to put in row major order
SPADE.assignToCluster <- function(tbl, cluster_data, cluster_assign)
	.Call("SPADE_assign",t(tbl),t(cluster_data),as.integer(cluster_assign))

SPADE.addClusterToFCSCIPHE <- function (
	
	in.fcs,
	
	clusters.fcs, 
	
	cols = NULL, 
	
	arcsinh_cofactor = NULL, 
	
	transforms=flowCore::linearTransform(a=1),
	
	comp = TRUE
){
	
	
	if (!is.null(arcsinh_cofactor)) 
	{
		warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
		transforms <- flowCore::arcsinhTransform(a = 0, b = 1/arcsinh_cofactor)
	}
	in_fcs <- in.fcs
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
	cluster_fcs <- clusters.fcs
	cluster_data <- exprs(cluster_fcs)
	cluster_params <- parameters(cluster_fcs)
	cluster_pd <- pData(cluster_params)
	c_idxs <- match(cols, cluster_pd$name)
	na.fail(c_idxs)
	assign <- SPADE.assignToCluster(SPADE.transform.matrix(in_data[,idxs], transforms), 
									SPADE.transform.matrix(cluster_data[,c_idxs], transforms), 
									cluster_data[, "cluster"])
	
	channel_number <- ncol(in_fcs) + 1
	channel_id <- paste("$P", channel_number, sep = "")
	channel_name <- "cluster"
	channel_range <- max(assign) + 1
	plist <- matrix(c(channel_name, channel_name, channel_range, 
					  0, channel_range - 1))
	rownames(plist) <- c("name", "desc", "range", "minRange", 
						 "maxRange")
	colnames(plist) <- c(channel_id)
	pd <- rbind(pd, t(plist))
	pData(params) <- pd
	out_data <- cbind(in_data, cluster = assign)
	out_frame <- flowFrame(out_data, params, description = description(in_fcs))
	keyval <- list()
	keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
	keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
	keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
	keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
	keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
	keyword(out_frame) <- keyval
	
	return(out_frame)
}
