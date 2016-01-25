fitTechNoise <- function(spikes, plot=FALSE, ...) 
# Fits a trend to the technical noise, i.e., the quarter-root variance (see the voom paper)
# of the log-CPMs against the average log CPM.
{
	my.var <- apply(log2(t(t(spikes+0.5)/colSums(spikes))), 1, FUN=var)^0.25
	covar <- edgeR::aveLogCPM(spikes)
	fit <- limma::loessFit(y=my.var, x=covar) 
	if (plot) {
		plot(covar, my.var, xlab="Average abundance", ylab="Quarter-root variance")
		o <- order(covar)
		lines(covar[o], fit$fitted[o], col="red", lwd=2)
	}
	FUN <- approxfun(x=covar, y=fit$fitted, rule=2)
	function(x) { FUN(x)^4 }
}

decomposeNoise <- function(gene.data, spike.data, FUN) 
# Decomposes the variance in the log-CPMs for the actual genes.
{
	spike.sizes <- colSums(spike.data)
	cur.data <- log2(t(t(gene.data+0.5)/spike.sizes))
	ave.ab <- edgeR::aveLogCPM(gene.data, lib.size=spike.sizes)
	tech.var <- FUN(ave.ab)
	total.var <- apply(cur.data, 1, FUN=var)
	bio.var <- pmax(0, total.var - tech.var)
	list(adjc=cur.data, ab=ave.ab, total=total.var, bio=bio.var, tech=tech.var)
}
