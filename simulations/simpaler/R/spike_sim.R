###########################################################################
# This code descibes the simulation process for generating variable spike-in quantities from real data.
# The aim is to explore whether the results are sensitive to spike-in variability, under the assumption of constant spike-ins.
# We re-sample the added volume for each well, and we scale the spike-ins to reflect this new volume.
# The analysis of the scaled data is then compared to the observed analysis (which is assumed to be truth, i.e., with spike-in volumes equal across wells).
# This is sufficient to test sensitivity to spike-in variability, without the need to specify a complex simulation of all the individual genes.
#
# We assume that sequencing is performed to saturation in each data set.
# This simplifies the simulation as it means that we do not have to account for undersampling, i.e., \rho_i = 1.
# Variability in the added volume will then directly translate to variability in the total spike-in count.
# It also means that we do not have to rescale the cellular counts to reflect variable undersampling.
# Doing so would be problematic, as direct scaling of the counts would distort the empirical mean-variance relationship of the count data.
#
# Note that the spike-in quantities are *more* variable when undersampling is not considered.
# This is because the undersampling factor \rho_i includes the spike-in volume as a term in the denominator.
# As a result, the total spike-in count is an expression where the volume is divided by itself (plus a constant).
# One can intuitively see that this is less variable than the volume alone, as the self-division will shrink the final value towards a constant.
# We do not find this effect objectionable with respect to the intent of our simulations.
# If variability in the spike-in volume has no effect without undersampling, then it should have even less effect with undersampling.
# 
# On to the practicalities; we compute the fitted NB mean for each gene across all wells, treating the size-scaled value as the true abundance of the spike-in for each well.
# We then sample a log-scaling factor from a normal distribution with the estimated \sigma^2_{vol}.
# We use this to scale the fitted mean in each well to simulate the effect of a change in the added spike-in volume.
# The count is obtained by matching the quantiles of the old and new means (a la pseudo-counts in classic edgeR; though, for large dispersions, we might as well scale the counts).
# We don't perform re-simulation to avoid adding more shot noise than what is already present.
# Otherwise, we wouldn't be able to tell if variability in the results was due to spike-in variability, or due to the sampling scheme (we'd need to do controls to prepare for that).
#
# Note that the dispersion estimation isn't affected by spike-in variability, as we rescale all spike-ins to the same library size prior to estimating the dispersion.

###########################################################################
# The assumption of constant spike-in volumes is generally necessary to generate realistic simulation data.
# For example, the most obvious approach would be to simulate counts from scratch, using a simulated distribution with some parameters.
# To estimate realistic values for those parameters, you'd need to normalize before you fit a model to the data.
# This normalization step would make the implicit assumption that the spike-in volumes are the same.
# Here, we just skip the step of getting the data parameters and resampling counts, by using the counts directly from the real data.
#
# Now, if there's any effect from violation of the constancy assumption, it would affect both simulation approaches.
# For example, for the de novo simulation, the estimated variance/dispersion parameter would be inflated due to (hidden) variability in the normalization factor.
# I don't think it's much of a problem; just let the unaccounted variability in the spike-ins be absorbed into the biological variability (of the total cellular RNA).
# Then, just consider the simulation to be as near-realistic as possible, and accept that there will be some inaccuracy.
# The same applies for our simulation approach, where variability in the spike-ins is shifted to biological variability when spike-in volumes are assumed to be constant.
# This may be excusable if the biological variability is a lot higher than the spike-in variability, based on our experiments.
#
# Hopefully, reassignment of the hidden spike-in variablity to the biological variability in the simulation design shouldn't have any effect.
# However, if it does, then we should also see an effect when we introduce our additional variability.
# This excludes pathological scenarios where the hidden variability saturates the system such that any additional variability has no effect.
# We can probably protect against this by showing that the effect increases with increasing variance.
#
# In fact, we provide some decent protection by rescaling the means to eliminate the putative variability.
# To wit, we calculate variance in the log-transformed spike-in library size across all wells.
# This contains the contribution from spike-in variance, as well as that from everything else, e.g., cell-specific RT biases.
# We subtract the former from the total and recalculate the library sizes with the new variability, using a quantile-quantile approach under normality.
# The new library sizes are treated as the values that would have been obtained with no spike-in variability.
# New variability is added to obtain new library sizes, and new counts are computed with quantile-quantile scaling as described above.
# Results using these new counts are compared to the original data to determine the effect of spike-in variability on the results.

###########################################################################
# Try to sell this in terms of a perspective of a realistic simulation, rather than as an alteration of real data.
# In particular, say something like we don't have to change the cellular counts because we're just interested in spike-in variability.
# Then, say that we simulate the spike-in counts (and spike-in variability) using the procedure above.
# This procedure gives us realistic parameter estimates for the mean, and already has realistic technical noise included!
# At all times, it should be clear that we're doing a simulation, and it should protect us from criticisms about the real data being different from our experiments.

###########################################################################

spikeParam <- function(spike.counts) 
# This computes the parameters of the spike-in counts, including the
# fitted values and the dispersions used for rescaling later.
{
	spike.counts <- as.matrix(spike.counts)
	totals <- colSums(spike.counts)
	dispersion <- edgeR::estimateDisp(edgeR::DGEList(spike.counts, lib.size=totals), 
		design=cbind(rep(1, length(totals))), tagwise=FALSE)$trended.dispersion # Fits are usually good.
	central.mu <- exp(edgeR::mglmOneGroup(spike.counts, dispersion=dispersion, offset=log(totals)))
	cur.mus <- outer(central.mu, totals)
	return(list(counts=spike.counts, lib.size=totals, mean=central.mu, fitted=cur.mus, dispersion=dispersion))
}

resampleSpikes <- function(param, var.log)
# Removes the putative spike-in variability from the fitted values of the
# original data, and then randomly adds it back in (assuming variance refers to
# the log-volume with base 'e'). Counts are rescaled such that quantiles are
# preserved between the original and new fitted values, given the fitted
# parameters to the spike-ins.
{
	# Removing the effect of spike-in variability from the original library sizes.
	log.lib <- log(param$lib.size)
	cur.var <- var(log.lib)
	if (cur.var < var.log) { warning("variance in original library sizes is lower than specified") }
	new.log.lib <- log.lib * sqrt(max(0, 1 - var.log/cur.var))
	new.log.lib <- new.log.lib - mean(new.log.lib) + mean(log.lib)
	new.lib.size <- exp(new.log.lib)

	# Adding spike-in variability back in, randomly.
	fc <- exp(rnorm(length(new.lib.size), 0, sd=sqrt(var.log)))
	new.fitted <- outer(param$mean, fc*new.lib.size)
	new.counts <- edgeR::q2qnbinom(param$counts, input=param$fitted, output=new.fitted, dispersion=param$dispersion)
	new.counts <- round(pmax(new.counts, 0)) # Making it reasonably count-like.
	dim(new.counts) <- dim(param$counts)
	return(new.counts)
}

