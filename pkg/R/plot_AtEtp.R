plot_AtEtp <- function(AtEtp_mcmc, data_m, data_d)
{
	if(class(AtEtp_mcmc)!='AtEtp_mc_model')
	{
		stop('The first parameter must be an object obtained from the AtEtp_mcmc function.')
	}
	
	model_cur <- AtEtp_mcmc

	pheno_m <- c(t(data_m[,1:2]))
	pheno_d <- c(t(data_d[,1:2]))
	T_m <- rep(data_m[,3], each=2)
	T_d <- rep(data_d[,3], each=2)

	order <- 3
	x <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=500)
	bb <- splineDesign(model_cur$knots_e, x = x, ord=order, outer.ok = TRUE)
	points_e <- exp(bb%*%model_cur$beta_e_mc)
	bb <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	points_a <- exp(bb%*%model_cur$beta_a_mc)

	plot(range(x), c(0,max(points_e, points_a)+1), type = "n", xlab = "Age", ylab = "Variance",main =  "Variance of the A and E component over age")
	lines(x, points_e, col = "pink", lwd = 2)
	fisher_e <- model_cur$cov_e_mc
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt((t(bb[i,])%*%fisher_e%*%bb[i,]))
		lower[i] <- sum(bb[i,]*model_cur$beta_e_mc) - 1.96*sd[i]
		upper[i] <- sum(bb[i,]*model_cur$beta_e_mc) + 1.96*sd[i]
	}

	lines(x, exp(lower), col = "yellow" ,lty = 2 , lwd = 0.6)
	lines(x, exp(upper), col = "yellow" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(exp(upper), rev(exp(lower))),col='grey',border = NA, lty=3, density=20)


	bb <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	lines(x, points_a, col = "red", lwd = 2)

	fisher_a <- model_cur$cov_a_mc
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt((t(bb[i,])%*%fisher_a%*%bb[i,]))
		lower[i] <- sum(bb[i,]*model_cur$beta_a_mc) - 1.96*sd[i]
		upper[i] <- sum(bb[i,]*model_cur$beta_a_mc) + 1.96*sd[i]
	}

	lines(x, exp(lower), col = "orange" ,lty = 2 , lwd = 0.6)
	lines(x, exp(upper), col = "orange" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(exp(upper), rev(exp(lower))),col='grey',border = NA, lty=3, density=20)


	legend(x[1], max(points_e, points_a), c('Additive genetic component','Unique environmental component'), col = c('red','pink'), lty=c(1,1), lwd=c(2,2))

}