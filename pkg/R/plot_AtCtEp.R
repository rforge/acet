plot_AtCtEp <- function(AtCtEp_mcmc, data_m, data_d)
{
	if(class(AtCtEp_mcmc)!='AtCtEp_mc_model')
	{
		stop('The first parameter must be an object obtained from the AtCtEp_mcmc function.')
	}
	
	model_cur <- AtCtEp_mcmc

	pheno_m <- c(t(data_m[,1:2]))
	pheno_d <- c(t(data_d[,1:2]))
	T_m <- rep(data_m[,3], each=2)
	T_d <- rep(data_d[,3], each=2)

	order <- 3
	x <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=500)
	bb <- splineDesign(model_cur$knots_c, x = x, ord=order, outer.ok = TRUE)
	points_c <- exp(bb%*%model_cur$beta_c_mc)
	bb <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	points_a <- exp(bb%*%model_cur$beta_a_mc)

	plot(range(x), c(0,max(points_c, points_a)+1), type = "n", xlab = "Age", ylab = "Variance",main =  "Variance of the A and C component over age")
	lines(x, points_c, col = "blue", lwd = 2)
	fisher_c <- model_cur$cov_c_mc
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt((t(bb[i,])%*%fisher_c%*%bb[i,]))
		lower[i] <- sum(bb[i,]*model_cur$beta_c_mc) - 1.96*sd[i]
		upper[i] <- sum(bb[i,]*model_cur$beta_c_mc) + 1.96*sd[i]
	}

	lines(x, exp(lower), col = "green" ,lty = 2 , lwd = 0.6)
	lines(x, exp(upper), col = "green" ,lty = 2 , lwd = 0.6)
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


	legend(x[1], max(points_c, points_a)+1, c('Additive genetic component','Common environmental component'), col = c('red','blue'), lty=c(1,1), lwd=c(2,2))

}