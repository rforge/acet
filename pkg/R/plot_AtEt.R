plot_AtEt <- function(AtEt, data_m, data_d)
{
	if(class(AtEt)!='AtEt_model')
	{
		stop('The first parameter must be an object obtained from the AtEt function.')
	}

	model_cur <- AtEt

	pheno_m <- c(t(data_m[,1:2]))
	pheno_d <- c(t(data_d[,1:2]))
	T_m <- rep(data_m[,3], each=2)
	T_d <- rep(data_d[,3], each=2)

	order <- 3
	x <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=500)
	bb <- splineDesign(model_cur$knots_e, x = x, ord=order, outer.ok = TRUE)
	points_e <- exp(bb%*%model_cur$beta_e)
	bb <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	points_a <- exp(bb%*%model_cur$beta_a)

	l_e <- length(model_cur$beta_e)
	l_a <- length(model_cur$beta_a)
	
	fisher <- solve(model_cur$hessian[1:(l_a+l_e),1:(l_a+l_e)])

	plot(range(x), c(0,max(points_e, points_a)+1), type = "n", xlab = "Age", ylab = "Variance",main =  "Variance of the A and E component over age")
	lines(x, points_e, col = "red", lwd = 2)

	fisher_a <- fisher[1:l_a,1:l_a]
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt(exp(2*sum(bb[i,]*model_cur$beta_a))*(t(bb[i,])%*%fisher_a%*%bb[i,]))
		lower[i] <- exp(sum(bb[i,]*model_cur$beta_a)) - 1.96*sd[i]
		upper[i] <- exp(sum(bb[i,]*model_cur$beta_a)) + 1.96*sd[i]
	}

	lines(x, lower, col = "orange" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "orange" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)


	bb <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	lines(x, points_a, col = "pink", lwd = 2)
	
	fisher_e <- fisher[(1+l_a):(l_a+l_e),(1+l_a):(l_a+l_e)]
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt(exp(2*sum(bb[i,]*model_cur$beta_e))*(t(bb[i,])%*%fisher_e%*%bb[i,]))
		lower[i] <- exp(sum(bb[i,]*model_cur$beta_e)) - 1.96*sd[i]
		upper[i] <- exp(sum(bb[i,]*model_cur$beta_e)) + 1.96*sd[i]
	}

	lines(x, lower, col = "yellow" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "yellow" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)


	legend(x[1], max(points_e, points_a)+1, c('Additive genetic component','Unique environmental component'), col = c('red','pink'), lty=c(1,1), lwd=c(2,2))

}