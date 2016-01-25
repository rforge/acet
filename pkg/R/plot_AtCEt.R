plot_AtCEt <- function(AtCEt, data_m, data_d)
{
	if(class(AtCEt)!='AtCEt_model')
	{
		stop('The first parameter must be an object obtained from the AtCEt function.')
	}
	
	model_cur <- AtCEt
	
	pheno_m <- c(t(data_m[,1:2]))
	pheno_d <- c(t(data_d[,1:2]))
	T_m <- rep(data_m[,3], each=2)
	T_d <- rep(data_d[,3], each=2)

	order <- 3
	x <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=500)
	bb_e <- splineDesign(model_cur$knots_e, x = x, ord=order, outer.ok = TRUE)
	points_e <- exp(bb_e%*%model_cur$beta_e)
	bb_a <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	points_a <- exp(bb_a%*%model_cur$beta_a)

	l_a <- length(model_cur$beta_a)
	l_e <- length(model_cur$beta_e)

	plot(range(x), c(0,max(points_e, points_a, model_cur$var)+1), type = "n", xlab = "Age", ylab = "Variance",main =  "Variance curves of the A, C and E components")
	lines(x, points_a, col = "red", lwd = 2)
	
	fisher <- solve(model_cur$hessian[2:(1+l_a+l_e),2:(1+l_a+l_e)])
	
	fisher_a <- fisher[1:l_a,1:l_a]
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt(exp(2*sum(bb_a[i,]*model_cur$beta_a))*(t(bb_a[i,])%*%fisher_a%*%bb_a[i,]))
		lower[i] <- exp(sum(bb_a[i,]*model_cur$beta_a)) - 1.96*sd[i]
		upper[i] <- exp(sum(bb_a[i,]*model_cur$beta_a)) + 1.96*sd[i]
	}

	lines(x, lower, col = "orange" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "orange" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)


	#bb <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	lines(x, points_e, col = "pink", lwd = 2)
	
	fisher_e <- fisher[(1+l_a):(l_e+l_a),(1+l_a):(l_a+l_e)]
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt(exp(2*sum(bb_e[i,]*model_cur$beta_e))*(t(bb_e[i,])%*%fisher_e%*%bb_e[i,]))
		lower[i] <- exp(sum(bb_e[i,]*model_cur$beta_e)) - 1.96*sd[i]
		upper[i] <- exp(sum(bb_e[i,]*model_cur$beta_e)) + 1.96*sd[i]
	}

	lines(x, lower, col = "yellow" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "yellow" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)
	
	points_c <- rep(model_cur$var, length(x))
	lower <- rep(model_cur$var-1.96/sqrt(model_cur$hessian[1]),length(x))
	upper <- rep(model_cur$var+1.96/sqrt(model_cur$hessian[1]),length(x))
	
	lines(x, points_c, col = "blue", lwd = 2)
	lines(x, lower, col = "green" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "green" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)

	legend(x[1], max(points_e, points_a, model_cur$var)+1, c('Additive genetic component','Common environmental component','Unique environmental component'), col = c('red','blue','pink'), lty=c(1,1,1), lwd=c(2,2,2))

}