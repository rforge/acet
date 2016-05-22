acetp_mcmc <- function(acetp, iter_num = 10000, sd = 0.1, burnin =1000)
{
	if(!(class(acetp) %in% c('AtCtEp_model','AtEtp_model','AtCtEtp_model')))
	{
		stop('The first parameter must be an acetp object.')
	}

	
	if(burnin >= iter_num)
	{
		stop('The number of burnins must be smaller than the number of MCMC iterations.')
	}


	if(class(acetp)=='AtCtEtp_model')
	{
		res <- AtCtEtp_mcmc(acetp, iter_num, sd, burnin)
		return(res)
	}

	if(class(acetp)=='AtEtp_model')
	{
		res <- AtEtp_mcmc(acetp, iter_num, sd, burnin)
		return(res)
	}
	
}




AtCtEtp_mcmc <-
function(AtCtEtp, iter_num = 10000, sd = 0.1, burnin =1000)
{

if(class(AtCtEtp)!='AtCtEtp_model')
{
	stop('The first parameter must be an object obtained from the AtCtEtp function.')
}

T_m <- AtCtEtp$T_m
num_m <- length(T_m)
T_d <- AtCtEtp$T_d
num_d <- length(T_d)

t_int <- max(c(T_m,T_d))-min(c(T_m,T_d))
l_m_1 <- (max(c(T_m,T_d))-T_m)/t_int
l_m_2 <- (T_m-min(c(T_m,T_d)))/t_int
l_d_1 <- (max(c(T_m,T_d))-T_d)/t_int
l_d_2 <- (T_d-min(c(T_m,T_d)))/t_int

order <- 3
if(length(AtCtEtp$beta_a)>2)
{
	B_des_a_m <- splineDesign(AtCtEtp$knot_a, x=AtCtEtp$T_m, ord=order)
	B_des_a_d <- splineDesign(AtCtEtp$knot_a, x=AtCtEtp$T_d, ord=order)
}else{
	if(length(AtCtEtp$beta_a)==2)
	{
		B_des_a_m <- matrix(NA, num_m, 2)
		B_des_a_m[,1] <- l_m_1
		B_des_a_m[,2] <- l_m_2
		B_des_a_d <- matrix(NA, num_d, 2)
		B_des_a_d[,1] <- l_d_1
		B_des_a_d[,2] <- l_d_2
	}else{
		B_des_a_m <- matrix(1, num_m, 1)
		B_des_a_d <- matrix(1, num_d, 1)
	}
}
if(length(AtCtEtp$beta_c)>2)
{
	B_des_c_m <- splineDesign(AtCtEtp$knot_c, x=AtCtEtp$T_m, ord=order)
	B_des_c_d <- splineDesign(AtCtEtp$knot_c, x=AtCtEtp$T_d, ord=order)
}else{
	if(length(AtCtEtp$beta_c)==2)
	{
		B_des_c_m <- matrix(NA, num_m, 2)
		B_des_c_m[,1] <- l_m_1
		B_des_c_m[,2] <- l_m_2
		B_des_c_d <- matrix(NA, num_d, 2)
		B_des_c_d[,1] <- l_d_1
		B_des_c_d[,2] <- l_d_2
	}else{
		B_des_c_m <- matrix(1, num_m, 1)
		B_des_c_d <- matrix(1, num_d, 1)
	}
}
if(length(AtCtEtp$beta_e)>2)
{
	B_des_e_m <- splineDesign(AtCtEtp$knot_e, x=AtCtEtp$T_m, ord=order)
	B_des_e_d <- splineDesign(AtCtEtp$knot_e, x=AtCtEtp$T_d, ord=order)
}else{
	if(length(AtCtEtp$beta_e)==2)
	{
		B_des_e_m <- matrix(NA, num_m, 2)
		B_des_e_m[,1] <- l_m_1
		B_des_e_m[,2] <- l_m_2
		B_des_e_d <- matrix(NA, num_d, 2)
		B_des_e_d[,1] <- l_d_1
		B_des_e_d[,2] <- l_d_2
	}else{
		B_des_e_m <- matrix(1, num_m, 1)
		B_des_e_d <- matrix(1, num_d, 1)
	}
}

result <- mcmc_epsp_AtCtEt(AtCtEtp$pheno_m, AtCtEtp$pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, AtCtEtp$var_b_a, AtCtEtp$var_b_c, AtCtEtp$var_b_e, AtCtEtp$D_a, AtCtEtp$D_c, AtCtEtp$D_e, iter_num, burnin, sd)

AtCtEtp_mc_mod <- list(beta_a_mc=result$beta_a_mc, beta_c_mc=result$beta_c_mc, beta_e_mc=result$beta_e_mc, cov_mc = result$cov, knots_a=AtCtEtp$knot_a, knots_c=AtCtEtp$knot_c, knots_e=AtCtEtp$knot_e, min_t = min(AtCtEtp$T_m, AtCtEtp$T_d), max_t = max(AtCtEtp$T_m, AtCtEtp$T_d))

class(AtCtEtp_mc_mod) <- 'AtCtEtp_mc_model'

return(AtCtEtp_mc_mod)
}
