gr_AtEt_epsp_g <- 
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d, var_b_a, var_b_e, D_a, D_e)
{
	beta_a <- param[1:ncol(D_a)]
	beta_e <- param[(ncol(D_a)+1):(ncol(D_a)+ncol(D_e))]
	
	num_m <- nrow(pheno_m)
	K_2m <- matrix(1,2,2)
	num_d <- nrow(pheno_d)
	K_2d <- matrix(0.5,2,2)
	K_1d <- matrix(c(0.5,0,0,0.5),2,2)
	
	n_a <- ncol(D_a)
	n_e <- ncol(D_e)
	b_a_m <- rep(0, n_a)
	b_e_m <- rep(0, n_e)

	for(i in 1:(0.5*num_m))
	{
		start <- 2*i - 1
		end <- start + 1
		BA_m <- exp(sum(beta_a*B_des_a_m[start,]))
		BE_m <- exp(sum(beta_e*B_des_e_m[start,]))
		V_m <- BA_m*K_2m+BE_m*diag(1,2)
		inv_V_m <- solve(V_m)
		inv_ph_m <- inv_V_m%*%pheno_m[start:end]
		ph_inv_m <- t(pheno_m[start:end])%*%inv_V_m

		temp_a <- BA_m*(sum(diag(inv_V_m%*%K_2m))-ph_inv_m%*%K_2m%*%inv_ph_m)
		
		for(j in 1:n_a)
		{
			b_a_m[j] <- b_a_m[j] + B_des_a_m[start,j]*temp_a 
		}
		
		temp_e <- BE_m*(sum(diag(inv_V_m))-ph_inv_m%*%inv_ph_m)

		for(j in 1:n_e)
		{
			b_e_m[j] <- b_e_m[j] + B_des_e_m[start,j]*temp_e
		}
	}
	
	b_a_d <- rep(0, n_a)
	b_e_d <- rep(0, n_e)
	
	for(i in 1:(0.5*num_d))
	{
		start <- 2*i - 1
		end <- start + 1
		BA_d <- exp(sum(beta_a*B_des_a_d[start,]))
		BE_d <- exp(sum(beta_e*B_des_e_d[start,]))
		V_d <- BA_d*(K_1d+K_2d)+BE_d*diag(1,2)
		inv_V_d <- solve(V_d)
		inv_ph_d <- inv_V_d%*%pheno_d[start:end]
		ph_inv_d <- t(pheno_d[start:end])%*%inv_V_d

		temp_a <- BA_d*(sum(diag(inv_V_d%*%(K_1d+K_2d)))-ph_inv_d%*%(K_1d+K_2d)%*%inv_ph_d)
		for(j in 1:ncol(D_a))
		{
			b_a_d[j] <- b_a_d[j] + B_des_a_d[start,j]*temp_a
		}
		temp_e <- BE_d*(sum(diag(inv_V_d))-ph_inv_d%*%inv_ph_d)
		for(j in 1:ncol(D_e))
		{
			b_e_d[j] <- b_e_d[j] + B_des_e_d[start,j]*temp_e 
		}

	}

	d_b_a <- b_a_m + b_a_d + ((D_a+t(D_a))%*%beta_a)[,1]/var_b_a
	d_b_e <- b_e_m + b_e_d + ((D_e+t(D_e))%*%beta_e)[,1]/var_b_e
		
	return(c(d_b_a[1:n_a], d_b_e[1:n_e]))
}
