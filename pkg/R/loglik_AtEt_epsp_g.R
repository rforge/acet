loglik_AtEt_epsp_g <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d, var_b_a, var_b_e, D_a, D_e)
{
beta_a <- param[1:ncol(D_a)]
beta_e <- param[(ncol(D_a)+1):(ncol(D_a)+ncol(D_e))]

num_m <- nrow(pheno_m)
K_2m <- matrix(1,2,2)
num_d <- nrow(pheno_d)
K_2d <- matrix(0.5,2,2)
K_1d <- matrix(c(0.5,0,0,0.5),2,2)

YSY_m <- 0
YSY_d <- 0
D_m <- 0
D_d <- 0

for(i in 1:(0.5*num_m))
	{
		start <- 2*i - 1
		end <- start + 1
		V_m <- exp(sum(beta_a*B_des_a_m[start,]))*K_2m+exp(sum(beta_e*B_des_e_m[start,]))*diag(1,2)
		inv_V_m <- solve(V_m)
		YSY_m <- YSY_m + t(pheno_m)[start:end]%*%inv_V_m%*%pheno_m[start:end]
		D_m <- D_m + log(V_m[1,1]^2-V_m[1,2]^2)
	}

	for(i in 1:(0.5*num_d))
	{
		start <- 2*i - 1
		end <- start + 1
		V_d <- exp(sum(beta_a*B_des_a_d[start,]))*(K_1d+K_2d)+exp(sum(beta_e*B_des_e_d[start,]))*diag(1,2)
		inv_V_d <- solve(V_d)
		YSY_d <- YSY_d + t(pheno_d)[start:end]%*%inv_V_d%*%pheno_d[start:end]
		D_d <- D_d + log(V_d[1,1]^2-V_d[1,2]^2)
	}
	
	nll <- D_m + YSY_m + D_d + YSY_d + (t(beta_a)%*%D_a%*%beta_a)/var_b_a + (t(beta_e)%*%D_e%*%beta_e)/var_b_e
	return(nll)
}
