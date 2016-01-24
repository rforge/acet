loglik_AtEt_epsp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_e_m, B_des_e_d, beta_e, D_e)
{
var_b_a <- param[1]
var_b_e <- param[2]

num_m <- nrow(pheno_m)
K_2m <- matrix(1,2,2)
num_d <- nrow(pheno_d)
K_2d <- matrix(0.5,2,2)
K_1d <- matrix(c(0.5,0,0,0.5),2,2)

YSY_m <- 0
	YSY_d <- 0
	D_m <- 0
	D_d <- 0
	gsd_a2m2a <- c()
	gsd_e2mm2a <- c()
	gsd_emme <- c()
	gsd_adda <- c()
	gsd_adde <- c()
	gsd_edde <- c()
	gsd_a2dd2a <- c()
	gsd_a2dd2e <- c()
	gsd_2d2 <- c()

	for(i in 1:(0.5*num_m))
	{
		start <- 2*i - 1
		end <- start + 1
		BA_m <- exp(sum(beta_a*B_des_a_m[start,]))
		BE_m <- exp(sum(beta_e*B_des_e_m[start,]))
		V_m <- BA_m*K_2m+BE_m*diag(1,2)
		inv_V_m <- solve(V_m)
		YSY_m <- YSY_m + t(pheno_m)[start:end]%*%inv_V_m%*%pheno_m[start:end]
		D_m <- D_m + log(V_m[1,1]^2-V_m[1,2]^2)
		temp <- sum(inv_V_m%*%K_2m%*%inv_V_m)
		gsd_a2m2a <- c(gsd_a2m2a, (BA_m^2)*temp)
		temp_mm <- inv_V_m%*%inv_V_m
		gsd_emme <- c(gsd_emme, (BE_m^2)*diag(temp_mm))
		gsd_e2mm2a <- c(gsd_e2mm2a, (BE_m*BA_m)*sum(temp_mm))
	}

	for(i in 1:(0.5*num_d))
	{
		start <- 2*i - 1
		end <- start + 1
		BA_d <- exp(sum(beta_a*B_des_a_d[start,]))
		BE_d <- exp(sum(beta_e*B_des_e_d[start,]))
		BA2_d <- BA_d^2
		BAE_d <- BA_d*BE_d
		BE2_d <- BE_d^2
		V_d <- BA_d*(K_1d+K_2d)+BE_d*diag(1,2)
		inv_V_d <- solve(V_d)
		YSY_d <- YSY_d + t(pheno_d)[start:end]%*%inv_V_d%*%pheno_d[start:end]
		D_d <- D_d + log(V_d[1,1]^2-V_d[1,2]^2)
		temp_dd <- inv_V_d%*%inv_V_d
		gsd_adda <- c(gsd_adda, BA2_d*diag(temp_dd))
		gsd_adde <- c(gsd_adde, BAE_d*diag(temp_dd))
		gsd_edde <- c(gsd_edde, BE2_d*diag(temp_dd))
		gsd_a2dd2a <- c(gsd_a2dd2a, BA2_d*0.5*sum(temp_dd))
		gsd_a2dd2e <- c(gsd_a2dd2e, BAE_d*0.5*sum(temp_dd))
		temp_si2i <- sum(inv_V_d%*%K_2d%*%inv_V_d)
		gsd_2d2 <- c(gsd_2d2, BA2_d*0.5*temp_si2i)
	}
	
	r_d_a <- nrow(D_a)
	r_d_e <- nrow(D_e)
	gsd_max <- matrix(NA, r_d_a+r_d_e, r_d_a+r_d_e)
	gsd_max[1:r_d_a,1:r_d_a] <- diag(diag(D_a+t(D_a)))/var_b_a + t(apply(t(B_des_a_m[seq(from=1, to=num_m, by=2),]),1,function(x) x*gsd_a2m2a))%*%B_des_a_m[seq(from=1, to=num_m, by=2),] + 0.25*t(apply(t(B_des_a_d),1,function(x) x*gsd_adda))%*%B_des_a_d + t(apply(t(B_des_a_d[seq(from=1, to=num_d, by=2),]),1,function(x) x*(gsd_a2dd2a+gsd_2d2)))%*%B_des_a_d[seq(from=1, to=num_d, by=2),]
	gsd_max[(r_d_a+1):(r_d_a+r_d_e), 1:r_d_a] <- t(apply(t(B_des_e_m[seq(from=1, to=num_m, by=2),]),1,function(x) x*gsd_e2mm2a))%*%B_des_a_m[seq(from=1, to=num_m, by=2),] + 0.5*t(apply(t(B_des_e_d),1,function(x) x*gsd_adde))%*%B_des_a_d + t(apply(t(B_des_e_d[seq(from=1, to=num_d, by=2),]),1,function(x) x*gsd_a2dd2e))%*%B_des_a_d[seq(from=1, to=num_d, by=2),]
	gsd_max[1:r_d_a,(r_d_a+1):(r_d_a+r_d_e)] <- t(gsd_max[(r_d_a+1):(r_d_a+r_d_e), 1:r_d_a])
	gsd_max[(r_d_a+1):(r_d_a+r_d_e),(r_d_a+1):(r_d_a+r_d_e)] <- diag(diag(D_e+t(D_e)))/var_b_e + t(apply(t(B_des_e_m),1,function(x) x*gsd_emme))%*%B_des_e_m + t(apply(t(B_des_e_d),1,function(x) x*gsd_edde))%*%B_des_e_d 

	gsd <- determinant(0.5*(gsd_max), logarithm = TRUE)

	# gsd <- determinant(0.5*(diag(diag(D+t(D)))/var_b+t(apply(t(B_des_m[seq(from=1, to=num_m, by=2),]),1,function(x) x*gsd_m))%*%B_des_m[seq(from=1, to=num_m, by=2),]+t(apply(t(B_des_d[seq(from=1, to=num_d, by=2),]),1, function(x) x*gsd_d))%*%B_des_d[seq(from=1, to=num_d, by=2),]), logarithm = TRUE)
	
	nll <- D_m + YSY_m + D_d + YSY_d + (length(beta_a)-2)*log(var_b_a) + (length(beta_e)-2)*log(var_b_e) + (t(beta_a)%*%D_a%*%beta_a)/var_b_a + (t(beta_e)%*%D_e%*%beta_e)/var_b_e + gsd$modulus
	return(nll)
}
