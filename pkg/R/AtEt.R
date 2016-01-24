AtEt <-
function(data_m, data_d, knot_a=5, knot_e=5)
{

pheno_m <- c(t(data_m[,1:2]))
pheno_d <- c(t(data_d[,1:2]))
T_m <- rep(data_m[,3], each=2)
T_d <- rep(data_d[,3], each=2)

order <- 3
#knot <- 8

knots_a <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_a)
interval_a <- knots_a[2] - knots_a[1]
knots_a <- c(c(min(T_m, T_d)-interval_a*2,min(T_m, T_d)-interval_a), knots_a)
knots_a <- c(knots_a, c(max(T_m, T_d)+interval_a,max(T_m, T_d)+interval_a*2))
B_des_a_m <- splineDesign(knots_a, x=T_m, ord=order)
B_des_a_d <- splineDesign(knots_a, x=T_d, ord=order)

knots_e <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_e)
interval_e <- knots_e[2] - knots_e[1]
knots_e <- c(c(min(T_m, T_d)-interval_e*2,min(T_m, T_d)-interval_e), knots_e)
knots_e <- c(knots_e, c(max(T_m, T_d)+interval_e,max(T_m, T_d)+interval_e*2))
B_des_e_m <- splineDesign(knots_e, x=T_m, ord=order)	
B_des_e_d <- splineDesign(knots_e, x=T_d, ord=order)


n_e <- ncol(B_des_e_m)
n_a <- ncol(B_des_a_m)

result <- optim(c(0,rep(0, n_a),rep(0, n_e)), loglik_AtCEt_esp, gr_AtCEt_esp, pheno_m = matrix(pheno_m), pheno_d = matrix(pheno_d), B_des_a_m = B_des_a_m, B_des_a_d = B_des_a_d, B_des_e_m = B_des_e_m, B_des_e_d = B_des_e_d, lower = c(rep(0,1),rep(-5, n_a),rep(-5, n_e)), upper = c(rep(0,1),rep(5, n_a),rep(5, n_e)), method = "L-BFGS-B", hessian = TRUE)

AtEt_model <- list(n_beta_a=n_a, n_beta_e=n_e, beta_a=result$par[2:(1+n_a)], beta_e=result$par[(2+n_a):(1+n_a+n_e)], hessian=result$hessian[2:(1+n_a+n_e),2:(1+n_a+n_e)], con=result$convergence, lik=result$value, knots_a =knots_a, knots_e = knots_e )

class(AtEt_model) <- 'AtEt_model'

return(AtEt_model)

}
