AtCtE <-
function(data_m, data_d, knot_a=5, knot_c=5)
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

knots_c <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_c)
interval_c <- knots_c[2] - knots_c[1]
knots_c <- c(c(min(T_m, T_d)-interval_c*2,min(T_m, T_d)-interval_c), knots_c)
knots_c <- c(knots_c, c(max(T_m, T_d)+interval_c,max(T_m, T_d)+interval_c*2))
B_des_c_m <- splineDesign(knots_c, x=T_m, ord=order)
B_des_c_d <- splineDesign(knots_c, x=T_d, ord=order)

n_c <- ncol(B_des_c_m)
n_a <- ncol(B_des_a_m)

result <- optim(c(1,rep(0, n_c),rep(0, n_a)), loglik_AtCtE_esp, gr_AtCtE_esp, pheno_m = matrix(pheno_m), pheno_d = matrix(pheno_d), B_des_a_m = B_des_a_m, B_des_a_d = B_des_a_d, B_des_c_m = B_des_c_m, B_des_c_d = B_des_c_d, lower = c(rep(0.00001,1),rep(-5, n_c),rep(-5, n_a)), upper = c(rep(100,1),rep(5, n_c),rep(5, n_a)), method = "L-BFGS-B", hessian = TRUE)

AtCtE_model <- list(n_beta_a=n_a, n_beta_c=n_c, var=result$par[1], beta_c=result$par[2:(1+n_c)], beta_a=result$par[(2+n_c):(1+n_c+n_a)], hessian=result$hessian, con=result$convergence, lik=result$value, knots_a =knots_a, knots_c = knots_c )

class(AtCtE_model) <- 'AtCtE_model'

return(AtCtE_model)

}
