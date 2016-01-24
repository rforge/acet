AtCtEp_mcmc <-
function(AtCtEp, iter_num = 10000, sd = 0.1, burnin =1000)
{

if(class(AtCtEp)!='AtCtEp_model')
{
	stop('The first parameter must be an object obtained from the AtCtEp function.')
}

order <- 3

B_des_a_m <- splineDesign(AtCtEp$knot_a, x=AtCtEp$T_m, ord=order)
B_des_a_d <- splineDesign(AtCtEp$knot_a, x=AtCtEp$T_d, ord=order)
B_des_c_m <- splineDesign(AtCtEp$knot_c, x=AtCtEp$T_m, ord=order)
B_des_c_d <- splineDesign(AtCtEp$knot_c, x=AtCtEp$T_d, ord=order)

result <- mcmc_epsp(AtCtEp$pheno_m, AtCtEp$pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, AtCtEp$var, AtCtEp$var_b_a, AtCtEp$var_b_c, AtCtEp$D_a, AtCtEp$D_c, iter_num, burnin, sd)

AtCtEp_mc_mod <- list(beta_a_mc=result$beta_a_mc, beta_c_mc=result$beta_c_mc, cov_a_mc = result$cov_a, cov_c_mc = result$cov_c, knots_a=AtCtEp$knot_a, knots_c=AtCtEp$knot_c)

class(AtCtEp_mc_mod) <- 'AtCtEp_mc_model'

return(AtCtEp_mc_mod)
}
