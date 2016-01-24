AtEtp_mcmc <-
function(AtEtp, iter_num = 10000, sd = 0.1, burnin =1000)
{

if(class(AtEtp)!='AtEtp_model')
{
	stop('The first parameter must be an object obtained from the AtEtp function.')
}

model_cur <- AtEtp

order <- 3

B_des_a_m <- splineDesign(model_cur$knot_a, x=T_m, ord=order)
B_des_a_d <- splineDesign(model_cur$knot_a, x=T_d, ord=order)
B_des_e_m <- splineDesign(model_cur$knot_e, x=T_m, ord=order)
B_des_e_d <- splineDesign(model_cur$knot_e, x=T_d, ord=order)

result <- mcmc_epsp_AtEtp(model_cur$pheno_m, model_cur$pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d, model_cur$var_b_a, model_cur$var_b_e, model_cur$D_a, model_cur$D_e, iter_num, burnin, sd)

AtEtp_mc_mod <- list(beta_a_mc=result$beta_a_mc, beta_e_mc=result$beta_e_mc, cov_a_mc = result$cov_a, cov_e_mc = result$cov_e, knots_a=model_cur$knots_a, knots_e=model_cur$knots_e)

class(AtEtp_mc_mod) <- 'AtEtp_mc_model'

return(AtEtp_mc_mod)
}
