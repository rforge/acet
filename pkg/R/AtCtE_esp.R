loglik_AtCtE_esp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d)
{
var <- param[1]
beta_c <- param[2:(1+ncol(B_des_c_m))]
beta_a <- param[(2+ncol(B_des_c_m)):(1+ncol(B_des_c_m)+ncol(B_des_a_m))]

nll <- .Call('loglik_AtCtE_esp_c', var, beta_c, beta_a, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d)

return(nll)
}

gr_AtCtE_esp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d)
{
var <- param[1]
beta_c <- param[2:(1+ncol(B_des_c_m))]
beta_a <- param[(2+ncol(B_des_c_m)):(1+ncol(B_des_c_m)+ncol(B_des_a_m))]

num_beta_c <- ncol(B_des_c_m)
num_beta_a <- ncol(B_des_a_m)

d <- .Call('gr_AtCtE_esp_c', var, beta_c, beta_a, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d)

return(d[1:(1+num_beta_c+num_beta_a)]) 
}
