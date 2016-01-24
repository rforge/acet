loglik_AtCEt_esp <-
  function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d)
  {
    var_c <- param[1]
    beta_a <- param[2:(1+ncol(B_des_a_m))]
    beta_e <- param[(2+ncol(B_des_a_m)):(1+ncol(B_des_a_m)+ncol(B_des_e_m))]
    
    nll <- .Call('loglik_AtCEt_esp_c', var_c, beta_a, beta_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d)
    
    return(nll)
  }

gr_AtCEt_esp <- function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d)
{
  var_c <- param[1]
  beta_a <- param[2:(1+ncol(B_des_a_m))]
  beta_e <- param[(2+ncol(B_des_a_m)):(1+ncol(B_des_a_m)+ncol(B_des_e_m))]
  
  num_beta_a <- ncol(B_des_a_m)
  num_beta_e <- ncol(B_des_e_m)
  
  d <- .Call('gr_AtCEt_esp_c', var_c, beta_a, beta_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d)

  return(d) 
}
