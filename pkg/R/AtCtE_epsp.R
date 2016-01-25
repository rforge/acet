loglik_AtCtE_epsp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c)
{
var <- param[1]
var_b_a <- param[2]
var_b_c <- param[3]

nll <- .Call('loglik_AtCtE_epsp_c', var, var_b_a, var_b_c, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c)

return(nll)
}

loglik_AtCtE_epsp_g <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, var, var_b_a, var_b_c, D_a, D_c)
{
beta_a <- param[1:ncol(D_a)]
beta_c <- param[(ncol(D_a)+1):(ncol(D_a)+ncol(D_c))]

nll <- .Call('loglik_AtCtE_epsp_g_c', beta_a, beta_c, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, var, var_b_a, var_b_c, D_a, D_c)

return(nll)
}

#gr_AtCtE_epsp <-
#function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c)
#{
#var <- param[1]
#var_b_a <- param[2]
#var_b_c <- param[3]

#d <- .Call('gr_AtCtE_epsp_c', var, var_b_a, var_b_c, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c)

#return(d)
#}

gr_AtCtE_epsp_g <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, var, var_b_a, var_b_c, D_a, D_c)
{
beta_a <- param[1:ncol(D_a)]
beta_c <- param[(ncol(D_a)+1):(ncol(D_a)+ncol(D_c))]

d <- .Call('gr_AtCtE_epsp_g_c', beta_a, beta_c, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, var, var_b_a, var_b_c, D_a, D_c)

return(d)
}

gr_AtCtE_epsp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c)
{
var <- param[1]
var_b_a <- param[2]
var_b_c <- param[3]

num_m <- nrow(pheno_m)
K_2m <- matrix(1,2,2)
K_3m <- matrix(1,2,2)
num_d <- nrow(pheno_d)
K_3d <- matrix(1,2,2)
K_2d <- matrix(0.5,2,2)
K_1d <- matrix(c(0.5,0,0,0.5),2,2)

var_YSY_m <- 0
var_D_m <- 0
var_YSY_d <- 0
var_D_d <- 0
gsd_a2m2a <- c()
gsd_c2m2a <- c()
gsd_c2m2c <- c()
gsd_amm2ma <- c()
gsd_am2mma <- c()
gsd_cmm2ma <- c()
gsd_cm2mma <- c()
gsd_cmm2mc <- c()
gsd_cm2mmc <- c()
gsd_dd <- c()
gsd_2dd2 <- c()
gsd_3dd3 <- c()
gsd_2d2 <- c()
gsd_3d2 <- c()
gsd_3d3 <- c()
gsd_ddd <- c()
gsd_2ddd2 <- c()
gsd_2dd2d2 <- c()
gsd_2d2dd2 <- c()
gsd_3ddd3 <- c()
gsd_3dd2d3 <- c()
gsd_3d2dd3 <- c()
gsd_3dd3d3 <- c()
gsd_3d3dd3 <- c()

for(i in 1:(0.5*num_m))
{
start <- 2*i - 1
end <- start + 1
BA_m <- exp(sum(beta_a*B_des_a_m[start,]))
BC_m <- exp(sum(beta_c*B_des_c_m[start,]))
BA2_m <- BA_m^2
BC2_m <- BC_m^2
BAC_m <- BA_m*BC_m
V_m <- var*diag(1,2)+BA_m*K_2m+BC_m*K_3m
inv_V_m <- solve(V_m)
inv_ph_m <- inv_V_m%*%pheno_m[start:end]
ph_inv_m <- t(pheno_m[start:end])%*%inv_V_m
var_YSY_m <- var_YSY_m + ph_inv_m%*%inv_ph_m
var_D_m <- var_D_m + sum(diag(inv_V_m))
temp <- inv_V_m%*%K_2m%*%inv_V_m
gsd_a2m2a <- c(gsd_a2m2a, BA2_m*sum(temp))
gsd_c2m2a <- c(gsd_c2m2a, BAC_m*sum(temp))
gsd_c2m2c <- c(gsd_c2m2c, BC2_m*sum(temp))
temp_mm2m <- sum(inv_V_m%*%temp)
temp_m2mm <- sum(temp%*%inv_V_m)
gsd_amm2ma <- c(gsd_amm2ma, BA2_m*temp_mm2m)
gsd_am2mma <- c(gsd_am2mma, BA2_m*temp_m2mm)
gsd_cmm2ma <- c(gsd_cmm2ma, BAC_m*temp_mm2m)
gsd_cm2mma <- c(gsd_cm2mma, BAC_m*temp_m2mm)
gsd_cmm2mc <- c(gsd_cmm2mc, BC2_m*temp_mm2m)
gsd_cm2mmc <- c(gsd_cm2mmc, BC2_m*temp_m2mm)
}

for(i in 1:(0.5*num_d))
{
start <- 2*i - 1
end <- start + 1
BA_d <- exp(sum(beta_a*B_des_a_d[start,]))
BC_d <- exp(sum(beta_c*B_des_c_d[start,]))
BA2_d <- BA_d^2
BC2_d <- BC_d^2
BAC_d <- BA_d*BC_d
V_d <- var*diag(1,2)+BA_d*(K_1d+K_2d)+BC_d*K_3d
inv_V_d <- solve(V_d)
inv_ph_d <- inv_V_d%*%pheno_d[start:end]
ph_inv_d <- t(pheno_d[start:end])%*%inv_V_d
var_YSY_d <- var_YSY_d + ph_inv_d%*%inv_ph_d
var_D_d <- var_D_d + sum(diag(inv_V_d))
temp_dd <- inv_V_d%*%inv_V_d
gsd_dd <- c(gsd_dd, BA2_d*diag(temp_dd))
gsd_2dd2 <- c(gsd_2dd2, BA2_d*0.5*sum(temp_dd))
gsd_3dd3 <- c(gsd_3dd3, BAC_d*sum(temp_dd))
temp_i2i <- inv_V_d%*%K_2d%*%inv_V_d
gsd_2d2 <- c(gsd_2d2, BA2_d*0.5*sum(temp_i2i))
gsd_3d2 <- c(gsd_3d2, BAC_d*sum(temp_i2i))
temp_i3i <- 2*temp_i2i
gsd_3d3 <- c(gsd_3d3, BC2_d*sum(temp_i3i))
temp_ddd <- temp_dd%*%inv_V_d
gsd_ddd <- c(gsd_ddd, BA2_d*diag(temp_ddd))
gsd_2ddd2 <- c(gsd_2ddd2, BA2_d*0.5*sum(temp_ddd))
gsd_3ddd3 <- c(gsd_3ddd3, BAC_d*sum(temp_ddd))
temp_dd2d <- inv_V_d%*%temp_i2i
temp_d2dd <- temp_i2i%*%inv_V_d
temp_dd3d <- inv_V_d%*%temp_i3i
temp_d3dd <- temp_i3i%*%inv_V_d
gsd_2dd2d2 <- c(gsd_2dd2d2, BA2_d*0.5*sum(temp_dd2d))
gsd_2d2dd2 <- c(gsd_2d2dd2, BA2_d*0.5*sum(temp_d2dd))
gsd_3dd2d3 <- c(gsd_3dd2d3, BAC_d*sum(temp_dd2d))
gsd_3d2dd3 <- c(gsd_3d2dd3, BAC_d*sum(temp_d2dd))
gsd_3dd3d3 <- c(gsd_3dd3d3, BC2_d*sum(temp_dd3d))
gsd_3d3dd3 <- c(gsd_3d3dd3, BC2_d*sum(temp_d3dd))
}

r_d_a <- nrow(D_a)
r_d_c <- nrow(D_c)
gsd_max <- matrix(NA, r_d_a+r_d_c, r_d_a+r_d_c)
gsd_max[1:r_d_a,1:r_d_a] <- diag(diag(D_a+t(D_a)))/var_b_a + t(apply(t(B_des_a_m[seq(from=1, to=num_m, by=2),]),1,function(x) x*gsd_a2m2a))%*%B_des_a_m[seq(from=1, to=num_m, by=2),] + 0.25*t(apply(t(B_des_a_d),1,function(x) x*gsd_dd))%*%B_des_a_d + t(apply(t(B_des_a_d[seq(from=1, to=num_d, by=2),]),1,function(x) x*(gsd_2dd2+gsd_2d2)))%*%B_des_a_d[seq(from=1, to=num_d, by=2),]
gsd_max[(r_d_a+1):(r_d_a+r_d_c), 1:r_d_a] <- t(apply(t(B_des_c_m[seq(from=1, to=num_m, by=2),]),1,function(x) x*gsd_c2m2a))%*%B_des_a_m[seq(from=1, to=num_m, by=2),] + t(apply(t(B_des_c_d[seq(from=1, to=num_d, by=2),]),1,function(x) x*(0.5*gsd_3dd3+gsd_3d2)))%*%B_des_a_d[seq(from=1, to=num_d, by=2),]
gsd_max[1:r_d_a,(r_d_a+1):(r_d_a+r_d_c)] <- t(gsd_max[(r_d_a+1):(r_d_a+r_d_c), 1:r_d_a])
gsd_max[(r_d_a+1):(r_d_a+r_d_c),(r_d_a+1):(r_d_a+r_d_c)] <- diag(diag(D_c+t(D_c)))/var_b_c + t(apply(t(B_des_c_m[seq(from=1, to=num_m, by=2),]),1,function(x) x*gsd_c2m2c))%*%B_des_c_m[seq(from=1, to=num_m, by=2),] + t(apply(t(B_des_c_d[seq(from=1, to=num_d, by=2),]),1,function(x) x*gsd_3d3))%*%B_des_c_d[seq(from=1, to=num_d, by=2),]

inv_gsd_max <- solve(gsd_max)
gsd_b_a <- matrix(0, r_d_a+r_d_c, r_d_a+r_d_c)
gsd_b_a[1:r_d_a,1:r_d_a] <- diag(diag(D_a+t(D_a)))
gsd_b_c <- matrix(0, r_d_a+r_d_c, r_d_a+r_d_c)
gsd_b_c[(r_d_a+1):(r_d_a+r_d_c),(r_d_a+1):(r_d_a+r_d_c)] <- diag(diag(D_c+t(D_c)))
gsd_e <- matrix(0, r_d_a+r_d_c, r_d_a+r_d_c)
gsd_e[1:r_d_a,1:r_d_a] <- t(apply(t(B_des_a_m[seq(from=1, to=num_m, by=2),]),1,function(x) x*(gsd_amm2ma+gsd_am2mma)))%*%B_des_a_m[seq(from=1, to=num_m, by=2),] + 0.5*t(apply(t(B_des_a_d),1,function(x) x*gsd_ddd))%*%B_des_a_d + t(apply(t(B_des_a_d[seq(from=1, to=num_d, by=2),]),1,function(x) x*(2*gsd_2ddd2+gsd_2dd2d2+gsd_2d2dd2)))%*%B_des_a_d[seq(from=1, to=num_d, by=2),]
gsd_e[(r_d_a+1):(r_d_a+r_d_c), 1:r_d_a] <- t(apply(t(B_des_c_m[seq(from=1, to=num_m, by=2),]),1,function(x) x*(gsd_cmm2ma+gsd_cm2mma)))%*%B_des_a_m[seq(from=1, to=num_m, by=2),] + t(apply(t(B_des_c_d[seq(from=1, to=num_d, by=2),]),1,function(x) x*(gsd_3ddd3+gsd_3dd2d3+gsd_3d2dd3)))%*%B_des_a_d[seq(from=1, to=num_d, by=2),]
gsd_e[1:r_d_a,(r_d_a+1):(r_d_a+r_d_c)] <- t(gsd_e[(r_d_a+1):(r_d_a+r_d_c), 1:r_d_a])
gsd_e[(r_d_a+1):(r_d_a+r_d_c),(r_d_a+1):(r_d_a+r_d_c)] <- t(apply(t(B_des_c_m[seq(from=1, to=num_m, by=2),]),1,function(x) x*(gsd_cmm2mc+gsd_cm2mmc)))%*%B_des_c_m[seq(from=1, to=num_m, by=2),] + t(apply(t(B_des_c_d[seq(from=1, to=num_d, by=2),]),1,function(x) x*(gsd_3dd3d3+gsd_3d3dd3)))%*%B_des_c_d[seq(from=1, to=num_d, by=2),]
gsd_e <- (-0.5)*gsd_e 


d_var_b_a <- ((-1)*(t(beta_a)%*%D_a%*%beta_a)-0.5*sum(diag(inv_gsd_max%*%gsd_b_a)))/(var_b_a^2) + (length(beta_a)-2)/var_b_a
d_var_b_c <- ((-1)*(t(beta_c)%*%D_c%*%beta_c)-0.5*sum(diag(inv_gsd_max%*%gsd_b_c)))/(var_b_c^2) + (length(beta_c)-2)/var_b_c
d_var <- var_D_m - var_YSY_m + var_D_d - var_YSY_d + sum(diag(inv_gsd_max%*%gsd_e))

return(c(d_var, d_var_b_a, d_var_b_c))
}
