## shtime_plant
# time_plots for 'plant'
library(deSolve)
source('flux_plant.R')
source('pars_plant.R')
##
#function shtime_plant(j)
  # created: 2000/09/26 by Bas Kooijman, modified 2009/01/05
  
  ## Syntax
  # <../shtime_plant.m *shtime_plant*> (j) 

  ## Description
  # time_plots for 'plant'
  #
  # Input:
  #
  # * j: optional scalar with plot number (default all)
  
  ## Remark
  # First run pars_plant See <plant.hmtl *plant*>
  
  # global X T_1 T_A T_L T_H T_AL T_AH n_N_ENR n_N_ER
  # global M_VSd M_VSm M_VRd M_VRm M_VSb M_VRb M_VSp M_ER0
  # global k_C k_O k_ECS k_ENS k_ES k_ECR k_ENR k_ER rho_NO
  # global J_L_K K_C K_O K_NH K_NO K_H
  # global j_L_Am j_C_Am j_O_Am j_NH_Am j_NO_Am
  # global y_ES_CH_NO y_CH_ES_NO y_ER_CH_NO y_CH_ER_NO y_ER_CH_NH
  # global y_VS_ES y_ES_VS y_VR_ER y_ER_VR y_ES_ER y_ER_ES
  # global y_ES_ENS y_ENS_ES y_ER_ENR y_ENR_ER y_ENS_ENR y_ECR_ECS
  # global kap_ECS kap_ECR kap_ENS kap_ENR kap_SS
  # global kap_SR kap_RS kap_TS kap_TR
  # global j_ES_MS j_ER_MR j_ES_JS j_PS_MS j_PR_MR y_PS_VS y_PR_VR
 
tmax = 500 # select time period

# State vector:
#   M = [M_PS, M_VS, M_ECS, M_ENS, M_ES, M_PR, M_VR, M_ECR, M_ENR, M_ER]
M0 = c(0, 1e-4, 1e-4, 1e-4, 1e-4, 0, 1e-4, 1e-4, 1e-4, 10) # initial value

t_Mt  = ode(M0, seq(0,tmax), flux_plant, p, method="ode45") # integrate
t <- t_Mt[,1]
colnames(t_Mt) <- NULL
Mt <- t_Mt[,2:11]
Mt2 <- t_Mt[,12:91]
head(t_Mt)
(t(t_Mt)[,1])

par(mfrow = c(2,2))
plot(t, Mt[,2], col = 'green', xlab='time, d', ylab='structures, mol', type = 'l', ylim = c(-150, 100)) 
points(t, -Mt[,7], col='red', type = 'l')
abline(1,0)

plot(t, Mt[,1], col = 'green', xlab='time, d', ylab='products, mol', type = 'l', ylim = c(-400, 400)) 
points(t, -Mt[,6], col='red', type = 'l')      
abline(1,0)

plot(t, Mt[,3], col = 'black', xlab='time, d', ylab='reserves, mol', type = 'l', ylim = c(-1500, 500)) 
points(t, Mt[,4], col='blue', type = 'l')  
points(t, Mt[,5], col='red', type = 'l')  
points(t, -Mt[,8], col='black', type = 'l')  
points(t, -Mt[,9], col='blue', type = 'l')  
points(t, -Mt[,10], col='red', type = 'l')  

m_ECR = Mt[,8]/Mt[,7]; m_ENR = Mt[,9]/Mt[,7];
m_ER = Mt[,10]/Mt[,7]; m_ERm = max(t(m_ECR),t(m_ENR));
plot(t, Mt[,3]/Mt[,2], col = 'black', xlab='time, d', ylab='reserve densities, mol/mol', type = 'l', ylim = c(-10, 15)) 
points(t, Mt[,4]/Mt[,2], col='blue', type = 'l')  
points(t, Mt[,5]/Mt[,2], col='red', type = 'l')  
points(t, -m_ECR, col='black', type = 'l')  
points(t, -m_ENR, col='blue', type = 'l')  
points(t, -pmin(m_ER,m_ERm), col='red', type = 'l')  


par(mfrow = c(2,2))
plot(t, Mt[,2] * w_VS, col = 'green', xlab='time, d', ylab='structures, g', type = 'l', ylim = c(-150 * 25, 100 * 25)) 
points(t, -Mt[,7] * w_VR, col='red', type = 'l')
abline(1,0)

plot(t, Mt[,1] * w_PS, col = 'green', xlab='time, d', ylab='products, g', type = 'l', ylim = c(-400 * 25, 400 * 25)) 
points(t, -Mt[,6] * w_PR, col='red', type = 'l')      
abline(1,0)

plot(t, Mt[,3] * w_ECS, col = 'black', xlab='time, d', ylab='reserves, g', type = 'l', ylim = c(-1500 * 25, 500 * 25)) 
points(t, Mt[,4] * w_ENS, col='blue', type = 'l')  
points(t, Mt[,5] * w_ES, col='red', type = 'l')  
points(t, -Mt[,8] * w_ECR, col='black', type = 'l')  
points(t, -Mt[,9] * w_ENR, col='blue', type = 'l')  
points(t, -Mt[,10] * w_ER, col='red', type = 'l')  
abline(1,0)

m_ECR = Mt[,8]/Mt[,7]; m_ENR = Mt[,9]/Mt[,7];
m_ER = Mt[,10]/Mt[,7]; m_ERm = max(t(m_ECR),t(m_ENR));
plot(t, Mt[,3]/Mt[,2], col = 'black', xlab='time, d', ylab='reserve densities, mol/mol', type = 'l', ylim = c(-1, 1.5)) 
points(t, Mt[,4]/Mt[,2], col='blue', type = 'l')  
points(t, Mt[,5]/Mt[,2], col='red', type = 'l')  
points(t, -m_ECR, col='black', type = 'l')  
points(t, -m_ENR, col='blue', type = 'l')  
points(t, -pmin(m_ER,m_ERm), col='red', type = 'l')  
abline(1,0)

