## Importation des packages

library(deSolve)
library(mcmc)

## importation des donnees

death <- read.csv("deathsCum.csv", header=TRUE, sep=",")

dates_brutes <- unlist(death[1], use.names = FALSE)
dates = strptime(dates_brutes, "%Y-%m-%d")

d_0_19 <- unlist(death[2], use.names = FALSE)
d_20_64 <- unlist(death[3], use.names = FALSE)
d_65plus <- unlist(death[4], use.names = FALSE)


##############################
# Plot des donnees de deces ##
##############################

# du 20 mars 2020 au 31 decembre 2020
nb_days = length(d_0_19)

axis.POSIXct(1, at=seq(from=strptime(dates[1],"%Y-%m-%d" ), to=strptime(dates[nb_days],"%Y-%m-%d" ), by="month"), format="%b", las=1)

par(mfrow = c(1,3))
plot(dates, d_0_19,xlab="Temps (mois)",ylab="décès (0-19ans)",col="green")
plot(dates, d_20_64,xlab="Temps (mois)", ylab="décès (20-64ans)", type = "o",col="blue")
plot(dates, d_65plus,xlab="Temps (mois)", ylab="décès (65+ ans)", type = "o",col="red")

##############################
#### Parametres du modele ####
##############################

# Population totale pour chaque classe d'age (0-19, 20-64, 65+)
# donnees issues de https://perspective.usherbrooke.ca/bilan/servlet/BMPagePyramide/GBR/2018/? (pour le nombre) et de
# https://www.larousse.fr/encyclopedie/divers/Grande-Bretagne_population/185857 pour la proportion de personnes en Angleterre
N = c(12998907, 30215341, 10147838)

# Nombre de reproduction de base
# confinement le 24 mars pour 3 semaines
# nouvement confinement du 5 novembre au 2 decembre
# reconfinement le 20 decembre
R_0 = c(2.5,1,2,1)

# duree moyenne de la periode infectieuse
gamma = c(1/8, 1/12, 1/15) 

contact_matrix = matrix(c(7, 5, 2, 5, 6, 4, 2, 4, 5), nrow = 3, ncol = 3)

vp = eigen(contact_matrix/gamma)
mean_vp = max(vp$values)

# Taux de transmission (passage au compartiment E)
beta1 = R_0[1]/mean_vp*c(1,1)      
beta2 = R_0[2]/mean_vp*c(1,1)
beta3 = R_0[3]/mean_vp*c(1,1)
beta4 = R_0[4]/mean_vp*c(1,1)

# Taux d'infection (passage au compartiment I)
i = c(0.01, 0.1, 0.5)

# taux de guerison (passage de I a R)
r = c(0.97, 0.70, 0.50)

param = list(beta1, beta2, beta3, beta4,i, r, gamma)

# Compartiments de base

E0 = c(1000,2000,1500)
I0 = c(1/((1-r[1])*gamma[1]),d_20_64[1]/((1-r[2])*gamma[2]),d_65plus[1]/((1-r[3])*gamma[3]))
D0 = c(1,d_20_64[1],d_65plus[1]) 
R0 = c(r[1]*gamma[1]*I0[1],r[2]*gamma[2]*I0[2],r[3]*gamma[3]*I0[3])
S0 = c(N[1] - E0[1] - I0[1] - D0[1] - R0[1], N[2]- E0[2] - I0[2] - D0[2] - R0[2], N[3]- E0[3] - I0[3] - D0[3] - R0[3])

Comparts = list(S0, E0, I0, D0, R0)


##############################
######## modele SEIRD ########
##############################

t = 0:(nb_days-1)

SEIRD = function(t, compart, param ){
  
  S1 = compart[1]
  S2 = compart[2]
  S3 = compart[3]
  E1 = compart[4]
  E2 = compart[5]
  E3 = compart[6]
  I1 = compart[7]
  I2 = compart[8]
  I3 = compart[9]
  R1 = compart[10]
  R2 = compart[11]
  R3 = compart[12]
  D1 = compart[13]
  D2 = compart[14]
  D3 = compart[15]
  
  beta <- if (t <= 4) {param[[1]]}
  else if (t <= 25) {param[[2]]}
  else if (t <= 221) {param[[3]]}
  else {param[[4]]}
  
  i = param[[5]] 
  r = param[[6]]
  gamma = param[[7]]
  
  dS1 = -(beta[1])*S1*(contact_matrix[1]*I1/(N[1])+contact_matrix[2]*I2/(N[2])+contact_matrix[3]*I3/(N[3]))
  dS2 = -(beta[1])*S2*(contact_matrix[4]*I1/(N[1])+contact_matrix[5]*I2/(N[2])+contact_matrix[6]*I3/(N[3]))
  dS3 = -(beta[1])*S3*(contact_matrix[7]*I1/(N[1])+contact_matrix[8]*I2/(N[2])+contact_matrix[9]*I3/(N[3]))
  
  dE1 = ((beta[1])*S1*(contact_matrix[1]*I1/(N[1])+contact_matrix[2]*I2/(N[2])+contact_matrix[3]*I3/(N[3])) - i[1] * E1)
  dE2 = ((beta[1])*S2*(contact_matrix[4]*I1/(N[1])+contact_matrix[5]*I2/(N[2])+contact_matrix[6]*I3/(N[3])) - i[2] * E2)
  dE3 = ((beta[1])*S3*(contact_matrix[7]*I1/(N[1])+contact_matrix[8]*I2/(N[2])+contact_matrix[9]*I3/(N[3])) - i[3] * E3)
  
  dI1 = i[1] * E1 - gamma[1]*I1
  dI2 = i[2] * E2 - gamma[2]*I2
  dI3 = i[3] * E3 - gamma[3]*I3
  
  dR1 = r[1]*gamma[1]*I1
  dR2 = r[2]*gamma[2]*I2
  dR3 = r[3]*gamma[3]*I3
  
  dD1 = (1-r[1])*gamma[1]*I1
  dD2 = (1-r[2])*gamma[2]*I2
  dD3 = (1-r[3])*gamma[3]*I3
  
  dX = c(dS1,dS2,dS3, dE1, dE2, dE3, dI1,dI2, dI3, dR1,dR2, dR3, dD1, dD2, dD3)
  return(list(dX))
}

theta0 = c(Comparts[[2]][1],Comparts[[2]][2], Comparts[[2]][3], #E
           Comparts[[3]][1],Comparts[[3]][2], Comparts[[3]][3], #I
           Comparts[[4]][1],Comparts[[5]][2], Comparts[[5]][3], #R
           Comparts[[5]][1],Comparts[[5]][2], Comparts[[5]][3], #D
           param[[1]][1], param[[2]][1], param[[3]][1],param[[4]][1], #beta
           param[[5]][1], param[[5]][2], param[[5]][3], #i
           param[[6]][1], param[[6]][2], param[[6]][3], #r
           param[[7]][1], param[[7]][2], param[[7]][3]) #gamma

CI_X = c(N[1]-theta0[1]-theta0[4]-theta0[7]-theta0[10],
         N[2]-theta0[2]-theta0[5]-theta0[8]-theta0[11],
         N[3]-theta0[3]-theta0[6]-theta0[9]-theta0[12],
         theta0[1],theta0[2],theta0[3],theta0[4],theta0[5], theta0[6],
         theta0[7],theta0[8],theta0[9],theta0[10],theta0[11], theta0[12])


param_X = list(theta0[13],theta0[14],theta0[15],theta0[16],c(theta0[17],theta0[18], theta0[19]),
               c(theta0[20],theta0[21], theta0[22]),c(theta0[23],theta0[24], theta0[25]))

X = ode(CI_X,t,SEIRD,param_X, method="bdf")

plot(X)


################################################
## fonction de calcul de la log-vraisemblance ##
################################################

log_l = function(log_t){
  
  theta = exp(log_t)
  
  if (theta[17] > 1 || theta[18] > 1 || theta[19] > 1 || theta[20] > 1 || theta[21] > 1 || theta[22] > 1 || theta[23] > 1 || theta[24] > 1 || theta[25] > 1) {return(-Inf)}
  
  CI = c(N[1]-theta[1]-theta[4]-theta[7]-theta[10],
           N[2]-theta[2]-theta[5]-theta[8]-theta[11],
           N[3]-theta[3]-theta[6]-theta[9]-theta[12],
           theta[1],theta[2],theta[3],theta[4],theta[5], theta[6],
           theta[7],theta[8],theta[9],theta[10],theta[11], theta[12])
  
  param = list(theta[13],theta[14],theta[15],theta[16],c(theta[17],theta[18], theta[19]),
                 c(theta[20],theta[21], theta[22]),c(theta[23],theta[24], theta[25]))
  
  X = ode(CI,t,SEIRD,param, method="bdf")
  
  d1 = X[,14]
  d2 = X[,15]
  d3 = X[,16]
  
  logL_d1 = dpois(d_0_19,d1,log=T)
  logL_d2 = dpois(d_20_64,d2,log=T)
  logL_d3 = dpois(d_65plus,d3,log=T)

  logL = sum(c(logL_d1,logL_d2,logL_d3))
  
  return(logL)
}

theta0 = log(c(Comparts[[2]][1],Comparts[[2]][2], Comparts[[2]][3], #E
               Comparts[[3]][1],Comparts[[3]][2], Comparts[[3]][3], #I
               Comparts[[4]][1],Comparts[[5]][2], Comparts[[5]][3], #R
               Comparts[[5]][1],Comparts[[5]][2], Comparts[[5]][3], #D
               param[[1]][1], param[[2]][1], param[[3]][1],param[[4]][1], #beta
               param[[5]][1], param[[5]][2], param[[5]][3], #i
               param[[6]][1], param[[6]][2], param[[6]][3], #r
               param[[7]][1], param[[7]][2], param[[7]][3])) #gamma)

## Optimisation de la log-vraisemblance

opt <- list(optim(theta0,log_l,control=list(fnscale=-1)))
opt[[1]]$par = exp(opt[[1]]$par)

## Resultats de l'optimisation

S01 = N[1]-opt[[1]]$par[1]-opt[[1]]$par[4] -opt[[1]]$par[7] -opt[[1]]$par[10];
S02 = N[2]-opt[[1]]$par[2]-opt[[1]]$par[5] -opt[[1]]$par[8] -opt[[1]]$par[11];
S03 = N[3]-opt[[1]]$par[3]-opt[[1]]$par[6] -opt[[1]]$par[9] -opt[[1]]$par[12];
E01 = opt[[1]]$par[1];
E02 = opt[[1]]$par[2];
E03 = opt[[1]]$par[3];
I01 = opt[[1]]$par[4];
I02 = opt[[1]]$par[5];
I03 = opt[[1]]$par[6];
R01 = opt[[1]]$par[7];
R02 = opt[[1]]$par[8];
R03 = opt[[1]]$par[9];
D01 = opt[[1]]$par[10];
D02 = opt[[1]]$par[11];
D03 = opt[[1]]$par[12];

CI = c(S01,S02, S03, E01, E02, E03, I01, I02, I03, R01, R02, R03, D01, D02, D03);

beta1 = opt[[1]]$par[13];
beta2 = opt[[1]]$par[14];
beta3 = opt[[1]]$par[15];
beta4 = opt[[1]]$par[16];

i = c(opt[[1]]$par[17], opt[[1]]$par[18], opt[[1]]$par[19]);
r = c(opt[[1]]$par[20], opt[[1]]$par[21], opt[[1]]$par[22]);
gamma = c(opt[[1]]$par[23], opt[[1]]$par[24], opt[[1]]$par[25]);

param = list(beta1,beta2,beta3,beta4,i,r,gamma) 

## Simulation du nouveau modele optimal

X_opt = ode(CI,t,SEIRD,param, method="bdf")

# Graphique

par(mfrow = c(1,3))

plot(dates, d_0_19,xlab="Temps (mois)",ylab="décès totaux (0-19ans)",col="blue")
lines(dates,X_opt[,14], col = "red", lwd = 2)
legend("topleft", legend = c("Donnees (0-19ans)", "Modele (0-19ans)"), col = c("blue", "red"), lty = c(0,1), lwd = c(1,2), pch = c(1,NA))

plot(dates, d_20_64,xlab="Temps (mois)", ylab="décès totaux (20-64ans)", type = "o",col="blue")
lines(dates,X_opt[,15], col = "red", lwd = 2)
legend("topleft", legend = c("Donnees (20-64ans)", "Modele (20-64ans)"), col = c("blue", "red"), lty = c(0,1), lwd = c(1,2), pch = c(1,NA))

plot(dates, d_65plus,xlab="Temps (mois)", ylab="décès totaux (65+ ans)", type = "o",col="blue")
lines(dates,X_opt[,16], col = "red", lwd = 2)
legend("topleft", legend = c("Donnees (65+ ans)", "Modele (65+ ans)"), col = c("blue", "red"), lty = c(0,1), lwd = c(1,2), pch = c(1,NA))


##########
## MCMC ##
##########
