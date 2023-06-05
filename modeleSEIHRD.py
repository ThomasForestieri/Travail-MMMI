import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import scipy.stats as ss
from scipy.optimize import minimize
import hospitalisation
import morts


########################
# Parametres du modele #
########################

N = np.array([12798907, 32215341, 10147838])  # Population totale pour chaque classe d'age (0-17, 18-64, 65+)
# donnees issues de https://perspective.usherbrooke.ca/bilan/servlet/BMPagePyramide/GBR/2018/? (pour le nombre) et de
# https://www.larousse.fr/encyclopedie/divers/Grande-Bretagne_population/185857 pour la proportion de personnes en Angleterre

beta = np.array([0.02, 0.04, 0.06])  # Taux de transmission pour chaque classe d'age
sigma = np.array([1/10, 1/10, 1/10])  # Taux d'incubation pour chaque classe d'age
gamma = np.array([1/100, 1/100, 1/100])  # Taux de guerison pour chaque classe d'age
alpha = np.array([1/100, 1/80, 1/60])  # Taux de progression des cas exposes vers les cas infectieux pour chaque classe d'age
delta = np.array([1/300, 1/80, 1/50])  # Taux d'hospitalisation pour chaque classe d'age
rho = np.array([1/1000, 1/400, 1/100])  # Taux de sortie de l'hopital pour chaque classe d'age

# Conditions initiales
E0 = np.array([1000, 500, 200])  # Cas exposes initiaux pour chaque classe d'age
I0 = np.array([500, 2000, 1000])  # Cas infectieux initiaux pour chaque classe d'age
R0 = np.array([0, 0, 0])  # Retablis initiaux pour chaque classe d'age
H0 = np.array([0, 0, 0])  # Hospitalisations initiales pour chaque classe d'age
D0 = np.array([0, 0, 0])  # Deces initiaux pour chaque classe d'age
S0 = np.array([N[i] - E0[i] - I0[i] - R0[i] - H0[i] for i in range(3)])  # Susceptibles initiaux pour chaque classe d'age

# Temps
t = np.linspace(0, 287, 287)  # du 19 mars 2020 au 31 decembre 2020

# Conditions initiales regroupees
y0 = np.concatenate((S0, E0, I0, R0, H0, D0))

# Modele SEIR avec hospitalisations
def seir_model(y, t, N, beta, sigma, gamma, alpha, delta, rho):
    S = y[:3]
    E = y[3:6]
    I = y[6:9]
    R = y[9:12]
    H = y[12:15]
    D = y[15:]
    dSdt = [-beta[i] * S[i] * (I[i] + alpha[i] * E[i]) / N[i] for i in range(3)]
    dEdt = [beta[i] * S[i] * (I[i] + alpha[i] * E[i]) / N[i] - sigma[i] * E[i] for i in range(3)]
    dIdt = [sigma[i] * E[i] - gamma[i] * I[i] - delta[i] * I[i] for i in range(3)]
    dHdt = [delta[i] * I[i] - rho[i] * H[i] for i in range(3)]
    dRdt = [gamma[i] * I[i] + rho[i] * H[i] for i in range(3)]
    dDdt = [rho[i] * H[i] for i in range(3)]
    return np.concatenate((dSdt, dEdt, dIdt, dRdt, dHdt, dDdt))

# Fonction de log-vraisemblance
def log_likelihood(theta):
    beta = theta[0:3]
    sigma = theta[3:6]
    gamma = theta[6:9]
    alpha = theta[9:12]
    delta = theta[12:15]
    rho = theta[15:]
    result = odeint(seir_model, y0, t, args=(N, beta, sigma, gamma, alpha, delta, rho))
    H = result[:, 12:15]
    D = result[:, 15:]
    logL = np.sum(ss.poisson.logpmf(H, H_data))
    logL2 = np.sum(ss.poisson.logpmf(D, M_data))
    return logL+logL2

# Donnees d'hospitalisation
H_data = np.zeros((287, 3))
i = 0
for date, (value, increase) in hospitalisation.valeurs_combinees["0-17"].items():
    H_data[i][0] = value
    i += 1

i = 0
for date, (value, increase) in hospitalisation.valeurs_combinees["18-64"].items():
    H_data[i][1] = value
    i += 1

i = 0
for date, (value, increase) in hospitalisation.valeurs_combinees["65+"].items():
    H_data[i][2] = value
    i += 1

# Donnees de mortalite
M_data = np.zeros((287, 3))
i = 0
for date, (value, total) in morts.valeurs_combinees["0-19"].items():
    M_data[i][0] = total
    i += 1

i = 0
for date, (value, total) in morts.valeurs_combinees["20-64"].items():
    M_data[i][1] = total
    i += 1

i = 0
for date, (value, total) in morts.valeurs_combinees["65+"].items():
    M_data[i][2] = total
    i += 1

    
# Optimisation de la log-vraisemblance
theta0 = np.array([beta, sigma, gamma, alpha, delta, rho])  # Valeurs initiales des parametres
opt = minimize(log_likelihood, np.ravel(theta0))
theta_opt = opt.x
beta_opt = theta_opt[0:3]
sigma_opt = theta_opt[3:6]
gamma_opt = theta_opt[6:9]
alpha_opt = theta_opt[9:12]
delta_opt = theta_opt[12:15]
rho_opt = theta_opt[15:]

# Resolution des equations differentielles avec les parametres optimaux
result_opt = odeint(seir_model, y0, t, args = (N, beta_opt, sigma_opt, gamma_opt, alpha_opt, delta_opt, rho_opt))
S_opt = result_opt[:, :3]
E_opt = result_opt[:, 3:6]
I_opt = result_opt[:, 6:9]
R_opt = result_opt[:, 9:12]
H_opt = result_opt[:, 12:15]
D_opt = result_opt[:, 15:]


##################################
# Graphique pour des estimations #
##################################

# Creation des sous-graphiques
fig1, axes = plt.subplots(3, 1, figsize=(6, 10))

# Date de reference (19 mars 2020)
date_reference = datetime(2020, 3, 19)

# Nombre de jours entre la date de reference et chaque jour dans t
jours_depuis_reference = [date_reference + timedelta(days=int(t_i)) for t_i in t]

# Conversion des dates en nombres decimaux
jours_decimal = [mdates.datestr2num(date.strftime('%Y-%m-%d')) for date in jours_depuis_reference]


# Trace des courbes pour chaque classe d'age
age_labels = ['0-17 ans', '18-64 ans', '65+ ans']
for i, ax in enumerate(axes):
    #ax.plot(jours_decimal, result_opt[:, i] / 1000, label='Susceptibles')
    #ax.plot(jours_decimal, result_opt[:, i+3] / 1000, label='Exposes')
    #ax.plot(jours_decimal, result_opt[:, i+6] / 1000, label='Infectieux')
    #ax.plot(jours_decimal, result_opt[:, i+9] / 1000, label='Retablis')
    ax.plot(jours_decimal, result_opt[:, i+12] / 1000, label='Hospitalisations')
    ax.plot(jours_decimal, result_opt[:, i+15] / 1000, label='Deces')
    ax.set_xlabel("Date")

    ax.set_xticks(jours_decimal)
    ax.set_xticklabels([mdates.num2date(date).strftime('%b %Y') for date in jours_decimal], rotation=45)

    ax.set_ylabel('Nb personnes (en milliers)')
    ax.set_title('Classe d\'age {}'.format(age_labels[i]))
    ax.legend()

    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))

fig1.tight_layout()


# Creation des sous-graphiques
fig2, axes2 = plt.subplots(3, figsize=(6, 10))

#######################################
####### Graphique des donnees #########
#######################################

# Dates pour l'axe des x
dates = list(hospitalisation.valeurs_combinees["0-17"].keys())
dates = [mdates.datestr2num(date) for date in dates]

# Tracer les graphiques
axes2[0].plot(dates, [hospitalisation.valeurs_combinees["0-17"][date][0]/1000 for date in hospitalisation.valeurs_combinees["0-17"].keys()], label="Hospitalisation 0-17 ans")
axes2[0].plot(dates, [morts.valeurs_combinees["0-19"][date][1]/1000 for date in morts.valeurs_combinees["0-19"].keys()], label="Deces 0-19 ans")
axes2[0].set_title("Categorie 0-17/19 ans (donnees)")
axes2[0].set_xlabel("Date")
axes2[0].set_xticks(dates)
axes2[0].set_xticklabels([mdates.num2date(date).strftime('%b %Y') for date in dates], rotation=45)
axes2[0].set_ylabel("Nb personnes (en milliers)")
axes2[0].legend()

axes2[0].xaxis.set_major_locator(mdates.MonthLocator())
axes2[0].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))

axes2[1].plot(dates, [hospitalisation.valeurs_combinees["18-64"][date][0]/1000 for date in hospitalisation.valeurs_combinees["18-64"].keys()], label="Hospitalisation 18-64 ans")
axes2[1].plot(dates, [morts.valeurs_combinees["20-64"][date][1]/1000 for date in morts.valeurs_combinees["20-64"].keys()], label="Deces 20-64 ans")
axes2[1].set_title("Categorie 18/20-64 ans (donnees)")
axes2[1].set_xlabel("Date")
axes2[1].set_xticks(dates)
axes2[1].set_xticklabels([mdates.num2date(date).strftime('%b %Y') for date in dates], rotation=45)
axes2[1].set_ylabel("Nb personnes (en milliers)")
axes2[1].legend()

axes2[1].xaxis.set_major_locator(mdates.MonthLocator())
axes2[1].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))

axes2[2].plot(dates, [hospitalisation.valeurs_combinees["65+"][date][0]/1000 for date in hospitalisation.valeurs_combinees["65+"].keys()], label="Hospitalisation 65+ ans")
axes2[2].plot(dates, [morts.valeurs_combinees["65+"][date][1]/1000 for date in morts.valeurs_combinees["65+"].keys()], label="Deces 65+ ans")
axes2[2].set_title("Categorie 65+ ans (donnees)")
axes2[2].set_xlabel("Date")
axes2[2].set_xticks(dates)
axes2[2].set_xticklabels([mdates.num2date(date).strftime('%b %Y') for date in dates], rotation=45)
axes2[2].set_ylabel("Nb de personnes (en milliers)")
axes2[2].legend()

axes2[2].xaxis.set_major_locator(mdates.MonthLocator())
axes2[2].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))

fig2.tight_layout()

plt.tight_layout()
plt.show()
