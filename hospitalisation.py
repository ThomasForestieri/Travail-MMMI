import json
import matplotlib.dates as mdates

# Charger les donnees a partir du fichier JSON
with open('dataHospitalisation.json') as f:
    data = json.load(f)

# Dictionnaire pour stocker les valeurs pour chaque tranche d'age et chaque jour
valeurs_par_tranche_age = {
    "0-17": {},
    "18-64": {},
    "65+": {}
}

# Parcourir les donnees et extraire les valeurs pour chaque tranche d'age et chaque jour
for entry in data["data"]:
    date = entry["date"]
    if "2020-03-19" <= date <= "2020-12-31":
        cum_admissions_by_age = entry["cumAdmissionsByAge"]
        for age_group in cum_admissions_by_age:
            age = age_group["age"]
            value = age_group["value"]
            if age == "0_to_5" or age == "6_to_17":
                if date not in valeurs_par_tranche_age["0-17"]:
                    valeurs_par_tranche_age["0-17"][date] = []
                valeurs_par_tranche_age["0-17"][date].append(value)
            elif age == "18_to_64":
                if date not in valeurs_par_tranche_age["18-64"]:
                    valeurs_par_tranche_age["18-64"][date] = []
                valeurs_par_tranche_age["18-64"][date].append(value)
            elif age == "65_to_84" or age == "85+":
                if date not in valeurs_par_tranche_age["65+"]:
                    valeurs_par_tranche_age["65+"][date] = []
                valeurs_par_tranche_age["65+"][date].append(value)

# Fonction pour calculer la quantite ajoutee par jour
def calculate_daily_increase(data):
    dates = sorted(data.keys())
    daily_increase = []
    for i in range(1, len(dates)):
        current_date = dates[i]
        previous_date = dates[i - 1]
        current_value = sum(data[current_date])
        previous_value = sum(data[previous_date])
        increase = current_value - previous_value
        daily_increase.append((current_date, increase))
    return daily_increase

# Calcul de la quantite ajoutee par jour pour chaque tranche d'age
daily_increase_0_17 = calculate_daily_increase(valeurs_par_tranche_age["0-17"])
daily_increase_18_64 = calculate_daily_increase(valeurs_par_tranche_age["18-64"])
daily_increase_65_plus = calculate_daily_increase(valeurs_par_tranche_age["65+"])

# Dictionnaire pour stocker les valeurs combinees de chaque tranche d'age
valeurs_combinees = {
    "0-17": {},
    "18-64": {},
    "65+": {}
}


for date, increase in daily_increase_0_17:
    # Ajout de la quantite ajoutee par jour pour la tranche d'age 0-17
    if date in valeurs_par_tranche_age["0-17"]:
        valeurs_combinees["0-17"][date] = (sum(valeurs_par_tranche_age["0-17"][date]), increase)


for date, increase in daily_increase_18_64:
    # Ajout de la quantite ajoutee par jour pour la tranche d'age 18-64
    if date in valeurs_par_tranche_age["18-64"]:
        valeurs_combinees["18-64"][date] = (sum(valeurs_par_tranche_age["18-64"][date]), increase)

for date, increase in daily_increase_65_plus:
    # Ajout de la quantite ajoutee par jour pour la tranche d'age 65+
    if date in valeurs_par_tranche_age["65+"]:
        valeurs_combinees["65+"][date] = (sum(valeurs_par_tranche_age["65+"][date]), increase)
