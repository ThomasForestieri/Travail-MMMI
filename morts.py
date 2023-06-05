import json
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Charger les donnees a partir du fichier JSON
with open('dataMorts.json') as f:
    data = json.load(f)

# Dictionnaire pour stocker les valeurs pour chaque tranche d'age et chaque jour
valeurs_par_tranche_age = {
    "0-19": {},
    "20-64": {},
    "65+": {}
}

# Parcourir les donnees et extraire les valeurs pour chaque tranche d'age et chaque jour
for entry in data["data"]:
    date = entry["date"]
    if "2020-03-19" <= date <= "2020-12-31":
        cum_admissions_by_age = entry["newDeaths28DaysByDeathDateAgeDemographics"]
        for age_group in cum_admissions_by_age:
            age = age_group["age"]
            value = age_group["deaths"]
            if age == "00_04" or age == "05_09" or age == "10_14" or age == "15_19":
                if date not in valeurs_par_tranche_age["0-19"]:
                    valeurs_par_tranche_age["0-19"][date] = []
                valeurs_par_tranche_age["0-19"][date].append(value)
            elif age == "20_24" or age == "25_29" or age == "30_34" or age =="35_39" or age == "40_44" or age == "45_49" or age =="50_54" or age =="55_59" or age =="60_64":
                if date not in valeurs_par_tranche_age["20-64"]:
                    valeurs_par_tranche_age["20-64"][date] = []
                valeurs_par_tranche_age["20-64"][date].append(value)
            elif age == "65_69" or age == "70_74" or age == "75_79" or age =="80_84" or age == "85_89" or age == "90+":
                if date not in valeurs_par_tranche_age["65+"]:
                    valeurs_par_tranche_age["65+"][date] = []
                valeurs_par_tranche_age["65+"][date].append(value)

# Fonction pour calculer la quantite totale par jour de morts
def calculate_daily_total(data):
    dates = sorted(data.keys())
    total = 0
    daily_total = []

    for i in range(1, len(dates)):
        current_date = dates[i]
        total +=  sum(data[current_date])
        daily_total.append((current_date, total))
    return daily_total

# Calcul de la quantite de mort par jour pour chaque tranche d'age
daily_total_0_19 = calculate_daily_total(valeurs_par_tranche_age["0-19"])
daily_total_20_64 = calculate_daily_total(valeurs_par_tranche_age["20-64"])
daily_total_65_plus = calculate_daily_total(valeurs_par_tranche_age["65+"])

# Dictionnaire pour stocker les valeurs combinees de chaque tranche d'age
valeurs_combinees = {
    "0-19": {},
    "20-64": {},
    "65+": {}
}


for date, total in daily_total_65_plus:    
    # Ajout de la quantite ajoutee par jour pour la tranche d'age 65+
    if date in valeurs_par_tranche_age["65+"]:
        valeurs_combinees["65+"][date] = (sum(valeurs_par_tranche_age["65+"][date]), total)


for date, total in daily_total_0_19:
    # Ajout de la quantite ajoutee par jour pour la tranche d'age 0-17
    if date in valeurs_par_tranche_age["0-19"]:
        valeurs_combinees["0-19"][date] = (sum(valeurs_par_tranche_age["0-19"][date]), total)
    

for date, total in daily_total_20_64:
    # Ajout de la quantite ajoutee par jour pour la tranche d'age 18-64
    if date in valeurs_par_tranche_age["20-64"]:
        valeurs_combinees["20-64"][date] = (sum(valeurs_par_tranche_age["20-64"][date]), total)

""""
# Affichage des valeurs combinees pour chaque tranche d'age
print("Tranche d'age de 0 a 19 ans :")
for date, (value, total) in valeurs_combinees["0-19"].items():
    print("     Date:", date)
    print("     Valeur:", value)
    print("     Quantite totale:", total)
    print()


print("Tranche d'age de 20 a 64 ans :")
for date, (value, total) in valeurs_combinees["20-64"].items():
    print("     Date:", date)
    print("     Valeur:", value)
    print("     Quantite totale:", total)
    print()


print("Tranche d'age de 65 ans et plus :")
for date, (value, total) in valeurs_combinees["65+"].items():
    print("     Date:", date)
    print("     Valeur:", value)
    print("     Quantite totale:", total)
    print()
"""
""""
# Dates pour l'axe des x
dates = list(valeurs_combinees["0-19"].keys())
dates = [mdates.datestr2num(date) for date in dates]

# Tracer les graphiques
plt.plot(dates, [valeurs_combinees["0-19"][date][1] for date in valeurs_combinees["0-19"].keys()], label="0-19 ans")
plt.plot(dates, [valeurs_combinees["20-64"][date][1] for date in valeurs_combinees["20-64"].keys()], label="20-64 ans")
plt.plot(dates, [valeurs_combinees["65+"][date][1] for date in valeurs_combinees["65+"].keys()], label="65+ ans")

# Configurer le graphique
plt.title("Deces entre le 19 mars 2020 et le 31 decembre 2020")
plt.xlabel("Date")
plt.xticks(rotation=45)
plt.ylabel("Nombre de morts")
plt.legend()

plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))

# Afficher le graphique
plt.show()"""