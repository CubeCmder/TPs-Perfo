import numpy as np
try:
    from Decollage import *
    from Longueur_Piste import *
except:
    pass

Poids = input("Entrer la valeur du poids de l'avion (lbs) : ")
Hp = input("Entrer la valeur de l'altitude pression (ft) : ")
ISA = input("Entrer la valeur de la déviation par rapport à ISA (C) : ")

valeurs_decollage = Decollage(float(Poids), float(Hp), float(ISA))

print("\nRésultats pour le décollage :")
print("-" * 60)
print("V1 minimum (ft/s) : " + str("{:.5g}".format(valeurs_decollage[0])))
print("V1 maximum (ft/s) : " + str("{:.5g}".format(valeurs_decollage[1])))
print("VR (ft/s) : " + str("{:.5g}".format(valeurs_decollage[2])))
print("V2 (ft/s) : " + str("{:.5g}".format(valeurs_decollage[3])))
print("Vloei (ft/s) : " + str("{:.5g}".format(valeurs_decollage[4])))
print("Vloaeo (ft/s) : " + str("{:.5g}".format(valeurs_decollage[5])))
print("V35aeo (ft/s) : " + str("{:.5g}".format(valeurs_decollage[6])))
print("delta_t (vlo-vr)oei (s) : " + str("{:.5g}".format(valeurs_decollage[7])))
print("delta_t (v35-vlo)oei (s) : " + str("{:.5g}".format(valeurs_decollage[8])))
print("delta_t (vlo-vr)aeo (s) : " + str("{:.5g}".format(valeurs_decollage[9])))
print("delta_t (v35-vlo)aeo (s) : " + str("{:.5g}".format(valeurs_decollage[10])))
print("Dist (vlo-vr)oei (ft) : " + str("{:.5g}".format(valeurs_decollage[11])))
print("Dist (v35-vlo)oei (ft) : " + str("{:.5g}".format(valeurs_decollage[12])))
print("Dist (vlo-vr)aeo (ft) : " + str("{:.5g}".format(valeurs_decollage[13])))
print("Dist (v35-vlo)aeo (ft) : " + str("{:.5g}".format(valeurs_decollage[14])))
print("-" * 60)
print("\n")

V1VR = input("Entrer la valeur du rapport V1/VR : ")

valeurs_longueur_piste = Longueur_Piste(float(V1VR),float(Poids), float(Hp), float(ISA))

print("\nRésultats pour la longueur de la piste :")
print("-" * 60)
print("FTOD AEO (ft) : " + str("{:.5g}".format(valeurs_longueur_piste[0])))
print("TOD OEI (ft) : " + str("{:.5g}".format(valeurs_longueur_piste[1])))
print("ASD AEO (ft) : " + str("{:.5g}".format(valeurs_longueur_piste[2])))
print("Longueur minimum requise (ft) : " + str("{:.5g}".format(valeurs_longueur_piste[3])))
print("Énergie totale absorbée par les freins pendant l’arrêt (millions ft-lb) : " + str("{:.5g}".format(valeurs_longueur_piste[4])))
print("-" * 60)
print("\n")