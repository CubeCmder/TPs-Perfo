AER8375 - Analyse et performance des avions
Laboratoire 4
Auteurs : Robert Fercal et Titouan Dallemagne
Groupe 04
==================================================================
==================================================================
Le dossier est composé de 9 fichiers principaux :

1. Constantes.py
	- Ce fichier contient toutes les constantes qui sont utilisées pour les calculs.

2. Transformations.py
	- Ce fichier contient toutes les formules pour les transformations nécessaires (ex. Celsius --> Kelvin).

3. Atmosphere.py
	- Ce fichier permet de calculer les paramètres atmosphériques (pression, theta, sigma, etc.) en fournissant
	soit la valeur de l'hauteur pression (pi) ainsi que la température (Celsius) ou la déviation ISA (Celsius).

4. Vitesses.py
	- Ce fichier permet de calculer les paramètres reliés à la vitesse (Vc, V, Ve, T totale, etc) en fournissant
	la vitesse désirée (Vc, V, Ve ou le nombre de Mach) ainsi que le poids de l'avion (lbs), la pression, la température, 
	delta, theta, sigma, rho, S, MAC et K. Beaucoup de ces valeurs proviennent du fichier Atmosphere.py.

5. Forces.py
	- Ces fichiers permet de calculer les paramètres reliés à aux forces (Cd, Cl, L/D, Cdp, Cdi, etc.) en fournissant
	les paramètres tels que l'hauteur pression, la température (T ou ISA), la vitesse, le poids de l'avion, le CG, etc.
	Beaucoup de ces valeurs utilisées par la fonction viennent du fichier Atmosphere.py.

6. Decollage.py
	- Ce fichier permet de calculer les valeurs en lien avec le décollage, tel que V1 minimum et maximum, 
	la vitesse de rotation VR, etc.

7. Longueur_Piste.py
	- Ce fichier permet de calculer les valeurs en lien avec la longueur de piste requise pour le décollage.

8. Resultats.py
	- Ce fichier permet à l'utilisateur de trouver les paramètres recherchés en fournissant les valeurs demandés.
	Vous n'avez que à lancer le programme et suivre les instructions dans la console. Aucun changement dans le code
	n'est nécessaire.

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Il est primordial de mettre les paramètres tels que spécifiés entre parenthèses 
	et des valeurs numériques lorsque des valeurs numériques sont demandés. Certains
	handling exceptions ont été implémentés, comme pour le choix de T ou ISA. Le
	programme ne vous laisse pas mettre autre chose que ces valeurs là. Il est donc 
	toujours important de suivre les parenthèses qui se trouvent après les input
	question.
 
	Entrer la valeur du poids de l'avion (lbs) : 45000
	Entrer la valeur de l'altitude pression (ft) : 0
	Entrer la valeur de la déviation par rapport à ISA (C) : 20

	Résultats pour le décollage :
	------------------------------------------------------------
	V1 minimum (ft/s) : 165.81
	V1 maximum (ft/s) : 221.92
	VR (ft/s) : 221.92
	V2 (ft/s) : 231.83
	Vloei (ft/s) : 226.63
	Vloaeo (ft/s) : 233.3
	V35aeo (ft/s) : 250.93
	delta_t (vlo-vr)oei (s) : 2.026
	delta_t (v35-vlo)oei (s) : 5.1802
	delta_t (vlo-vr)aeo (s) : 1.3
	delta_t (v35-vlo)aeo (s) : 3.8
	Dist (vlo-vr)oei (ft) : 454.39
	Dist (v35-vlo)oei (ft) : 1187.5
	Dist (vlo-vr)aeo (ft) : 295.89
	Dist (v35-vlo)aeo (ft) : 920.05
	------------------------------------------------------------


	Entrer la valeur du rapport V1/VR : 0.9

	Résultats pour la longueur de la piste :
	------------------------------------------------------------
	FTOD AEO (ft) : 4508.1
	TOD OEI (ft) : 5236.6
	ASD AEO (ft) : 4125.6
	Longueur minimum requise (ft) : 4508.1
	Énergie totale absorbée par les freins pendant l’arrêt (millions ft-lb) : 27.265
	------------------------------------------------------------

	L'utilisateur n'a pas à modifier lui-même aucun des codes. Il suffit seulement
	de lancer le fichier Resultats.py et suivre les instructions.

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

9. TP4_Q4.py
	- 

	Lancer simplement le code sans modifier quoi que ce soit. Les valeurs demandées apparaîtront
	dans la console.