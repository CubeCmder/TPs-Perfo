AER8375 - Analyse et performance des avions
Laboratoire 2
Auteurs : Robert Fercal et Titouan Dallemagne
Groupe 04
==================================================================
==================================================================
Le dossier est composé de 8 fichiers principaux :

1. Constantes.py
	- Ce fichier contient toutes les constantes qui sont utilisées pour les calculs.
2. Transformations.py
	- Ce fichier contient toutes les formules pour les transformations nécessaires (ex. Celsius --> Kelvin).
3. Atmosphere.py
	- Ce fichier permet de calculer les paramètres atmospheriques (pression, theta, sigma, etc.) en fournissant
	soit la valeur de l'hauteur pression (pi) ainsi que la température (Celsius) ou la déviation ISA (Celsius).
4. Forces.py !!!!!!!!!!!!!!!!! A MODIFIER !!!!!!!!!!!!!!!!
	- Ce fichier permet de calculer les paramètres reliés à aux forces (Cd, Cl, L/D, Cdp, Cdi, etc.) en fournissant
	les paramèetres tels que l'hauteur pression, la température (T ou ISA), la vitesse, le poids de l'avion, le CG, etc.
	Beacuop de ces valeurs utilisées par la fonction viennent du fichier Atmosphere.py. L'utilisateur n'a qu'à lancer
	le fichier Resultats.py et suivre les prompt. Aucune modification du code ne doit pas être faite par l'utilisateur.
5. Resultats.py
	- Ce fichier permet à l'utilisateur de trouver les paramètres recherchés en fournissant les valeurs demandés.
	Vous n'avez que à lancer le programme et suivre les instructions dans la console. Aucun changement dans le code
	n'est nécessaire.

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Il est primordial de mettre les paramètres tels que spécifiés entre paranthèses 
	et des valeurs numériques lorsque des valeurs numériques sont demandés. Certains
	handling exceptions ont été implémentés, comme pour le choix de T ou ISA. Le
	programme ne vous laisse pas mettre autre chose que ces valeurs là. Il est donc 
	toujours important de suivre les paranthèses qui se trouvent après les input
	question.
 
	Exemple : 

	Entrer la valeur de l'hauteur pression (ft) : 5000
	Choisir le paramètre à utiliser pour les calculs (T ou ISA): ISA
	Entrer la valeur de la déviation par rapport à delta ISA (C) : 0

	Résultats atmosphériques :
	------------------------------------------------------------
	Température (K) : 278.2
	Température (C) : 5.094
	Theta : 0.9656
	Delta : 0.8321
	Sigma : 0.8617
	Pression (lb/ft^2) : 1761
	Rho (slugs/ft^3) : 0.002048
	------------------------------------------------------------


	Entrer le choix de vitesse (Vc, Ve, V, Mach ou Vsr) : Vsr
	Entrer la vitesse basée sur votre choix précédent (kts ou sans unité) : 1.23
	Entrer le choix du régime (MCR, MTO, GA, Idle, MCT, MCL ou Poussee_totale) : MTO
	Entrer le nombre de moteurs (AEO, OEI) - Mettre n'importe lequel si vous avez choisi Poussee_totale : AEO
	Entrer le choix du régime pour la montée (Mach_constant ou CAS_constant) : CAS_constant
	Entrer le poids de l'avion (lbs) : 35000
	Entrer l'angle des flaps (0, 20, 45) : 20
	Entrer la position des flaps (UP, DOWN) : UP
	Entrer la position du centre de gravité (%MAC) : 0.15
	Choisir entre entrer la valeur de Nz ou f : Nz
	Entrer la valeur de Nz (g) : 1

	Résultats forces :
	------------------------------------------------------------
	CD : 0.096442
	CL : 1.2228
	Portance (lb) : 35000
	Finesse (L/D) : 12.679
	CDp : 0.0465
	Dp : 1330.9
	CDi : 0.049942
	Di : 1330.9
	CDcomp : 0
	Dcomp : 0
	CDcntl : 0
	Dcntl : 0
	CDwm : 0
	Dwm : 0
	Poussée totale (lb) : 12452
	Trainée (lb) : 2760.4
	Angle d'attaque (deg) : 9.5783
	Nz sw : 1.4158
	Phi sw : 45.063
	Nz buffet : 0
	Gamma montée : 27.017
	Taux de montée : 3758.2
	Taux de montée pression : 3758.2
	Facteur d'accéleration : 0.024969
	Accélération selon l’axe de la trajectoire de vol : 0.0067458
	------------------------------------------------------------

	L'utilisateur n'a pas à modifier lui-même aucun des codes. Il suffit seulement
	de lancer le fichier Resultats.py et suivre les instructions. Les résultats de 0, comme pour
	le Nz buffet, représentent en fait la valeur NA, soit aucun Nz buffet.

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

6. TP1_Q5.py
	- Permet de générer le graphique et la valeur de l'intersection qui se trouve dans la réponse de 
	la question 5 du mini-rapport.
	Lancer simplement le code sans modifer quoi que ce soit.
7. TP1_Q6.py
	- Permet de générer le graphique et la valeur de l'intersection qui se trouve dans la réponse de 
	la question 6 du mini-rapport.
	Lancer simplement le code sans modifer quoi que ce soit.
8. TP1_Q7.py
	- Permet de générer le graphique et la valeur de l'intersection qui se trouve dans la réponse de 
	la question 7 du mini-rapport.
	Lancer simplement le code sans modifer quoi que ce soit.