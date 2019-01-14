# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:34:38 2019

@author: Alexis
"""
import numpy as np
class Voiture():
    def __init__(self, nom, puissance, couple, poids, acceleration, prix, pollution):
        self.nom = nom
        self.puissance = puissance
        self.couple = couple
        self.poids = poids
        self.acceleration = acceleration
        self.prix = prix
        self.pollution = pollution
    
    def get_params(self):
        return  np.array([self.puissance, self.couple, self.poids, self.acceleration, self.prix, self.pollution])

def import_data(data):
    voitures = []
    for i in range(len(data)):
        voitures.append(Voiture(data[i][0], data[i][1], data[i][2], data[i][3], data[i][4], data[i][5], data[i][6]))
    return voitures

def pareto(voitures):
    pareto = [0 for i in range(len(voitures))]
    for i in range(len(voitures)):
        voiture = voitures[i]
        is_pareto = True
        for voiture2 in voitures:
            if sum([0 if voiture.get_params()[t] >= voiture2.get_params()[t] else 1 for t in range(len(voiture2.get_params()))]) == 0 and not np.array_equal(voiture.get_params(),voiture2.get_params()):
                is_pareto = False
                break
        
        pareto[i] = is_pareto
    
    return pareto

def somme_pondere(voitures, coefficients, minimum = True):
    index = 0
    valeur = float("inf") if minimum else -float("inf")
    for i in range(len(voitures)):
        for j in range(len(coefficients)):
            somme = sum(voitures[i].get_params() * coefficients)
            if minimum and somme <= valeur:
                valeur = somme
                index = i
            
            elif not minimum and somme>= valeur:
                valeur = somme
                index = i
        
    return index
                
            
def find_nadir(voitures):
    valeur = []
    for i in range(len(voitures[0].get_params())):
        valeur.append(voitures[somme_pondere(voitures, [1 if j == i else 0 for j in range(len(voitures[0].get_params()))], False)].get_params()[i])
    return np.array(valeur)

def find_ideal(voitures):
    valeur = []
    for i in range(len(voitures[0].get_params())):
        valeur.append(voitures[somme_pondere(voitures, [1 if j == i else 0 for j in range(len(voitures[0].get_params()))])].get_params()[i])
    return np.array(valeur)

def distanceTchebychev(x, x0, coeff):
    return max(np.abs(x - x0) * coeff)

def distanceTchebychevAugmentee(x, x0, coeff, epsilon=0.1):
    return distanceTchebychev(x, x0, coeff) + epsilon * np.sum(np.abs(x - x0) * coeff)

def distance(x, x1, x0, alpha, epsilon=0.1):
    coeff = np.divide(np.ones(len(x1)), abs(x1 - x0)) * alpha
    return distanceTchebychevAugmentee(x, x1, coeff, epsilon)

def solutiondistanceTchebychevAugmentee(voitures, condition=False, epsilon=0.1):
    alpha = np.ones(len(voitures[0].get_params()))
    
    ideal = find_ideal(voitures)
    nadir = find_nadir(voitures)
        
    min = float("inf")
    sol = voitures[0]
    for voiture in voitures:
        continuer = 1
        if condition:
            for i, max in condition:
                if voiture.get_params()[i] >= max:
                    continuer = 0
        if continuer:
            val = distance(voiture.get_params(), ideal, nadir, alpha, epsilon)
            if val <= min:
                min = val
                sol = voiture    
    return sol


data = [["alpha romeo mito veloce", -170, -250, 1145, 73, 25890, 124],
        ["audi s1 quattro", -231, -370, 1315, 58, 35590, 166],
        ["abarth 595", -145, -206, 1035, 78, 18600, 139],
        ["ford fiesta 1.0 ecoboost 140 ch black/red edition", -140, -210, 1128, 90, 18800, 104],
        ["ford fiesta st", -182, -240, 1163, 69, 24950, 138],
        ["mini cooper 5", -192, -280, 1160, 68, 25600, 134.5],
        ["mini cooper s jcw", -231, -320, 1205, 63, 32550, 155],
        ["opel adam s", -150, -220, 1103, 85, 19000, 139],
        ["opel corsa (e) opc", -207, -280, 1218, 68, 23900, 174],
        ["peugeot 208 gti", -208, -300, 1160, 65, 26000, 125],
        ["renault clio 4 rs", -200, -240, 1204, 67, 28100, 133],
        ["seat ibiza cupra 192 ch", -192, -320, 1269, 67, 22870, 139],
        ["suzuki swift sport", -136, -160, 1040, 87, 25680, 147],
        ["volkswagen polo bluegt", -150, -250, 1212, 79, 23620, 110],
        ["volkswagen polo gti", -192, -320, 1269, 67, 26900, 139]]

voitures = import_data(data)
index_pareto = pareto(voitures)
pareto = [voitures[i] for i in range(len(voitures)) if index_pareto[i]]

solution1 = solutiondistanceTchebychevAugmentee(pareto)
solution2 = solutiondistanceTchebychevAugmentee(pareto, [(0, -208), (1, -320)])

    