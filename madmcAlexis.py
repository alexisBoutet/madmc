# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:34:38 2019

@author: Alexis
"""
import numpy as np
import random

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

def pareto(liste):
    pareto = [0 for i in range(len(liste))]
    for i in range(len(liste)):
        valeur = liste[i]
        is_pareto = True
        for valeur2 in liste:
            if sum([0 if valeur[t] >= valeur2[t] else 1 for t in range(len(valeur2))]) == 0 and not np.array_equal(valeur,valeur2):
                is_pareto = False
                break
        
        pareto[i] = is_pareto
    
    return pareto

def somme_pondere(liste, coefficients, minimum = True):
    index = 0
    valeur = float("inf") if minimum else -float("inf")
    for i in range(len(liste)):
        somme = sum(liste[i] * coefficients)
        if minimum and somme <= valeur:
            valeur = somme
            index = i
        
        elif not minimum and somme>= valeur:
            valeur = somme
            index = i
        
    return index
                
            
def find_nadir(liste):
    valeur = []
    for i in range(len(liste[0])):
        valeur.append(liste[somme_pondere(liste, [1 if j == i else 0 for j in range(len(liste[0]))], False)][i])
    return np.array(valeur)

def find_ideal(liste):
    valeur = []
    for i in range(len(liste[0])):
        valeur.append(liste[somme_pondere(liste, [1 if j == i else 0 for j in range(len(liste[0]))])][i])
    return np.array(valeur)

def distanceTchebychev(x, x0, coeff):
    return max(np.abs(x - x0) * coeff)

def distanceTchebychevAugmentee(x, x0, coeff, epsilon=0.1):
    return distanceTchebychev(x, x0, coeff) + epsilon * np.sum(np.abs(x - x0) * coeff)

def distance(x, x1, x0, alpha, epsilon=0.1):
    coeff = np.divide(np.ones(len(x1)), abs(x1 - x0)) * alpha
    return distanceTchebychevAugmentee(x, x1, coeff, epsilon)

def solutiondistanceTchebychevAugmentee(liste, condition=False, epsilon=0.1):
    alpha = np.ones(len(liste[0]))
    
    ideal = find_ideal(liste)
    nadir = find_nadir(liste)
        
    min = float("inf")
    sol = 0
    for i in range(len(liste)):
        valeur = liste[i]
        continuer = 1
        if condition:
            for j, max in condition:
                if valeur[j] >= max:
                    continuer = 0
        if continuer:
            val = distance(valeur, ideal, nadir, alpha, epsilon)
            if val <= min:
                min = val
                sol = i    
    return sol

def randomCoeff(voitures):
    p = len(voitures[0].get_params())
    coeff = [random.gauss(1 / np.sum([abs(voiture.get_params()[k]) for voiture in voitures]),2) for k in range(p)]
    return np.array([c if c > 0 else 0 for c in coeff])/np.sum(coeff)
    
def interaction(voitures):
    n, p = len(voitures[0].get_params()), len(voitures)
    index_pareto = pareto([v.get_params() for v in voitures])
    voitures = [voitures[i] for i in range(p) if index_pareto[i]]
    liste = [voiture.get_params() for voiture in voitures]
    ideal = find_ideal(liste)
    nadir = find_nadir(liste)    
    
    coefficient = randomCoeff(voitures)
    print(coefficient)
    
    i = somme_pondere(liste, coefficient)

    condition = []
    continuer = 1
    nb = 0
    while continuer:
        nb+=1
        j = solutiondistanceTchebychevAugmentee(liste, condition)
        print(i,j)
        if j==i:
            print(nb)
            return voitures[i]
        else:
            indexMaxEcart = np.argmax(np.array([(voitures[j].get_params()[k] -voitures[i].get_params()[k])/(nadir[k]-ideal[k]) for k in range(n)]))
            condition.append((indexMaxEcart, voitures[j].get_params()[indexMaxEcart]))
    
    
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
        ["volkswagen polo gti", -192, -320, 1269, 67, 26900, 139],
        ["alfa romeo giuluetta veloce", -240, -340, 1320, 60, 35600, 162],
        ["audi d3", -310, -380, 1425, 52, 50000, 162],
        ["audi rs3", -367, -465, 1425, 43, 57900, 189],
        ["bmw m140i", -340, -500, 1450, 48, 45600, 179],
        ["ford focus st", -250, -360, 1362, 65, 29200, 159],
        ["ford focus rs 350", -350, -440, 1524, 47, 39600, 175],
        ["kia pro_cee'd gt", -204, -265, 1359, 77, 29290, 170],
        ["mercedes a45 amg", -381, -475, 1480, 42, 55350, 162],
        ["mini clubman couper s (f54)", -192, -280, 1360, 72, 29800, 147],
        ["peugeot 308 gt 205", -205, -285, 1200, 75, 31350, 130],
        ["peugeot 308 gti 270", -270, -330, 1205, 60, 37700, 139],
        ["renault megane 4 gt", -205, -280, 1392, 71, 32200, 134],
        ["seat leon fr", -184, -380, 1295, 75, 29720, 113],
        ["seat leon cupra", -300, -380, 1375, 57, 35000, 149],
        ["volkswagen golf gtd", -184, -380, 1475, 79, 36570, 115],
        ["volkswagen golf gte", -204, -350, 1524, 76, 38900, 35],
        ["volkswagen golf gte", -237.5, -380, 1307, 64, 34530, 147],
        ["volkswagen golf r", -310, -380, 1401, 51, 43980, 180],
        ["aston martin rapide", -470, -600, 1990, 181000, 53, 300],
        ["audi 54", -330, -440, 1780, 62740, 50, 190],
        ["bmw m5 (f10)", -560, -680, 1875, 117100, 43, 235],
        ["cadillac cts-v", -564, -747, 1742, 68127, 41, 365],
        ["mercedes c63 amg", -457, -600, 1730, 89300, 44, 280],
        ["mitsubishi lancer evo x", -295, -366, 1600, 52850, 66, 250],
        ["opel insignia opc", -325, -435, 1750, 45350, 60, 249],
        ["porsche panamera gts", -430, -520, 1920, 119575, 45, 256],
        ["subaru wrx sti", -300, -407, 1507, 44950, 52, 242]]

voitures = import_data(data)
#i = interaction(voitures)


