# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:34:38 2019

@author: Alexis
"""
import numpy as np
import random as r
import operator
import scipy.optimize
import madmcAlexis as aux
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
        ["volkswagen polo gti", -192, -320, 1269, 67, 26900, 139],
        ["ALFA ROMEO Giulietta Veloce", -240, -340, 1320, 60, 35.600, 162],
        ["AUDI S3", -310, -380, 1425, 52, 50000, 162],
        ["AUDI RS3", -367, -465, 1425, 43, 57900, 189],
        ["BMW M140i", -340, -500, 1450, 48, 45600, 179],
        ["FORD Focus ST", -250, -360, 1362, 65, 29200, 159],
        ["FORD Focus RS 350", -350, -440, 1524, 47, 39600, 175],
        ]

voitures = import_data(data)
index_pareto = pareto(voitures)
pareto = [voitures[i] for i in range(len(voitures)) if index_pareto[i]]

solution1 = solutiondistanceTchebychevAugmentee(pareto)
solution2 = solutiondistanceTchebychevAugmentee(pareto, [(0, -208), (1, -320)])



#### Partie 2 ###

def normalize(voitures, nadir, ideal):
    voitures_normal = []
    normaliseur = [nadir[i] - ideal[i] for i in range(6)]
    for i in range(len(voitures)):
        voiture = voitures[i].get_params()
        voitures_normal.append(
            [
                (voiture[0] - ideal[0] )/normaliseur[0], 
                (voiture[1] - ideal[1] )/normaliseur[1], 
                (voiture[2] - ideal[2] )/normaliseur[2], 
                (voiture[3] - ideal[3] )/normaliseur[3], 
                (voiture[4] - ideal[4] )/normaliseur[4], 
                (voiture[5] - ideal[5] )/normaliseur[5]
            ]
            
            )
    return voitures_normal

def query(x,y,w_DM):
    fx = [x[i] * w_DM[i] for i in range(len(x))]
    fy = [y[i] * w_DM[i] for i in range(len(y))]
    if (fx > fy) :
        return (x,y)
    else:
        return (y,x)


def PL_pmr(x,y,P):

    c = np.zeros(6) # objective vector 
    for i in range(6):
        c[i] = - x[i] + y[i]
    A_ub = np.zeros( ( len(P) , 6 ) ) # constraint array using P, i.e. answers of DM to queries
    for i in range(len(P)) :
        for j in range(6) :
            A_ub[i][j] = P[i][1][j] - P[i][0][j]
    b_ub = np.zeros(len(P))
    A_eq = np.ones((1,6))
    b_eq = np.ones(1)

    # print(A_ub)

    # bounds variable bounds 
    w0_bounds = (0,1)
    w1_bounds = (0,1)
    w2_bounds = (0,1)
    w3_bounds = (0,1)
    w4_bounds = (0,1)
    w5_bounds = (0,1)

    bounds = (w0_bounds, w1_bounds, w2_bounds, w3_bounds, w4_bounds, w5_bounds)
    
    res = scipy.optimize.linprog(c, A_ub, b_ub, A_eq, b_eq, bounds)
    return -res.fun

def inP(x,y,P):
    return P.count((x,y)) > 0

def neq_list(l1, l2):
    return len(l1) != len(l2) or any( l1[i] != l2[i] for i in range(len(l1)) )


def increment(voitures, eps, DMcoeff):
    # print(DMcoeff)
    P = []
    n_query = 0
    rank = 0
    min_pmr = float("inf")
    old_x_star_i = -1
    x_star = None
    y_star = None
    argmin = []
    argmax_tab = []
    y_max = None
    while min_pmr >= eps :
        # print("--- new iteration ---")
        for x in voitures:
            argmax = []
            max_pmr = -float("inf")
            for y in voitures:
                if (neq_list(x,y)):
                    pmr = PL_pmr(x,y,P)
                    if (pmr > max_pmr):
                        max_pmr = pmr
                        # y_max = y
                        argmax = [y]
                    elif (pmr == max_pmr and not inP(x, y, P) ): # avoid chosing the same y* every time
                    #     y_max = y
                        argmax.append(y)
            if (max_pmr < min_pmr):
                min_pmr = max_pmr
                # x_star = x
                # y_star = y_max
                argmin = [x]
                argmax_tab = [argmax]
            elif (max_pmr == min_pmr and len(argmax) > 0 ):
            #     x_star = x
            #     y_star = y_max
                argmin.append(x)
                argmax_tab.append(argmax)
        x_star_i = r.randint(0, len(argmin)-1)
        if (old_x_star_i != x_star_i):
            rank = n_query
        x_star = argmin[x_star_i]
        argmax = argmax_tab[x_star_i]
        y_star = argmax[r.randint(0, len(argmax)-1)]
        P.append(query(x_star,y_star,DMcoeff))
        n_query += 1
        # print(aux.data[voitures.index(x_star)][0])
    #     print(min_pmr)
    # print("------")
    # print(aux.data[voitures.index(x_star)][0])
    # print('nb questions')
    # print(n_query)    
    # print(P)
    # print("------")
    return (n_query,rank,voitures.index(x_star))

def ideal(voitures):
    ideal = [float("inf"),float("inf"),float("inf"),float("inf"),float("inf"),float("inf")]
    for i in range(len(voitures)):
        for j in range(6):
            if (voitures[i].get_params()[j] < ideal[j]):
                ideal[j] = voitures[i].get_params()[j]
    return ideal

def nadir(pareto):
    nadir = [-float("inf"),-float("inf"),-float("inf"),-float("inf"),-float("inf"),-float("inf")]
    for i in range(len(pareto)):
        for j in range(6):
            if (pareto[i].get_params()[j] > nadir[j]):
                nadir[j] = pareto[i].get_params()[j]
    return nadir

def randomDM():
    vect = [ r.random() for i in range(6) ]
    s = sum(vect)
    return [ vi / s for vi in vect ]

# # DMcoeff = [.2,.1,.1,.2,.3,.1]
# # DMcoeff = [.0,.0,.5,.25,.0,.25]
# DMcoeff = randomDM()
# # pareto = [ voitures[i] for i in range(len(voitures)) if index_pareto[i]]
# voitures = import_data(aux.data)
# ideal = ideal(voitures)
# nadir = nadir(voitures)
# print(nadir)
# print(ideal)
# voitures = normalize(voitures, nadir, ideal)

# increment(voitures, 0.001)


def run():
    DMcoeff = randomDM()
    voitures = import_data(aux.data)
    ideal_list = ideal(voitures)
    nadir_list = nadir(voitures)
    voitures = normalize(voitures, nadir_list, ideal_list)

    return increment(voitures, 0.001,DMcoeff)

def runs(n):
    results = []
    ranks = []
    indexes = []
    for i in range(n):
        res = run()
        print(res)
        results.append(res[0])
        ranks.append(res[1])
        indexes.append(res[2])
    print(results)
    print(ranks)
    print(indexes)

runs(10)
