# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 18:02:11 2019

@author: Alexis
"""

import numpy as np
import random
import pulp
from itertools import product
from time import time

from madmc import *

class Objet():
    def __init__(self, poids, utilites):
        self.poids = poids
        self.utilites = utilites
    
    def print(self):
        print("utilite", self.utilites)
        print("poids", self.poids)


def generateSac(n, p):
    objets = [Objet(random.randint(1, 10), [random.randint(1, 10) for i in range(n)]) for j in range(p)]
    sac = sum([objet.poids for objet in objets])/2
    return objets, sac

def randomCoeff(p, i=1, j=100):
    coeff = [random.randint(i, j) for k in range(p)]
    return np.array([c if c > 0 else 0 for c in coeff])/np.sum(coeff)


def sommePondereSac(objets, poidsMax, coeff):
    n, p = len(objets[0].utilites), len(objets)
    model = pulp.LpProblem("Cost minimising blending problem", pulp.LpMaximize)
    
    variable = []
    variable += [pulp.LpVariable('x'+str(i), cat='Binary') for i in range(p)]

    model += pulp.lpSum([coeff[j] * objets[i].utilites[j] * variable[i] for i, j in product(range(p), range(n))])
    model += pulp.lpSum([objets[i].poids * variable[i] for i in range(p)]) <= poidsMax

    model.solve()
    return model
    
def idealNadirSac(objets, poidsMax, condition=False, epsilon=0.1):
    n, p = len(objets[0].utilites), len(objets)
    values = []
    for j in range(n):
        
        model = pulp.LpProblem("Cost minimising blending problem", pulp.LpMaximize)

        variable = []
        variable += [pulp.LpVariable('x'+str(i), cat='Binary') for i in range(p)]
        model += pulp.lpSum([objets[i].utilites[j] * variable[i] for i in range(p)])
        model += pulp.lpSum([objets[i].poids * variable[i] for i in range(p)]) <= poidsMax
        
        if condition:
            for i, max in condition:
                model += pulp.lpSum([objets[i].utilites[i] * variable[j] for j in range(p)]) >= max
            
        model.solve()
        values.append([sum([variable[j].varValue * objets[j].utilites[i] for j in range(p)]) for i in range(n)])
    if len(values[0]) == 1:
        ideal = [v[0] for v in values]
        nadir = [0 for v in values]
    else:
        ideal = [np.max(v) for v in values]
        nadir = [np.min([values[i][j] for j in range(n)]) for i in range(n)]
    return ideal, nadir
        
    
def resolutionPLSac(objets, poidsMax, condition=False, epsilon=0.1):
    n, p = len(objets[0].utilites), len(objets)
    
    ideal, nadir = idealNadirSac(objets, poidsMax, condition=False)    
    
    model = pulp.LpProblem("Cost minimising blending problem", pulp.LpMinimize)
    
    variable = []
    variable.append(pulp.LpVariable('z', cat='Continuous'))
    variable += [pulp.LpVariable('x'+str(i), cat='Binary') for i in range(p)]
    
    model += 1 * variable[0], "m"
    
    # Contrainte de sac à dos
    model += pulp.lpSum([objets[i].poids * variable[i+1] for i in range(p)]) <= poidsMax
     
    for i in range(n):        
        model += pulp.lpSum([-((nadir[i]-objets[j].utilites[i])/(nadir[i] - ideal[i] +1e-5) + epsilon * sum([(nadir[k]-objets[j].utilites[k])/(nadir[k] - ideal[k] +1e-5) for k in range(n)])) * variable[j+1] for j in range(p)]) <= variable[0]
    
    if condition:
        for i, max in condition:
            model += pulp.lpSum([objets[j].utilites[i] * variable[j+1] for j in range(p)]) >= max + 1
    t = time()
    model.solve()
    
    return model, time() - t

def interaction(objets, poidsMax):
    n, p = len(objets[0].utilites), len(objets)
    
    ideal, nadir = idealNadirSac(objets, poidsMax, condition=False)
    coeff = randomCoeff(n)
    
    m = sommePondereSac(objets, poidsMax, coeff)
    sacCible = [sum([objets[i].utilites[j] * m.variables()[i].varValue for i in range(p)])for j in range(n)]
    condition = []
    continuer = 1
    t = 0
    while continuer:
        t +=1
        model = resolutionPLSac(objets, poidsMax, condition)
        sacUtilites = [sum([objets[i].utilites[j] * model.variables()[i].varValue for i in range(p)])for j in range(n)]

        if sacCible==sacUtilites:
            return model, t
        else:
            indexMaxEcart = np.argmax(np.array([(sacCible[k] - sacUtilites[k])/(ideal[k] - nadir[k] +1e-5) for k in range(n)]))
            condition.append((indexMaxEcart, sacUtilites[indexMaxEcart]))


    
couples = [(5, 10), (5,20), (10,50), (10, 100), (10, 500)]
def test(i, j):
    for k in range(i):
        t = time()
        temps = 0
        
        for l in range(200):
            objets, poidsMax = generateSac(k+1, j)
            r, tps = resolutionPLSac(objets, poidsMax)
            temps += tps
        print("Test sur "+str(k+1)+" critères et "+str(j)+" objets")
        print(time()-t, temps)
        print()
        
for i, j in couples:
    test(i, j)