#Imports
import numpy as np
import matplotlib.pyplot as plt
from constants import *


#Functions that calculate entropies in different ways:
#Given a frequency array, alpha value (the summation index for each S_i, in AA's it's 21), 
# return a numpy array with the entropies S_i
def shannon_entropy (probabilities_matrix):
    # S_i = - sum(p_i * log(p_i))
    # probabilities is matrix where len(probabilities_matrix) = number of positions
    S_i = np.zeros(len(probabilities_matrix))
    for i in range(len(probabilities_matrix)):
        for prob in probabilities_matrix[i]:
            if prob > 0:
                S_i[i] = S_i[i] - (prob * np.log(prob))
    return S_i

#Given an entropy array, do things with it (other functions):

def lineal_entropy (probabilities_matrix):
    # L_i = - sum(p_i * (1 - p_i))
    L_i = np.zeros(len(probabilities_matrix))
    for i in range(len(probabilities_matrix)):
        for prob in probabilities_matrix[i]:
            L_i[i] = L_i[i] - (prob*(1 - prob))
    return L_i

def renyi_entropy (probabilities_matrix, q):
    # R_q = log(sum(p_i**q))/(1-q)
    # S = lim (R_q) con q -> 1
    R_i = np.zeros(len(probabilities_matrix))
    if q >= 0 and q != 1:
        for i in range(len(probabilities_matrix)):
            for prob in probabilities_matrix[i]:
                if prob > 0:
                    R_i[i] = R_i[i] + prob** q
                    # If to avoid 0**0 = 1                    
            if R_i[i]>0:
                R_i[i] = (1/(1-q))*np.log(R_i[i])
            else: R_i[i]= -1
    return R_i

def tsallis_entropy (probabilities_matrix, q):
    # T_q = (1-sum(p_i**q))/(q-1)
    T_i = np.zeros(len(probabilities_matrix))
    if q >= 0 and q != 1:
        for i in range(len(probabilities_matrix)):
            for prob in probabilities_matrix[i]:
                if prob > 0:
                    T_i[i] = T_i[i] + prob** q
            T_i[i] = (1/(q-1))*(1-T_i[i])
    return T_i



#Selector functions that, given a "word code", calculates a kind of entropy -- or f*** it, 
# calculate every single one could be an option too.