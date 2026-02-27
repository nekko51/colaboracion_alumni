#Code relating to the study of each AA's properties
from constants import *
import numpy as np
#I would take an approach similar to colormapping (and entropy, for that matter) -- create various ways to divvy up the AA's

# ==========
# WAYS TO DIVIDE
# ==========
# Remember the dictionaries below won't be a bijection so, obviously, we can't make a "reverse" one like we did with AA's and numbers.
# There are two cysteine possibilites (CSS and CSH; or C55 and C5H), I've defaulted to 0 when they aren't in the same group
# Maybe it'd be a good idea to include both posibilites, one where we default to 0 and another one where we default to 1
# They come into conflict in "minuscule" and "polar" categories


# Hidrophobic
hydrophobic = {
    "-": 0,
    "A": 1,
    "C": 1,
    "D": 0,
    "E": 0,
    "F": 1,
    "G": 1,
    "H": 1,
    "I": 1,
    "K": 1,
    "L": 1,
    "M": 1,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 0,
    "T": 1,
    "V": 1,
    "W": 1,
    "Y": 1,
}

# Aromatic
aromatic = {
    "-": 0,
    "A": 0,
    "C": 0,
    "D": 0,
    "E": 0,
    "F": 1,
    "G": 0,
    "H": 1,
    "I": 0,
    "K": 0,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 0,
    "T": 0,
    "V": 0,
    "W": 1,
    "Y": 1,
}

# Aliphatic
aliphatic = {
    "-": 0,
    "A": 0,
    "C": 0,
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 1,
    "K": 0,
    "L": 1,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 0,
    "T": 0,
    "V": 1,
    "W": 0,
    "Y": 0,
}

# Polar
polar = {
    "-": 0,
    "A": 0,
    "C": 0, #C5H is and C55 isn't, so I'm defaulting to 0 just in case? (or S instead of 5)
    "D": 1,
    "E": 1,
    "F": 0,
    "G": 0,
    "H": 1,
    "I": 0,
    "K": 1,
    "L": 0,
    "M": 0,
    "N": 1,
    "P": 0,
    "Q": 1,
    "R": 1,
    "S": 1,
    "T": 1,
    "V": 0,
    "W": 1,
    "Y": 1,
}

# Small
small = {
    "-": 0,
    "A": 1,
    "C": 1, #Both C55 and C5H are small
    "D": 1,
    "E": 0,
    "F": 0,
    "G": 1,
    "H": 0,
    "I": 0,
    "K": 0,
    "L": 0,
    "M": 0,
    "N": 1,
    "P": 1,
    "Q": 0,
    "R": 0,
    "S": 1,
    "T": 1,
    "V": 1,
    "W": 0,
    "Y": 0,
}

# Minuscule
minuscule = {
    "-": 0,
    "A": 1,
    "C": 0, #defaulting to 0
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 1,
    "H": 0,
    "I": 0,
    "K": 0,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 1,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0,
}

# +1 for positively charged; -1 for negatively charged
charged = {
    "-": 0,
    "A": 0,
    "C": 0,
    "D": -1,
    "E": -1,
    "F": 0,
    "G": 0,
    "H": +1,
    "I": 0,
    "K": +1,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": +1,
    "S": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0,
}

empty = {
    "-": 0,
    "A": 0,
    "C": 0,
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 0,
    "K": 0,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0,
}

valid_divisions = ["hydrophobic", "aromatic", "aliphatic", "polar", "small", "minuscule", "charged", "empty"]

# ==========
# FUNCTIONS
# ==========
#Given a frequency array, and a string representing division kind (hydrophobic, aromatic, aliphatic, polar, small, miniscule, charged), 
# return another frequency array according to that segmentation
def GenerateFrequencyArray(old_frequency, division):
    division.lower()
    #Check if the specified division is valid
    if division not in valid_divisions:
        print(f"{division} not valid, assigning empty...")
        division = "empty"

    if division == "charged":
        new_frequency = np.zeros((CHAIN_LENGTH, 3))
    elif not division == "charged":
        new_frequency = np.zeros((CHAIN_LENGTH, 2))
    
    for i in range(AMINOACID_NUMBER):
        if eval(division).get(reverse_aminoacids.get(i)) == 1:
            new_frequency[:,1] = old_frequency[:,i]
        elif eval(division).get(reverse_aminoacids.get(i)) == 0:
            new_frequency[:,0] = old_frequency[:,i]
    
    return(new_frequency)
