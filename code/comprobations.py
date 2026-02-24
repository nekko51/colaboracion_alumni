#File aimed to do the previous comprobations of a .txt file containings seqs
from constants import *

#Given a file, determine whether all its lines are the same length or not:
def SameLength(f):
    f.seek(0)

    lengths = []
    for line in f:
        lengths.append(len(line))

    # We'll just take the length of the first line and compare the rest to that first one
    n = lengths[0]
    for i in range(len(lengths)):
        if lengths[i] != n:
            errorMsg = f"Line {i}'s length is {lengths[i]}, which is different from the first one's length ({n})"
            return False, errorMsg
            #Possible upgrade: create a bool variable that turns to false when the first error arises, 
            # but let the for loop print out all lines that aren't the same length as the first one.
    
    return True, ""


#Given a file, determine whether there appear (aminoácidos ambigüos) or not:
def AmbiguousAminoacids(f):
    f.seek(0)

    cnt = 0
    for line in f:
        cnt += 1
        for letter in line:
            if letter == "B" or letter == "J" or letter == "X" or letter == "Z" or letter == "U" or letter == "O":
                errorMsg = f"There appears an (aminoácido ambigüo) in line {cnt} with illegal letter {letter}"
                return False, errorMsg
                #Possible upgrade: create a bool variable that turns to false when the first error arises, 
                # but let the for loop print out all illegal aminoacids (same idea as before).
    return True, ""


#Possible alternative for the analysis of amb aminoacids. It should be more efficient.
#Lo malo de este código es que no detecta la posición en la que se ha encontrado el aminoácido erróneo, pero detecta más rápido
def AmbiguousAminoacids2 (f):
    f.seek(0)

    valid = set(aminoacids.keys())
    valid.add ("\n")
    sequence = set(f.read())

    ambiguous = sequence - valid

    if ambiguous:
        ambiguousAminoacidsSuccesful, aminoacidsErrorMsg = AmbiguousAminoacids(f)
        if ambiguousAminoacidsSuccesful == True:
            print(f"Weird... Ambiguous Aminoacids function nº2 returned false but nº1 returned true")
        return ambiguousAminoacidsSuccesful, aminoacidsErrorMsg
    return True, ""