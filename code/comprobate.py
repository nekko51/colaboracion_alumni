import comprobations as comp
from constants import *

#Given a filename, run previous checks and add it to "comprobated.txt":
def PreviousComprobations(filename):
    #Try to open file:
    try:
        f = open(filename, "r")
    except FileNotFoundError:
        print(f"Couldn't open {filename}")
        return False

    #Run previous checks:
    sameLengthSuccesful, lengthErrorMsg = comp.SameLength(f)
    if sameLengthSuccesful == False:
        print(lengthErrorMsg)
        return False
    
    ambiguousAminoacidsSuccesful, aminoacidsErrorMsg = comp.AmbiguousAminoacids2(f)
    if ambiguousAminoacidsSuccesful == False:
        print(aminoacidsErrorMsg)
        return False
    
    return True


def Comprobate(filename):
    if PreviousComprobations(filename) == True:
        #Try to open file:
        try:
            f = open("comprobations/comprobated.txt", "r")
        except FileNotFoundError:
            print(f"Couldn't open comprobations/comprobated.txt in r mode")
            return False
        
        already_present = False
        for line in f:
            line = line.strip()
            if line == filename:
                already_present = True
        f.close()

        try:
            f = open("comprobations/comprobated.txt", "a")
        except FileNotFoundError:
            print(f"Couldn't open comprobations/comprobated.txt in a mode")
            return False

        if already_present == False:
            f.write(f"{filename}\n")
            print(f"Added {filename} to comprobated directory")
        else:
            print(f"Didn't add {filename} to directory since it's already present")

Comprobate("seqs/learn_human.txt")
Comprobate("seqs/learn_mouse.txt")
Comprobate("seqs/test_human.txt")
Comprobate("seqs/test_mouse.txt")