#Imports
import numpy as np
import matplotlib.pyplot as plt

#Global variables and aminoacid dictionary
CHAIN_LENGTH = 298
AMINOACID_NUMBER = 21
aminoacids = {
    "-": 0,
    "A": 1,
    "C": 2,
    "D": 3,
    "E": 4,
    "F": 5,
    "G": 6,
    "H" : 7,
    "I": 8,
    "K": 9,
    "L": 10,
    "M": 11,
    "N": 12,
    "P": 13,
    "Q": 14,
    "R": 15,
    "S": 16,
    "T": 17,
    "V": 18,
    "W": 19,
    "Y": 20,
    #B, J, X, Z, U y O no las añado de momento
}

reverse_aminoacids = {
    0: '-',
    1: 'A',
    2: 'C',
    3: 'D',
    4: 'E',
    5: 'F',
    6: 'G',
    7: 'H',
    8: 'I',
    9: 'K',
    10: 'L',
    11: 'M',
    12: 'N',
    13: 'P',
    14: 'Q',
    15: 'R',
    16: 'S',
    17: 'T',
    18: 'V',
    19: 'W',
    20: 'Y'
}

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


#Given a filename, run all previous checks before starting to work on it:
def PreviousComprobations(filename):
    #Try to open file:
    try:
        f = open(filename, "r")
    except FileNotFoundError:
        print(f"Couldn't open {filename}")
        return False

    #Run previous checks:
    sameLengthSuccesful, lengthErrorMsg = SameLength(f)
    if sameLengthSuccesful == False:
        print(lengthErrorMsg)
        return False
    
    ambiguousAminoacidsSuccesful, aminoacidsErrorMsg = AmbiguousAminoacids(f)
    if ambiguousAminoacidsSuccesful == False:
        print(aminoacidsErrorMsg)
        return False
    
    return True


#Calculate each aminoacid's frequency in every position
def AminoacidFrequency(filename):
    #Try to open file:
    try:
        f = open(filename, "r")
    except FileNotFoundError:
        print(f"Couldn't open {filename}")
        return False
    
    appearance = np.zeros((CHAIN_LENGTH, AMINOACID_NUMBER), dtype = int)
    total_chains = 0

    for line in f:          #Would prefer ranges (for i in range)
        line = line.strip()
        total_chains += 1
        pos = 0
        for letter in line:
            appearance[pos, aminoacids.get(letter)] += 1
            pos += 1

    #Calculate frequencies
    frequency = appearance/total_chains
    return appearance, frequency






# filename = "testing_adulterated1.seqs"
filename = "testing_human.seqs"
if PreviousComprobations(filename) == False:
    print(f"There was an error, stopping the program...")
else:
    appearance, frequency = AminoacidFrequency(filename)
    for i in range(CHAIN_LENGTH):
        sum = 0
        for j in range(AMINOACID_NUMBER):
            sum += frequency[i, j]
        print(sum)
    
    #We could make the color map MEAN something; that is, give all aminoacids that are hydrophobic a red-ish color, etc. 
    # - that way we could uncover "hidden" patterns relating to the chemical properties of each aminoacid in our chains
    x = range(CHAIN_LENGTH)
    random_color = np.zeros((AMINOACID_NUMBER,3))
    for i in range(AMINOACID_NUMBER):
        random_color[i] = np.random.rand(3,)
        plt.plot(range(CHAIN_LENGTH), frequency[:, i], 
                 marker=".", ms = 4, color=random_color[i],
                 linestyle = "solid",
                 label = reverse_aminoacids.get(i))
    
    #Plot customization:
    plt.title("Human Chains", fontsize=18, color ="black")
    plt.xlabel("Chain Position", fontsize = 12)
    plt.ylabel("Aminoacid Frequency", fontsize = 12)
    plt.xticks(range(0,301,25))
    plt.legend(
    bbox_to_anchor=(1.05, 1),
    loc='best',
    fontsize=4,
    ncol=2,
    title="Aminoacids"
    )
    plt.tight_layout()

    # plt.show()
    plt.savefig("human.png", dpi = 600)
    plt.close()

filename = "testing_mouse.seqs"
if PreviousComprobations(filename) == False:
    print(f"There was an error, stopping the program...")
else:
    appearance, frequency = AminoacidFrequency(filename)
    for i in range(CHAIN_LENGTH):
        sum = 0
        for j in range(AMINOACID_NUMBER):
            sum += frequency[i, j]
        print(sum)
    
    #We could make the color map MEAN something; that is, give all aminoacids that are hydrophobic a red-ish color, etc. 
    # - that way we could uncover "hidden" patterns relating to the chemical properties of each aminoacid in our chains
    x = range(CHAIN_LENGTH)
    for i in range(AMINOACID_NUMBER):
        plt.plot(range(CHAIN_LENGTH), frequency[:, i], 
                 marker=".", ms = 4, color=random_color[i],
                 linestyle = "solid",
                 label = reverse_aminoacids.get(i))
    
    #Plot customization:
    plt.title("Mouse Chains", fontsize=18, color ="black")
    plt.xlabel("Chain Position", fontsize = 12)
    plt.ylabel("Aminoacid Frequency", fontsize = 12)
    plt.xticks(range(0,301,25))
    plt.legend(
    bbox_to_anchor=(1.05, 1),
    loc='best',
    fontsize=4,
    ncol=2,
    title="Aminoacids"
    )
    plt.tight_layout()

    # plt.show()
    plt.savefig("mouse.png", dpi = 600)
    plt.close()