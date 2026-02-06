#Imports
import numpy as np
import matplotlib.pyplot as plt
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
    
    ambiguousAminoacidsSuccesful, aminoacidsErrorMsg = AmbiguousAminoacids2(f)
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


#Main function; checks whether a filename contains chains of aminoacids (of the same length) without
# any ambiguous aminoacids, and generates a plot for the frequency of each aminoacid in each position.
def GeneratePlot(filename, title, filename_output, colors):
    if PreviousComprobations(filename) == False:
        print(f"There was an error, stopping plot for {filename}")
        return
    appearance, frequency = AminoacidFrequency(filename)

    #We could make the color map MEAN something; that is, give all aminoacids that are hydrophobic a red-ish color, etc. 
    # - that way we could uncover "hidden" patterns relating to the chemical properties of each aminoacid in our chains
    
    #Plot generation:
    plt.figure(figsize=(16, 8))
    x = range(CHAIN_LENGTH)
    for i in range(AMINOACID_NUMBER):
        plt.plot(x, frequency[:, i], 
                    marker=".", ms = 4, color=colors[i],
                    linestyle = "solid", alpha = 0.5,
                    label = reverse_aminoacids.get(i))
    
    #Plot customization:
    plt.title(title, fontsize=18, weight='bold')
    plt.xlabel("Chain Position", fontsize = 12)
    plt.ylabel("Aminoacid Frequency", fontsize = 12)
    plt.xticks(range(0, CHAIN_LENGTH+1, 25))
    #Legend:
    plt.legend(
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        fontsize=15,
        ncol=1,
        title="Aminoacids",
        title_fontsize = 10
        )
    #Saving
    plt.tight_layout()
    plt.savefig(filename_output, dpi = 480)
    plt.close()
    print(f"Image {filename_output} created")


def GenerateComparativePlot(filename1, filename2, title, filename_output, colors):
    if PreviousComprobations(filename1) == False:
        print(f"There was an error, stopping plot for {filename1}")
        return
    appearance1, frequency1 = AminoacidFrequency(filename1)
    if PreviousComprobations(filename2) == False:
        print(f"There was an error, stopping plot for {filename2}")
        return
    appearance2, frequency2 = AminoacidFrequency(filename2)

    frequency = np.absolute(frequency1-frequency2)
    #Plot generation:
    plt.figure(figsize=(16, 8))
    x = range(CHAIN_LENGTH)
    for i in range(AMINOACID_NUMBER):
        plt.plot(x, frequency[:, i], 
                    marker=".", ms = 4, color=colors[i],
                    linestyle = "solid", alpha = 0.5,
                    label = reverse_aminoacids.get(i))
    
    #Plot customization:
    plt.title(title, fontsize=18, weight='bold')
    plt.xlabel("Chain Position", fontsize = 12)
    plt.ylabel("Aminoacid Frequency", fontsize = 12)
    plt.xticks(range(0, CHAIN_LENGTH+1, 25))
    #Legend:
    plt.legend(
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        fontsize=15,
        ncol=1,
        title="Aminoacids",
        title_fontsize = 10
        )
    #Saving
    plt.tight_layout()
    plt.savefig(filename_output, dpi = 1200)
    plt.close()
    print(f"Image {filename_output} created")

#Colormap importation - pay attention to the structure in "colormap.txt" when changing it; or rewrite this section entirely. PAY ATTENTION MARTA >:(
colormap = np.loadtxt("colormap/colormap.txt")

GeneratePlot("seqs/testing_human.seqs", "Human Chains", "images/human.png", colormap)
GeneratePlot("seqs/testing_mouse.seqs", "Mouse Chains", "images/mouse.png", colormap)
GenerateComparativePlot("seqs/testing_human.seqs", "seqs/testing_mouse.seqs", "Relative Frequency", "images/relative.png", colormap)

#   Next steps:
#     LP - Low priority; MP - Medium Priority; HP - High Priority; UHP - Critical Priority
#       [VISUAL - LP] Make it so the name of the aminoacid shows up in its correspondig point on the graph 
#           if its frequency is higher than some value (e.g. 0.5)
#       [HP] Create a fixed, meaningful colormap (for example color corresponding to chem properties)
#       [DATA - MP] Maybe make sub-plots (every 25-50 positions, make a new graph) to compare more easily
#       [HP] Come up with more ways to compare the two plots.
#       [TOYING AROUND - LP] Try poking around and tweaking the sequences so that there exists an offset between the data 
#           (don't know what that'd be good for, but I'll just throw it in the list i guess)
#       [MP] Maybe rethink the functions so that we can work with the "frequency" and "appearance" arrays
#           (e.g. to generate various plots with different colormaps efficiently, without calculating them for every plot)