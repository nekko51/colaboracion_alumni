#Imports
import numpy as np
import matplotlib.pyplot as plt
from constants import *
import colormap as cm

#Given a filename, load valid, comprobated .txt files containing seqs
def LoadValid():
    valid = np.loadtxt("comprobations/comprobated.txt", dtype=str)
    return(valid)

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
    return(frequency) #We could also return "appearance" if needed

#Graph, given a frequency array, the frequencies with the specified filename output and dpi
def GeneratePlot(title, filename_output, colormap, frequency, dpi):
    #Plot generation:
    plt.figure(figsize=(16, 8))
    x = range(1, CHAIN_LENGTH+1)
    for i in range(AMINOACID_NUMBER):
        plt.plot(x, frequency[:, i], 
                    marker=".", ms = 4, color=colormap[i],
                    linestyle = "solid", alpha = 0.5,
                    label = reverse_aminoacids.get(i))
    
    #Plot customization:
    plt.title(title, fontsize=18, weight='bold')
    plt.xlabel("Chain Position", fontsize = 12)
    plt.ylabel("Aminoacid Frequency", fontsize = 12)
    plt.ylim(-0.05, 1.05)
    plt.xticks(range(1, CHAIN_LENGTH+1, 25))
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
    plt.savefig(filename_output, dpi = dpi)
    plt.close()
    print(f"Image {filename_output} created")

#Graph, given a frequency array, the frequencies of each aminoacid on its own (and join them all in a single graph)
def GenerateMiniPlot(title, filename_output, colormap, frequency, dpi):
    fig, axs = plt.subplots(3,7, sharey=True)
    fig.suptitle(title)
    axes = axs.flatten()
    x = range(1, CHAIN_LENGTH+1)

    for i in range(0, AMINOACID_NUMBER):
        ax = axes[i]
        ax.bar(x, frequency[:, i], color = colormap[i])
        ax.set_title(reverse_aminoacids.get(i))
    
    plt.tight_layout()
    plt.savefig(filename_output, dpi = dpi)
    plt.close
    print(f"Image {filename_output} created")

#Main function; checks whether a filename contains chains of aminoacids (of the same length) without
# any ambiguous aminoacids, and calls the graphing function.
def GenerateImage(filename, title, filename_output, colormap, dpi, valid):
    if filename not in valid:
        print(f"{filename} was not found in the \"comprobated.txt\" directory, stopping plot generation...")
        return
    frequency = AminoacidFrequency(filename)
    GeneratePlot(title, filename_output, colormap, frequency, dpi)

def GenerateComparativeImage(filename1, filename2, title, filename_output, comparative_filename_output, colormap, dpi, valid):
    if filename1 not in valid:
        print(f"{filename1} was not found in the \"comprobated.txt\" directory, stopping plot generation...")
        return
    frequency1 = AminoacidFrequency(filename1)
    if filename2 not in valid:
        print(f"{filename2} was not found in the \"comprobated.txt\" directory, stopping plot generation...")
        return
    frequency2 = AminoacidFrequency(filename2)

    frequency = np.absolute(frequency1-frequency2)

    GeneratePlot(title, filename_output, colormap, frequency, dpi)
    GenerateMiniPlot(title, comparative_filename_output, colormap, frequency, dpi)

#Colormap selection
colormap = cm.ColormapSelection("-ACDEFGHIKLMNPQRSTVWY", "rasmol")

#Load valid files
valid = LoadValid()
print(valid)

GenerateImage("seqs/learn_human.txt", "Human Chains", "images/learn_human.png", colormap, 480, valid)
GenerateImage("seqs/learn_mouse.txt", "Mouse Chains", "images/learn_mouse.png", colormap, 480, valid)
GenerateComparativeImage("seqs/learn_human.txt", "seqs/learn_mouse.txt", "Relative Frequency", 
                         "images/learn_relative.png", "images/learn_z_mini.png", colormap, 1200, valid)
GenerateImage("seqs/test_human.txt", "Human Chains", "images/test_human.png", colormap, 480, valid)
GenerateImage("seqs/test_mouse.txt", "Mouse Chains", "images/test_mouse.png", colormap, 480, valid)
GenerateComparativeImage("seqs/test_human.txt", "seqs/test_mouse.txt", "Relative Frequency", 
                         "images/test_relative.png", "images/test_z_mini.png", colormap, 1200, valid)
print(f"\nDone.")

#   Next steps:
#     LP - Low priority; MP - Medium Priority; HP - High Priority; UHP - Critical Priority
#       [VISUAL - LP] Make it so the name of the aminoacid shows up in its correspondig point on the graph 
#           if its frequency is higher than some value (e.g. 0.5)
#       DISCARDED - [DATA - MP] Maybe make sub-plots (every 25-50 positions, make a new graph) to compare more easily
#       [HP] Come up with more ways to compare the two plots.
#       [TOYING AROUND - LP] Try poking around and tweaking the sequences so that there exists an offset between the data 
#           (don't know what that'd be good for, but I'll just throw it in the list i guess)
#       [MP] Maybe rethink the functions so that we can work with the "frequency" and "appearance" arrays
#           (e.g. to generate various plots with different colormaps efficiently, without calculating them for every plot)