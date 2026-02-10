#Definition of colormap selection and generation
import numpy as np

def ColormapSelection(order, type):
    if type == "rasmol":
        return(ColormapRasmol(order))
    
    #If no type is specified, we'll assign it the (current) random one.
    else:
        return(ColormapRandom())

def ColormapRandom():
    colors = np.loadtxt("colormap/colormap.txt")
    return(colors)

def ColormapRasmol (order):
    #We'll use the color dictionary provided by Rasmol
    rasmol_rgb = {
        'A': [200, 200, 200],  # Dark Grey
        'C': [230, 230, 0],    # Yellow
        'D': [230, 10, 10],    # Bright Red 
        'E': [230, 10, 10],    # Bright Red
        'F': [50, 50, 170],    # Mid Blue
        'G': [128, 128, 185],  # Medium Grey
        'H': [130, 130, 210],  # Pale Blue
        'I': [15, 130, 15],    # Green
        'K': [20, 90, 255],    # Blue
        'L': [15, 130, 15],    # Green
        'M': [230, 230, 0],    # Yellow
        'N': [0, 220, 220],    # Cyan
        'P': [220, 150, 130],  # Flesh
        'Q': [0, 220, 220],    # Cyan
        'R': [20, 90, 255],    # Blue
        'S': [250, 150, 0],    # Orange
        'T': [250, 150, 0],    # Orange
        'V': [15, 130, 15],    # Green
        'W': [180, 90, 180],   # Purple
        'Y': [50, 50, 170],    # Mid Blue
        '-': [190, 160, 110]   # Tan (Using "Other" for the dash)
    }

    colors = []  #We'll store each aa in this list

    for aa in (order):
        rgb = rasmol_rgb.get (aa)  #We pick the color of the aa from the dictionary

        if rgb:
            rgb_normalizado = [valor / 255.0 for valor in rgb]  #We must normalize it because Matplotlib
            colors.append(rgb_normalizado)                     
        else:
            colors.append([1.0, 1.0, 1.0])  #If there's any error, we assign it black.

    return(np.array(colors))