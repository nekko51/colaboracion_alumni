#Definición de la función que a cada aminoacid le asigna su color siguiente la clasificación RasMol
import numpy as np

def coloring_function (orden):

    #Usamos el diccionario de colores otorgado por RasMol

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
        '-': [190, 160, 110]   # Tan (Usando la categoría 'Others' para el hueco)
    }

    colores = []  #Creamos una lista vacía donde vamos a incluir el color de cada aminoácido

    for aa in (orden):
        rgb = rasmol_rgb.get (aa)  #Cogemos un aminoacid del diccionario

        if rgb:
            rgb_normalizado = [valor / 255.0 for valor in rgb]  #Debemos normalizarlo pues Matplotlib
            colores.append(rgb_normalizado)                     #espera valores entre 0 y 1.
        else:
            colores.append([1.0, 1.0, 1.0])  #Si hay algún defecto se le asigna el color negro

    return np.array(colores)