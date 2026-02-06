from constants import AMINOACID_NUMBER
import numpy as np

colors = np.zeros((AMINOACID_NUMBER, 3))
for i in range(AMINOACID_NUMBER):
    colors[i] = np.random.rand(3)

np.savetxt("colormap/colormap.txt", colors)
