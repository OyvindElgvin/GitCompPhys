

from matplotlib.pyplot import *
import numpy as np

cols = []

f = open("2dA_eigvalues","r")
f.readline()
f.readline()
for line in f:
    words = line.split()
    cols.append(float(words[0]))


#cols.sort()
#print(cols)



R = np.loadtxt("2dR", skiprows=2)
#print(R)

GS_index = np.argmin(np.array(cols))

GS_vec = R[:, GS_index]
print("jacobi 2d lowest eigenvalu:")
print(cols[GS_index])
#print(GS_index, GS_vec)
