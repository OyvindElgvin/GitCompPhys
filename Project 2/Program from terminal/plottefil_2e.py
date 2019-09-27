
from matplotlib.pyplot import *
import numpy as np

cols = []

f = open("2eA_eigvalues","r")
f.readline()
f.readline()
for line in f:
    words = line.split()
    cols.append(float(words[0]))


#cols.sort()
#print(cols)



R = np.loadtxt("2eR", skiprows=2)
#print(R)

GS_index = np.argmin(np.array(cols))

GS_vec = R[:, GS_index]
Rmax = np.linspace(4,50,100)

print("jacobi 2e lowest eigenvalu:")
#print(cols[GS_index])
print(GS_index, GS_vec)
print(Rmax)


plot(Rmax, GS_vec) per omega
show()

"""
legg inn forskjelling omega over her
"""

"""
cols.sort()
print("test of eigenvalues")
print(cols[0])
print(cols[1])
print(cols[2])
print(cols[3])
"""
