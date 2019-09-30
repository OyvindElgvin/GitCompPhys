
from matplotlib.pyplot import *
import numpy as np

Rmin = 0.0
Rmax = 8
N = 200

# sort eigenvectors and eigenvalues

r = np.linspace(Rmin,Rmax,N)
EigValues = np.loadtxt("2bA_eigvalues", skiprows=2)
EigVectors = np.loadtxt("2bR_eigenvectors", skiprows=2)
permute = EigValues.argsort()
EigValues = EigValues[permute]
EigVectors = EigVectors[:,permute]

# plots the results for the three lowest lying eigenstates
print("three lowest eigenvalues for 2b:")
for i in range(3):
    print (EigValues[i])

# sorting the eigenvalues and use that order to sort the eigenvectors
FirstEigvector = EigVectors[:,0]
SecondEigvector = EigVectors[:,1]
ThirdEigvector = EigVectors[:,2]

# plotting the three lowest eigenvectors
plot(r, FirstEigvector**2, label="Ground state")
plot(r, SecondEigvector**2, label="Second state")
plot(r, ThirdEigvector**2, label="Third state")
xlabel(r'$r$')
ylabel(r'Radial probability $r^2|R(r)|^2$')
title(r'Radial probability distributions for three lowest-lying states')
legend()
savefig('2b_eigenvector.pdf')
show()

# plotting iterations vs N
iterations = [7, 17664, 70859, 160258]
N = [4, 100, 200, 300]
plot(N,iterations)
xlabel("N")
ylabel("Number of iterations")
title("Iterations as function of N")
savefig('2b_Iterations_vs_N.pdf')
show()

# plotting time vs N
seconds = [0, 0.87475, 12.6909, 62.4196]
N = [4, 100, 200, 300]
plot(N,seconds)
xlabel("N")
ylabel("Seconds")
title("Time as function of N")
savefig('2b_time_vs_N.pdf')
show()
