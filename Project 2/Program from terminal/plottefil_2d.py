

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
#print("jacobi 2d lowest eigenvalu:")
#print(cols[GS_index])
#print(GS_index, GS_vec)

# makes a plot of eigenvalues vs rho_max

Rmax =           [       4        ,       6         ,    8           ,      10         ,       12        ,         14      ,        16       ,      18         ,        20]
eigenvalue_200 = [2.99990030309549, 2.99971872151564, 2.9994999099550, 2.99921853010167, 2.99887454386168, 2.99846790460321, 2.99799855711194, 2.99746643753911, 2.99687147334007]
eigenvalue_150 = [2.99980189487923, 2.99949990995511, 2.9991108264135, 2.99861041570032, 2.99799855711192, 2.99727510270977, 2.99643987707005, 2.99549267698529, 2.99443327111675]
eigenvalue_100 = [2.99952186709985, 2.99887454386175, 2.9979985571119, 2.99687147334005, 2.99549267698525, 2.99386141093105, 2.99197677357179, 2.98983771527599, 2.98744303420867]
eigenvalue_80  = [2.99923899185503, 2.99824107317421, 2.9968714733400, 2.99510856223066, 2.99295082575143, 2.99039639489375, 2.98744303420868, 2.98408812774886, 2.98032866223699]
eigenvalue_60  = [2.99862862261759, 2.99687147334012, 2.9944332711167, 2.99129207930722, 2.98744303420867, 2.98288008419831, 2.97759592063636, 2.97158189142357, 2.96482789432465]


plot(Rmax, eigenvalue_200, label="N = 200")
plot(Rmax, eigenvalue_150, label="N = 150")
plot(Rmax, eigenvalue_100, label="N = 100")
plot(Rmax, eigenvalue_80, label="N = 80")
plot(Rmax, eigenvalue_60, label="N = 60")
title("eig func of rho_max")
xlabel("Rho_max")
ylabel("Eigenvalues")
legend()
show()
