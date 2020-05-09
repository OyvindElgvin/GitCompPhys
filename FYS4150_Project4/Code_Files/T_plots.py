from Py_Functions import readmatrices
from numpy import array, zeros, ones, amax, where, polyfit
import matplotlib.pyplot as plt

N = 1000000000
L_values = array([40,60,80,100])
T = readmatrices("10_9medstortempspredning/Results_4e_2_N_%s_L_40_R_0.txt" % N)[0][0][:]

E_values = []; M_values = []; Cv_values = []; X_values = []

#Extracting values to arrays in values lists, i.e. [array1,array2,array3,array4]
for i in range(4):
	L = L_values[i]
	File_0 = readmatrices("10_9medstortempspredning/Results_4e_2_N_%s_L_%s_R_0.txt" % (N,L))[0]
	File_1 = readmatrices("10_9medstortempspredning/Results_4e_2_N_%s_L_%s_R_1.txt" % (N,L))[0]
	File_2 = readmatrices("10_9medstortempspredning/Results_4e_2_N_%s_L_%s_R_2.txt" % (N,L))[0]
	File_3 = readmatrices("10_9medstortempspredning/Results_4e_2_N_%s_L_%s_R_3.txt" % (N,L))[0]

	Emat  = File_0[1][:] + File_1[1][:] + File_2[1][:] + File_3[1][:]
	Mmat  = File_0[3][:] + File_1[3][:] + File_2[3][:] + File_3[3][:]
	Cvmat = File_0[4][:] + File_1[4][:] + File_2[4][:] + File_3[4][:]
	Xmat  = File_0[6][:] + File_1[6][:] + File_2[6][:] + File_3[6][:]

	E = []; M = []; Cv = []; X = []
	for j in range(len(T)):
		E.append(Emat[j][-1])
		M.append(Mmat[j][-1])
		Cv.append(Cvmat[j][-1])
		X.append(Xmat[j][-1])
	
	E_values.append(array(E))
	M_values.append(array(M))
	Cv_values.append(array(Cv))
	X_values.append(array(X))

#Calculating critical temperature
Tc_values = zeros(4)
ind = 0
for Cv in Cv_values:
	max_value = amax(Cv)
	max_index = int(where(Cv == max_value)[0][0])
	Tc_values[ind] = T[max_index]
	ind += 1

print Tc_values

L = array([40,60,80,100])
coeffs = polyfit(L,Tc_values,1)
T_inf = coeffs[-1]

print("Critical temperature using regression method: "+str(T_inf))
print "Relative deviation = ", 1-T_inf/2.269


R = []
i = 0
jmax = 4
for i in range(3):
	for j in range(i+1,jmax):
		R.append((L[i]*Tc_values[i]-L[j]*Tc_values[j])/(L[i]-L[j]))
	
R = array(R)
T_inf = sum(R)/6.

print("Critical temperature using pairs of L and averaging Cv: "+str(T_inf))
print "Relative deviation = ", 1-T_inf/2.269


#plotting and saving

values = [E_values,M_values,Cv_values,X_values]
titles = ["Expectation value of energy \n as function of temperature","Expectation value of magetic moment \n as function of temperature","Specific heat capacity \n as function of temperature", "Susceptibility as function \n of temperature"]
axes = ["$\langle E \\rangle$","$\langle |M| \\rangle$","$C_v$","$\chi_{abs}$"]
filenames = ["E_of_T","M_of_T","Cv_of_T","X_of_T"]
for i in range(4):
	plt.figure() ; plt.title(titles[i])	
	for j in range(4):
		plt.plot(T,values[i][j],label="L = %d" % L_values[j])
		plt.xlabel("T")
		plt.ylabel(axes[i])
	plt.legend()
	#plt.savefig("%s_mpi_2_2.6.pdf" % filenames[i])
plt.show()
