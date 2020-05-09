#from Py_Functions import readarrays, plot_median_d, plot_prob_distribution_d, log_arrays, nlog_arrays
import Py_Functions as pf
import Py_Functions_Improved as pfi
from numpy import array, zeros, linspace, log, log10, exp, sort, polyfit, polyval
import matplotlib.pyplot as plt
import matplotlib.image as pmg

import scipy.stats as stats
'''
N = range(1,501)
Money = []

#filename = "Money_distributions_savings_L_0.250000"
filename = "Test_sampling"

for R in range(4):
	M_R = readarrays("../build-Laptop_Project-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/%s_%s.txt" % (filename,R))[0]
	for i in range(0,len(M_R)):
		Money += list(M_R[i])


w, m_intervals, patches = plt.hist(Money,len(N),density=True)

'''
'''
coeff = polyfit(m_intervals[:-1],w,10)
W = polyval(coeff,m_intervals[:-1])
plt.plot(m_intervals[:-1],W,'r')

fit_alpha, fit_loc, fit_beta=stats.gamma.fit(w)

plt.plot()
'''

'''
#plt.plot(m_intervals[:-1],w,'r')

#plt.savefig("histogram_a.png")

lw = zeros(len(w))
for i in range(len(lw)):
	if w[i] != 0:
		lw[i] = log(w[i])
	else:
		lw[i] = w[i]

plt.figure()
plt.plot(log10(m_intervals[:-1]),log10(w))
#plt.savefig("log_w_a.png")

slope = (lw[100]-lw[1])/(m_intervals[100]-lw[1])
print slope

plt.show()


'''

#pfi.wealthtax_hist()
pfi.VAT_hist()

D = 3


filenames = []
filestart = "../Results/Median_D_%s" % D

Nvalues = ["500", "1000"]
Lvalues = ["0.250000", "0.500000", "0.900000"]
avalues = ["0.500000", "1.000000", "1.500000", "2.000000"]
gvalues = ["1.000000", "2.000000", "3.000000", "4.000000"]

#Meidans taxes
yn = raw_input("Do you want to plot medians for taxes? y/n \n")
if yn == "y":
	filenames = []
	D = raw_input("Exponent of D? \n")
	filestart = "../Results/Median_Taxes_D_%s_N_500" % D
	for T in ["0.100000","0.250000","0.500000","0.900000"]:
		filenames.append(filestart + "_t_" + T + ".txt")

	for f in filenames:
		plot_median_d(f,save=False)



#Medians task d
yn = raw_input("Do you want to plot medians for task d)? y/n \n")
if yn == "y":
	filenames = []
	D = raw_input("Exponent of D? \n")
	filestart = "../Results/Median_D_%s" % D
	for N in Nvalues:
		for L in ["0.000000","0.500000"]:
			for a in avalues:
				filenames.append(filestart + "_N_" + N + "_L_" + L + "_a_" + a + ".txt")

	for f in filenames:
		plot_median_d(f,save=False)

#Probability distribution task d
yn = raw_input("Do you want to plot probability distributions for d)? y/n \n")
if yn == "y":
	filenames = []
	D = raw_input("Exponent of D? \n")
	filestart = "../Results/Test_Money_distributions_D_%s" % D
	for N in Nvalues:
		for L in ["0.000000","0.500000"]:
			for a in avalues:
				filenames.append(filestart + "_N_" + N + "_L_" + L + "_a_" + a + ".txt")


	W = [] ; M = []
	for i in range(len(filenames)):
		f = filenames[i]
		w, m = plot_prob_distribution_d(f,save=False)
		W.append(w) ; M.append(m[:-1])

	i = 0
	while i < len(W):

		pf.plot_loglogW_d(M[i:i+4],W[i:i+4],filenames[i:i+4],save=False)
		pf.Pareto_dist_d(M[i:i+4],W[i:i+4],filenames[i:i+4],save=False)
		pf.Gibbs_dist_d(M[i:i+4],W[i:i+4],filenames[i:i+4],save=False)
		i += 4

#Median task e
yn = raw_input("Do you want to plot medians for task e)? y/n \n")
if yn == "y":
	filenames = []
	D = raw_input("Exponent of D? \n")
	filestart = "../Results/Median_D_%s_N_1000" % D
	for L in ["0.000000","0.500000"]:
		for a in ["1.000000","2.000000"]:
			for g in gvalues:
				filenames.append(filestart + "_L_" + L + "_a_" + a + "_g_" + g + ".txt")
	filenames.append("../Results/Median_D_8_N_1000_L_0.500000_a_2.000000_g_1.000000.txt")
	#for f in filenames:
	#	pf.plot_median_e(f,save=True)
	pf.plot_median_e(filenames[-1],save=True)



yn = raw_input("Do you want to plot probability distributions for e)? y/n \n")
if yn == "y":
	filenames = []
	D = raw_input("Exponent of D? \n")
	filestart = "../Results/TEST_Money_distributions_D_%s_N_1000" % D
	for L in ["0.000000","0.500000"]:
		for a in ["1.000000","2.000000"]:
			for g in gvalues:
				filenames.append(filestart + "_L_" + L + "_a_" + a + "_g_" + g + ".txt")
				#if N == "1000" and L == "0.000000" and a == "2.000000":
				#	filenames.append(filestart + "_N_" + N + "_L_" + L + "_a_" + "9.000000" + ".txt")


	W = [] ; M = []
	for i in range(len(filenames)):
		f = filenames[i]
		w, m = pf.plot_prob_distribution_e(f,save=False)
		W.append(w) ; M.append(m[:-1])

	i = 0
	while i < len(W):

		pf.plot_loglogW_e(M[i:i+4],W[i:i+4],filenames[i:i+4],save=False)
		pf.Pareto_dist_e(M[i:i+4],W[i:i+4],filenames[i:i+4],save=True)
		#pf.Gibbs_dist_e(M[i:i+4],W[i:i+4],filenames[i:i+4],save=False)
		i += 4
