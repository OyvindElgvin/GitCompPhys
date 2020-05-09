from Py_Functions_DY_lambda import readarrays, plot_median_c, nlog_arrays, log_arrays, P
import Py_Functions_DY_lambda as pf
from numpy import array, zeros, linspace, log, log10, exp, sort, polyfit, polyval
from matplotlib.pyplot import *
import matplotlib.image as pmg
import scipy.stats as stats


abc = input("Which task would you like to run? a/b/c \n")

"""
task a)
"""
if (abc == "a"):
	# Gibbs distribution
	N = range(1,501)
	Money = []
	M_R = readarrays("../Results_O/OE_a_Money_distributions_D_5_N_500_0.txt")[0]
	for i in range(0,len(M_R)):
		Money += list(M_R[i])
	unused, m_intervals, patches = hist(Money,len(N)*2,density=True, label="Numerical")
	plot(m_intervals, exp(-m_intervals), label = "Analytical", linewidth = 0.9)
	title("Wealth distribution with $N$ = 500")
	xlabel("Wealth m ")
	ylabel("$w_m$")
	legend()
	savefig("../Plots/Plots_a_b_c/a_histogram.png", dpi=300)
	show()

	# Medianer
	plot_median_c("../Results_O/OE_a_Median_D_5_N_500.txt")


"""
task b)
"""
if abc == "b":
	# Linear regression of log plot
	N = range(1,501)
	Money = []
	M_R = readarrays("../Results_O/OE_a_Money_distributions_D_5_N_500_0.txt")[0]
	for i in range(0,len(M_R)):
		Money += list(M_R[i])
	fig = figure()
	w, m_intervals, patches = hist(Money,len(N),density=True)
	close(fig)
	m = m_intervals[:-1]
	p, w = nlog_arrays(m, w)
	x, y = polyfit(exp(p[:-60]),w[:-60],1)
	plot(exp(p),w)
	plot(exp(p), x*exp(p)+y, label = "Gibbs ditribution = {0:.4f} $\cdot$ m + {1:.3f} ".format(x, y))
	title("Equilibrium, $w_m$, as function of $m$")
	xlabel("Wealth m")
	ylabel("$log(w_m)$")
	legend()
	savefig("../Plots/Plots_a_b_c/b_logplot_nlog.png", dpi=300)
	show()






"""
task c
"""
if abc == "c":

	# Medianer
	D = 5
	filenames = []
	filestart = "../Results_O/OE_c_Median_D_%s" % D
	Nvalues = ["500"]
	Lvalues = ["0.250000", "0.500000", "0.900000"]
	for L in Lvalues:
		filenames.append(filestart + "_N_" + Nvalues[0] + "_L_" + L + ".txt")

	for f in filenames:
		plot_median_c(f)


	# Three seperate plots of lambda superimposed on histograms
	N = range(1,501)
	Money0 = []
	filenames = []
	filestart = "../Results_O/OE_c_Money_distributions_D_5_N_500"
	Lvalues = ["0.250000", "0.500000", "0.900000"]
	for L in Lvalues:
		filenames.append(filestart + "_L_" + L + "_0" + ".txt")
	for R in range(len(Lvalues)):
		Money0 = []
		M_R = readarrays("%s" % filenames[R])[0]
		for i in range(0,len(M_R)):
			Money0 += list(M_R[i])
		figure()
		unused, m_intervals0, patches = hist(Money0,len(N),density=True, label="Numerical")
		plot(m_intervals0, P(Lvalues[R],m_intervals0), label = "Analytical")
		title("Distribution of wealth for $\lambda = ${0:.2f}".format(float(Lvalues[R])))
		xlabel("Wealth m $(m_0 = 1)$")
		ylabel("$P(m)$")
		legend()
		savefig("../Plots/Plots_a_b_c/c_money_dist_hist" + "_L_" + Lvalues[R] + ".png", dpi=300)
		show()
		close()


	# All lambdas in the same plot
	# Numerical Lambda = 0
	N = range(1,501)
	Money1 = []
	M_R = readarrays("../Results_O/OE_a_Money_distributions_D_5_N_500_0.txt")[0]
	for i in range(0,len(M_R)):
		Money1 += list(M_R[i])
	fig = figure()
	w1, m_intervals1, patches = hist(Money1,len(N),density=True, label="Numerical $\lambda = 0$")
	close(fig)
	plot(m_intervals1[:-1], w1, '.', markersize=2, label="Numerical $\lambda = 0$")
	# Numerical Lambda = 0.25, 0.50, 0.90
	filenames = []
	filestart = "../Results_O/OE_c_Money_distributions_D_5_N_500"
	Lvalues = ["0.250000", "0.500000", "0.900000"]
	for L in Lvalues:
		filenames.append(filestart + "_L_" + L + "_0" + ".txt")
	for R in range(len(Lvalues)):
		Money2 = []
		M_R = readarrays("%s" % filenames[R])[0]
		for i in range(0,len(M_R)):
			Money2 += list(M_R[i])
		fig = figure()
		w2, m_intervals2, patches = hist(Money2,len(N),density=True, label="Numerical")
		close(fig)
		plot(m_intervals2[:-1], w2, '.', markersize=2, label= "Numerical $\lambda = ${0:.2f}".format(float(Lvalues[R])))
	# Analytical lambda = 0
	N = range(1,501)
	Money3 = []
	M_R = readarrays("../Results_O/OE_a_Money_distributions_D_5_N_500_0.txt")[0]
	for i in range(0,len(M_R)):
		Money3 += list(M_R[i])
	fig = figure()
	un_used, m_intervals3, patches = hist(Money3,len(N),density=True, label="Numerical")
	close(fig)
	plot(m_intervals3, exp(-m_intervals3), label = "Analytical $\lambda = 0$")
	# Analytical Lambda = 0.25, 0.5, 0.9
	for i in range(len(Lvalues)):
		plot(m_intervals1, P(Lvalues[i],m_intervals1), label = "Analytical $\lambda = ${0:.2f}".format(float(Lvalues[i])))
	title("Distribution of wealth for different $\lambda$")
	xlabel("Wealth m")
	ylabel("$P_n(m)$")
	legend()
	xlim(0,2)
	savefig("../Plots/Plots_a_b_c/c_money_dist_L_123.png", dpi=300)
	show()


	# plotting tail end for each lambda
	N = range(1,501)
	filenames = []
	filestart = "../Results_O/OE_c_Money_distributions_D_5_N_500"
	Lvalues = ["0.250000", "0.500000", "0.900000"]
	for L in Lvalues:
		filenames.append(filestart + "_L_" + L + "_0" + ".txt")
	for R in range(len(Lvalues)):
		Money2 = []
		M_R = readarrays("%s" % filenames[R])[0]
		for i in range(0,len(M_R)):
			Money2 += list(M_R[i])
		fig = figure()
		ww, m_intervals2, patches = hist(Money2,len(N),density=True)
		close(fig)
		m2 = m_intervals2[:-1]
		p2, www = log_arrays(m2, ww)


		# loglog plot
		start1 = (100,140,230)
		stop1 = (-170,-170,-110)
		# tail for lambda = 0.25 with fitting
		if R == 0:
			slope, intercept, r_value, p_value, std_err = stats.linregress(p2[start1[0]:stop1[0]], www[start1[0]:stop1[0]])
			plot(p2,www, label = "$w_m$")
			plot(p2, slope*p2+intercept, label = "{0:.2f} $\cdot$ log10(m) + {1:.2f}".format(float(slope),float(intercept)))
			title("Linear regression of $w_m$ vs $m$ with $\lambda$ = {0:.2f}".format(float(Lvalues[0])))
			xlabel("log$_{10}$(m)")
			ylabel("log$_{10}(w_m)$")
			legend()
			xlim(-2.5,1)
			ylim(-4.5,0)
			savefig("../Plots/Plots_a_b_c/c_money_dist_L_" + Lvalues[0] + "_tail_loglog.png", dpi=300)
			show()
		# tail for lambda = 0.50 with fitting
		if R == 1:
			slope1, intercept1, r_value, p_value, std_err = stats.linregress(p2[start1[1]:stop1[1]], www[start1[1]:stop1[1]])
			plot(p2,www, label = "$w_m$")
			plot(p2, slope1*p2+intercept1, label = "{0:.2f} $\cdot$ log10(m) + {1:.2f}".format(float(slope1),float(intercept1)))
			title("Linear regression of $w_m$ vs $m$ with $\lambda$ = {0:.2f}".format(float(Lvalues[1])))
			xlabel("log$_{10}$(m)")
			ylabel("log$_{10}(w_m)$")
			legend()
			xlim(-1.5,0.75)
			ylim(-4.5,0)
			savefig("../Plots/Plots_a_b_c/c_money_dist_L_" + Lvalues[1] + "_tail_loglog.png", dpi=300)
			show()
		# tail for lambda = 0.90 with fitting
		if R == 2:
			slope2, intercept2, r_value, p_value, std_err = stats.linregress(p2[start1[2]:stop1[2]], www[start1[2]:stop1[2]])
			plot(p2,www, label = "$w_m$")
			plot(p2, slope2*p2+intercept2, label = "{0:.2f} $\cdot$ log10(m) + {1:.2f}".format(float(slope2),float(intercept2)))
			title("Linear regression of $w_m$ vs $m$ with $\lambda$ = {0:.2f}".format(float(Lvalues[2])))
			xlabel("log$_{10}$(m)")
			ylabel("log$_{10}(w_m)$")
			legend()
			xlim(-0.5,0.35)
			ylim(-4,1)
			savefig("../Plots/Plots_a_b_c/c_money_dist_L_" + Lvalues[2] + "_tail_loglog.png", dpi=300)
			show()
