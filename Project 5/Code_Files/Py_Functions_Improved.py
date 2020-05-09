from numpy import array, zeros, linspace, log, log10, exp, polyfit
import matplotlib.pyplot as plt
from timeit import default_timer as timer

def readarrays(filename):
	#start = timer()
	values = open(filename, "r")
	#print values.read()
	lines = values.readlines()
	#end = timer()
	#print (end-start)


	#Counting
	C = 0
	D = 0
	Dims = []
	A = []

	#start = timer()
	for i in lines:
		if i != "\n":
			D += 1
		if i == "\n":
			C += 1
			Dims.append(D)
			A.append(zeros(D))
			D = 0
	#end = timer()
	#print (end-start)

	#start = timer()
	#Filling
	F = 0
	G = 0
	for i in lines:
		if i != "\n":
			A[F][G] = i
			G += 1
		if i == "\n":
			F += 1
			G = 0
	#end = timer()
	#print (end-start)
	values.close()
	return A,len(A)

def readmatrices(filename):
	#start = timer()
	values = open(filename, "r")
	#print values.read()
	lines = values.readlines()
	#end = timer()
	#print (end-start)


	#Counting
	C = 0
	D = 0
	Dims = []
	A = []

	#start = timer()

	for i in range(len(lines)):
		if lines[i] != "\n":
			D += 1
		if lines[i] == "\n":
			C += 1
			Dims.append(D)
			A.append(zeros(shape=(D,len(lines[i-1].split()))))
			D = 0


	#start = timer()
	#Filling
	F = 0
	G = 0
	for i in lines:
		if i != "\n":
			for j in range(len(i.split())):
				A[F][G][j] = i.split()[j]
			G += 1
		if i == "\n":
			F += 1
			G = 0
	#end = timer()
	#print (end-start)
	values.close()

	return A,len(A)

def log_arrays(M,W):
	logM = [] ; logW = []
	for i in range(len(M)):
		if W[i] != 0.:
			logM.append(log10(M[i]))
			logW.append(log10(W[i]))
	return array(logM), array(logW)

def nlog_arrays(M,W):
	logM = [] ; logW = []
	for i in range(len(M)):
		if W[i] != 0.:
			logM.append(log(M[i]))
			logW.append(log(W[i]))
	return array(logM), array(logW)

def extract_parametres(string):
	f = string.split("_")
	params = []

	for i in range(1,len(f),2):
		if f[i-1] == "D" or f[i-1] == "N":
			params.append(int(f[i]))
		elif i != len(f)-1:
			params.append(float(f[i]))
		else:
			params.append(float(f[i][:-4]))
	#D = int(f[3]) ; N = int(f[5]) ; L = float(f[7]) ; a = float(f[9][:-4])
	return params

def plot_median(filename,savefile,D,N=1000,L=0,a=0,g=0,save=False):
	median = readarrays(filename)[0][0]

	MC = log10( (10**(D-2)) * array(range(0,len(median))) )
	MC[0] = 0.
	figtitle = "Median Plots \n"
	savetitle = savefile
	var = [D,N,L,a,g]
	varnames = ["D","N","L","a","g"]
	for i in range(len(var)):
		if var[i] != 0:
			figtitle += varnames[i] + "%.1f" % var[i] + ", "
			savefile += "_" + varnames[i] + "_" + "%.1f" % var[i] + " "

	plt.figure()
	plt.title(figtitle)
	plt.plot(N,median,'.',N,median,'r')
	plt.xlabel("$\log_{10}(MC)$")
	plt.ylabel("$\mu_{1/2}$")


	'''
	fig, (ax1, ax2) = plt.subplots(1,2,sharey=True,gridspec_kw={'width_ratios': [1,5]})
	fig.suptitle(filename)
	ax1.plot(N,median,'b.',N,median,'r')
	ax2.plot(N,median,'b.',N,median,'r')


	ax1.set_xlim(-0.01,0.01)
	ax2.set_xlim(2.9,5)
	ax1.set_xticks([0])

	ax1.spines['right'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	ax1.yaxis.tick_left() ; ax1.tick_params(labelright='off')
	ax2.yaxis.tick_right() ; ax2.set_yticks([])

	d = 0.015
	kwargs = dict(transform=ax1.transAxes, color='k',clip_on=False)
	ax1.plot( (1-d,1+d), (-d,+d), **kwargs )
	ax1.plot( (1-d,1+d), (1-d,1+d), **kwargs )

	kwargs.update(transform=ax2.transAxes)
	ax2.plot( (-d,+d), (1-d,1+d), **kwargs )
	ax2.plot( (-d,+d), (-d,+d), **kwargs )

	ax2.set_xlabel("$\log_{10}(MC)$")
	ax1.set_ylabel("$\mu_{1/2}$")
	'''
	if save:
		plt.savefig(savetitle + ".png")
	else:
		plt.show()

	plt.close()

def plot_prob_distribution(filename,savefile,D,N=1000,L=0,a=0,g=0,save=False):
	Money = []

	figtitle = "Median Plots \n"
	savetitle = savefile
	var = [D,N,L,a,g]
	varnames = ["D","N","L","a","g"]
	for i in range(len(var)):
		if var[i] != 0:
			figtitle += varnames[i] + " = %.1f" % var[i] + ", "
			savetitle += "_" + varnames[i] + "_" + "%.1f" % var[i]
	print filename
	print savetitle

	for R in range(4):
		M_R = readarrays("%s_%s.txt" % (filename,R))[0]
		for i in range(0,len(M_R)):
			Money += list(M_R[i])

	Money = array(Money)

	plt.figure()
	w, m_intervals, patches = plt.hist(Money,N,density=True)


	plt.title(figtitle)
	plt.xlabel("Amount of money")
	plt.ylabel("$w(m)$")

	if save:
		plt.savefig(savetitle +  ".png")
	else:
		plt.show()
	plt.close()
	return w, m_intervals

def VAT_hist():
	W = [] ; M = []
	Tvalues = ["0.100000","0.250000","0.500000","0.900000"]
	for T in Tvalues:
		f = "../Results/Taxes_VAT_D_5_N_500_t_" + T

		w, m = plot_prob_distribution(f,"../Plots/Plots_tax/VAT_Histogram",5,N=500,L=float(T),save=True)

		W.append(w) ; M.append(m)
	plt.figure()
	for i in range(len(M)):
		plt.plot(M[i][:-1],W[i],label="t = %.2f" % float(Tvalues[i]))
	plt.legend()
	plt.xlabel("Amount of money")
	plt.ylabel("$w(m)$")
	plt.axis([0,6,0,1])
	plt.title("VAT histogram as plots")
	plt.savefig("../Plots/Plots_tax/VAT_hist_plot.png")

def wealthtax_hist():
	W = [] ; M = []
	Tvalues = ["0.100000","0.250000","0.500000","0.900000"]
	for T in Tvalues:
		f = "../Results/Taxes_Wealth_D_5_N_500_t_" + T

		w, m = plot_prob_distribution(f,"../Plots/Plots_tax/Wealthtax_Histogram",5,N=500,L=float(T),save=True)

		W.append(w) ; M.append(m)
	plt.figure()
	for i in range(len(M)):
		plt.plot(M[i][:-1],array(W[i])/2500.,label="t = %.2f" % float(Tvalues[i]))
	plt.legend()
	plt.title("Wealth Tax Histogram as Plots")
	plt.xlabel("Amount of money")
	plt.ylabel("$w(m)$")
	plt.axis([0.99,1.01,0,1])
	plt.savefig("../Plots/Plots_tax/WealthTax_hist_plot.png")
