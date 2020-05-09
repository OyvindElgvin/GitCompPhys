from numpy import array, zeros, linspace, log, log10, exp, polyfit
from matplotlib.pyplot import *
from timeit import default_timer as timer
from scipy.special import gamma, factorial

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


def plot_median_c(filename):
	median = readarrays(filename)[0][0]
	D = filename.split("_")
	D = int(D[5])
	filename = filename.split("/")
	filename = filename[2].split(".txt")
	filename = filename[0]
	N = log10( (10**(D-2)) * array(range(0,len(median))) )
	N[0] = 0.0

	figure()
	title("Median for one experiment")
	plot(N,median,'.',N,median,'r')
	xlabel("$\log_{10}(MC)$")
	ylabel("$\mu_{1/2}$")
	savefig("../Plots/Plots_a_b_c/" + filename + ".png", dpi=300)
	show()
	close()


def P(L,x):
	n = 1 + 3*float(L) / (1-float(L))
	a = n**n / gamma(n)
	p = a * x**(n-1) * exp(-n * x)
	return p



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
