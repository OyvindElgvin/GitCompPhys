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

def plot_median_d(filename,save=False):
	median = readarrays(filename)[0][0]
	D = filename.split("_") ; D = int(D[3])
	#filename = filename.split("/") ; filename = filename[2].split(".txt") ; filename = filename[0]
	N = log10( (10**(D-2)) * array(range(0,len(median))) )
	N[0] = 0.
	#print filename
	
	plt.figure()
	#plt.title(filename)
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
		Folder = "Plots_d_D_%s/" % D
		#plt.savefig("../Plots/Plots_d/" + Folder + filename + ".png")	
	else:
		plt.show()

	plt.close()

def plot_prob_distribution_d(filename,save=False):
	
	if filename[filename.index("N")+ 2] == "5":
		N = filename[filename.index("N")+ 2 : filename.index("N")+ 5]
	else: 
		N = filename[filename.index("N")+ 2 : filename.index("N")+ 6]
	N = int(N)
	
	Money = []

	for R in range(4):
		M_R = readarrays("%s_%s.txt" % (filename[:-4],R))[0]
		for i in range(0,len(M_R)):
			Money += list(M_R[i])
	
	Money = array(Money)

	plt.figure()
	w, m_intervals, patches = plt.hist(Money,N,density=True)


	filename = filename.split("/") ; filename = filename[2].split(".txt") ; filename = filename[0]
	plt.title(filename)
	plt.xlabel("Amount of money")
	plt.ylabel("$w(m)$")
	
	if save:
		Folder = "Plots_d_D_%s/" % D
		plt.savefig("../Plots/Plots_d/" + Folder + filename + ".png")
	#else:
	#	plt.show()	
	plt.close()
	return w, m_intervals

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

def extract_parametres(filename):
	f = filename.split("/")[2]
	f = f.split("_")
	params = []

	for i in range(4,len(f),2):
		if f[i-1] == "D" or f[i-1] == "N":		
			params.append(int(f[i]))
		elif i != len(f)-1:
			params.append(float(f[i]))
		else:
			params.append(float(f[i][:-4]))
	#D = int(f[3]) ; N = int(f[5]) ; L = float(f[7]) ; a = float(f[9][:-4])
	return params


def plot_loglogW_d(M,W,filenames,save=False):

	D, N, L, a = extract_parametres(filenames[0])
	plt.figure()
	plt.title("log-plot of w(m), N = %d, L = %f" % (N,L) )
	for i in range(len(filenames)):
		D, N, L, a = extract_parametres(filenames[i])

		logM, logW = log_arrays(M[i],W[i])
		plt.plot(logM,logW,".",markersize=0.8,label="$\\alpha = %f $" % a)
	plt.legend()
	plt.axis([-1,2,-5,0.5])
	plt.xlabel("$\log_{10}(m)$") ; plt.ylabel("$\log_{10}(w)$")
	if save:
		Folder = "Plots_d_D_%s/" % D
		plt.savefig("../Plots/Plots_d/" + Folder + "Money Distribution_D_%d_N_%d_L_%.1f" % (D,N,L) + ".png")
	else:
		plt.show()
	plt.close()	



def Pareto_dist_d(M,W,filenames,save=False):
	D, N, L, a = extract_parametres(filenames[0])
	fig, ax = plt.subplots(2,2)	
	fig.suptitle("Linear Regression of log-log plot for $\omega$ \n N = %d, L = %.1f" % (N,L) )
	fig.text(0.5, 0.04, '$\log{10}(m)$', ha='center')
	fig.text(0.04, 0.5, '$\log_{10}(\omega)$', va='center', rotation='vertical')	


	SP = []
	for i in range(len(ax)):
		for j in range(len(ax[0])):
			SP.append(ax[i][j])
	
	for i in range(len(filenames)):
		D, N, L, a = extract_parametres(filenames[i])
		logM,logW = log_arrays(M[i],W[i])
		
		SI = 0 #start index for linreg
		mtemp = 10		
		eps = 0.1		
		for j in range(len(logM)):
						
			if abs(logM[j]) < abs(mtemp):
				SI = j
				mtemp = logM[j]

		print SI

		EI = len(logM)#-200 # end index
		coeff = polyfit(logM[SI:EI],logW[SI:EI],1)
		Linreg = coeff[0]*logM + coeff[1]
		
		
		SP[i].plot(logM,logW,label="$\omega (m ; \\alpha = %.1f )$" % a)
		SP[i].plot(logM,Linreg,label="$ %.2f \cdot \log_{10}(m) + %.2f $" % (coeff[0],coeff[1]) )
		SP[i].legend()
		SP[i].set_xlim(-1,2)
		SP[i].set_ylim(-5,0.5)

	if save:
		Folder = "Plots_d_D_%s/" % D
		plt.savefig("../Plots/Plots_d/" + Folder + "Linear_Regression_Pareto_D_%d_N_%d_L_%.1f" % (D,N,L) + ".png")
	else:
		plt.show()
	plt.close()	
	return

def Gibbs_dist_d(M,W,filenames,save=False):
	D, N, L, a = extract_parametres(filenames[0])
	fig, ax = plt.subplots(2,2)	
	fig.suptitle("Linear Regression of log plot for $\omega$ \n N = %d, L = %.1f" % (N,L) )
	fig.text(0.5, 0.04, '$(m)$', ha='center')
	fig.text(0.04, 0.5, '$\log(\omega)$', va='center', rotation='vertical')	


	SP = []
	for i in range(len(ax)):
		for j in range(len(ax[0])):
			SP.append(ax[i][j])
	
	for i in range(len(filenames)):
		D, N, L, a = extract_parametres(filenames[i])
		nlogM,nlogW = nlog_arrays(M[i],W[i])
		
		m = exp(nlogM)
		SI = 0 #start index for linreg
		mtemp = 10		
		eps = 0.1		
		'''		
		for j in range(len(logM)):
						
			if abs(logM[j]) < abs(mtemp):
				SI = j
				mtemp = logM[j]

		#print SI
		'''
		SI = 0
		EI = len(m)#-200 # end index
		coeff = polyfit(m[SI:EI],nlogW[SI:EI],1)
		Linreg = coeff[0]*m + coeff[1]
		
		
		SP[i].plot(m,nlogW,label="$\omega (m ; \\alpha = %.1f )$" % a)
		SP[i].plot(m,Linreg,label="$ %.2f \cdot m + %.2f $" % (coeff[0],coeff[1]) )
		SP[i].legend()
		#SP[i].set_xlim(-1,2)
		#SP[i].set_ylim(-5,0.5)

	if save:
		Folder = "Plots_d_D_%s/" % D
		plt.savefig("../Plots/Plots_d/" + Folder + "Linear_Regression_Gibbs_D_%d_N_%d_L_%.1f" % (D,N,L) + ".png")
	else:
		plt.show()
	plt.close()	
	return


def plot_median_e(filename,save=False):
	median = readarrays(filename)[0][0]
	D = filename.split("_") ; D = int(D[2])
	filename = filename.split("/") ; filename = filename[2].split(".txt") ; filename = filename[0]
	N = log10( (10**(D-2)) * array(range(0,len(median))) )
	N[0] = 0.
	print filename
	
	'''
	plt.figure()
	plt.title(filename)
	plt.plot(N,median,'.',N,median,'r')
	plt.xlabel("$\log_{10}(MC)$")
	plt.ylabel("$\mu_{1/2}$")
	'''
	
	
	fig, (ax1, ax2) = plt.subplots(1,2,sharey=True,gridspec_kw={'width_ratios': [1,5]})
	fig.suptitle(filename)
	ax1.plot(N,median,'b.',N,median,'r')
	ax2.plot(N,median,'b.',N,median,'r')
	

	ax1.set_xlim(-0.01,0.01)
	ax2.set_xlim(3.5,8)
	ax1.set_xticks([0])

	ax1.spines['right'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	ax1.yaxis.tick_left() ; ax1.tick_params(labelright='off')
	ax2.yaxis.tick_right() #; ax2.set_yticks([])

	d = 0.015
	kwargs = dict(transform=ax1.transAxes, color='k',clip_on=False)
	ax1.plot( (1-d,1+d), (-d,+d), **kwargs )
	ax1.plot( (1-d,1+d), (1-d,1+d), **kwargs )

	kwargs.update(transform=ax2.transAxes)
	ax2.plot( (-d,+d), (1-d,1+d), **kwargs )
	ax2.plot( (-d,+d), (-d,+d), **kwargs )

	ax2.set_xlabel("$\log_{10}(MC)$")
	ax1.set_ylabel("$\mu_{1/2}$")
	
	if save:
		Folder = "Plots_e_D_%s/" % D
		plt.savefig("../Plots/Plots_e/" + Folder + filename + ".png")	
	else:
		plt.show()

	plt.close()

def plot_prob_distribution_e(filename,save=False):
	D,N,L,a,g = extract_parametres(filename)
	
	Money = []
	
	for R in range(4):
		M_R = readarrays("%s_%s.txt" % (filename[:-4],R))[0]
		for i in range(0,len(M_R)):
			Money += list(M_R[i])
	
	Money = array(Money)
	nCols = int(N/0.5)
	plt.figure()
	w, m_intervals, patches = plt.hist(Money,nCols,density=True)

	plt.title("Probability Distribution \n N = %d, L=%.1f, $\\alpha = %.1f$, $\gamma = %.1f$" % (N,L,a,g))
	plt.xlabel("Amount of money")
	plt.ylabel("$w(m)$")
	
	if save:
		Folder = "Plots_e_D_%s/" % D
		plt.savefig("../Plots/Plots_e/" + Folder + "prob_dist_N_%d_L_%.1f_a_%.1f_g_%.1f" % (N,L,a,g) + ".png")
	#else:
	#	plt.show()	
	plt.close()
	return w, m_intervals




def plot_loglogW_e(M,W,filenames,save=False):

	D, N, L, a, g = extract_parametres(filenames[0])
	plt.figure()
	plt.title("log-plot of w(m), N = %d, L = %f, $\\alpha = %.1f$" % (N,L,a) )
	for i in range(len(filenames)):
		g = extract_parametres(filenames[i])[-1]

		logM, logW = log_arrays(M[i],W[i])
		plt.plot(logM,logW,".",markersize=0.8,label="$\gamma = %.1f $" % g)
	plt.legend()
	plt.axis([-3,2,-5,2])
	plt.xlabel("$\log_{10}(m)$") ; plt.ylabel("$\log_{10}(w)$")
	if save:
		Folder = "Plots_e_D_%s/" % D
		plt.savefig("../Plots/Plots_e/" + Folder + "Money Distribution_D_%d_N_%d_L_%.1f_a_%.1f" % (D,N,L,a) + ".png")
	else:
		plt.show()
	plt.close()	
	return

def Pareto_dist_e(M,W,filenames,save=False):
	D, N, L, a, g = extract_parametres(filenames[0])
	fig, ax = plt.subplots(2,2)	
	fig.suptitle("Linear Regression of log-log plot for $\omega$ \n N = %d, L = %.1f, $\\alpha$ = %.1f" % (N,L,a) )
	fig.text(0.5, 0.04, '$\log{10}(m)$', ha='center')
	fig.text(0.04, 0.5, '$\log_{10}(\omega)$', va='center', rotation='vertical')	


	SP = []
	for i in range(len(ax)):
		for j in range(len(ax[0])):
			SP.append(ax[i][j])
	
	for i in range(len(filenames)):
		D, N, L, a, g = extract_parametres(filenames[i])
		logM,logW = log_arrays(M[i],W[i])
		
		SI = 0 #start index for linreg
		mtemp = 10		
		eps = 0.1		
				
		if a == 1:		
			for j in range(len(logM)):
				if abs(logM[j]) < abs(mtemp):
					SI = j
					mtemp = logM[j]
			EI = len(logM)-500 # end index
		else:
			if L == 0:
				if g != 4:
					SI = 5
					EI = SI + 10 # end index
				else:
					SI = 14
					EI = SI + 20
			else:
				if g != 4:
					SI = 10
					EI = SI + 30
				else:
					SI = 30
					EI = SI + 10
		#print SI
		
		
		coeff = polyfit(logM[SI:EI],logW[SI:EI],1)
		Linreg = coeff[0]*logM + coeff[1]
		
		
		SP[i].plot(logM,logW,"o",markersize=0.8,label="$\omega (m ; \gamma = %.1f )$" % g)
		SP[i].plot(logM,Linreg,label="$ %.2f \cdot \log_{10}(m) + %.2f $" % (coeff[0],coeff[1]) )
		SP[i].legend()
		SP[i].set_xlim(-10,2)
		SP[i].set_ylim(-5,2)

	if save:
		Folder = "Plots_e_D_%s/" % D
		plt.savefig("../Plots/Plots_e/" + Folder + "Linear_Regression_Pareto_D_%d_N_%d_L_%.1f_%.1f" % (D,N,L,a) + ".png")
	else:
		plt.show()
	plt.close()	
	return


def Gibbs_dist_e(M,W,filenames,save=False):
	D, N, L, a, g = extract_parametres(filenames[0])
	fig, ax = plt.subplots(2,2)	
	fig.suptitle("Linear Regression of log plot for $\omega$ \n N = %d, L = %.1f" % (N,L) )
	fig.text(0.5, 0.04, '$(m)$', ha='center')
	fig.text(0.04, 0.5, '$\log(\omega)$', va='center', rotation='vertical')	


	SP = []
	for i in range(len(ax)):
		for j in range(len(ax[0])):
			SP.append(ax[i][j])
	
	for i in range(len(filenames)):
		D, N, L, a, g = extract_parametres(filenames[i])
		nlogM,nlogW = nlog_arrays(M[i],W[i])
		
		m = exp(nlogM)
		SI = 0 #start index for linreg
		mtemp = 10		
		eps = 0.1		
		'''		
		for j in range(len(logM)):
						
			if abs(logM[j]) < abs(mtemp):
				SI = j
				mtemp = logM[j]

		#print SI
		'''
		SI = 0
		EI = len(m)-500 # end index
		coeff = polyfit(m[SI:EI],nlogW[SI:EI],1)
		Linreg = coeff[0]*m + coeff[1]
		
		
		SP[i].plot(m,nlogW,label="$\omega (m ; \\alpha = %.1f )$" % a)
		SP[i].plot(m,Linreg,label="$ %.2f \cdot m + %.2f $" % (coeff[0],coeff[1]) )
		SP[i].legend()
		#SP[i].set_xlim(-1,2)
		#SP[i].set_ylim(-5,0.5)

	if save:
		Folder = "Plots_d_D_%s/" % D
		plt.savefig("../Plots/Plots_d/" + Folder + "Linear_Regression_Gibbs_D_%d_N_%d_L_%.1f" % (D,N,L) + ".png")
	else:
		plt.show()
	plt.close()	
	return

