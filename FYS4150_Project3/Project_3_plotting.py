from numpy import zeros, ones, array, sqrt, linspace, exp, meshgrid, pi, log10, polyfit, polyval
import matplotlib.pyplot as plt
from Py_Functions import readarrays

Exact = (5*pi**2)/256.


# Plotting integrand as function of r1 and r2

yn = raw_input("Do you want to plot integrand? y/n \n")
if yn == "y":
	def integrand(r1,r2):
		f = exp(-4*(r1+r2))
		return f

	r1 = linspace(0,0.5,1000)
	r2 = linspace(0,0.5,1000)

	X,Y = meshgrid(r1,r2)

	I = integrand(X,Y) 
	plt.figure()
	plt.title("Plot of $ e^{-4(r_1+r_2)}$ as function of $r_1$ and $r_2$",y=1.02)
	plt.xlabel("$r_1$")
	plt.ylabel("$r_2$")
	plt.contourf(X,Y,I)
	plt.colorbar()
	plt.savefig("integrand_plot.png")
	plt.show()


#Tabulating results from Laguerre

yn = raw_input("Do you want to tabulate Laguerre results? y/n \n")
if yn == "y":
	Values = open("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/Results_Laguerre.txt")

	lines = Values.readlines()

	Values.close()

	N_values = []
	I_values = []

	for i in range(len(lines)):
		line_i = lines[i].split()
		N_values.append(float(line_i[2]))
		I_values.append(float(line_i[5]))

	Runtimes = readarrays("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/Gauss_Laguerre_Runtimes.txt")[0][1]

	N_values = array(N_values)
	I_values = array(I_values)
	Error_Laguerre = abs(I_values-ones(len(I_values))*Exact)

	Lag_error = open("Laguerre_Errors.txt","w+")
	Lag_error.write("Results of Integral using Gaussian Laguerre \n\n")
	Lag_error.write("          |  Computed Value  |   Error   |      Runtime     |\n")
	for i in range(len(N_values)):
		Lag_error.write(" N = %3d  |     %8.5f     |  %7.5f  |   %10.5f     |\n" % (N_values[i], I_values[i],Error_Laguerre[i],Runtimes[i]))

	Lag_error.close()


#Tabulating results from Legendre

yn = raw_input("Do you want to tabulate Legendre results? y/n \n")
if yn == "y":

	Values = open("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/Results_Legendre.txt")

	lines = Values.readlines()

	Values.close()

	L_values = []
	N_values = []
	I_values = []
	Error_Legendre = []

	N_values.append([])
	I_values.append([])
	Error_Legendre.append([])
	L_counter = 1
	N_counter = 0

	L_values.append(float(lines[0].split()[4]))
	for i in range(1,len(lines)):
		if L_counter == 1:
			N_counter += 1
		if lines[i].split()[4] != lines[i-1].split()[4]:
			N_values.append([])
			I_values.append([])
			Error_Legendre.append([])

			L_values.append(float(lines[i].split()[4]))
			L_counter += 1

	i = 0

	for l in range(L_counter):
		for n in range(N_counter):
			N_values[l].append(int(lines[i].split()[7]))
			I_values[l].append(float(lines[i].split()[10]))
			Error_Legendre[l].append(abs(float(lines[i].split()[10])-Exact))
			i += 1

	Leg_error = open("Legendre_Errors.txt","w+")
	Leg_error.write("Results of Integral using Gaussian Legendre \n\n")
	Leg_error.write("                      |  Computed Value  |   Error   |\n")
	for l in range(L_counter):
		Leg_error.write(" L = %d  |" % L_values[l])
		for i in range(N_counter):
			if i == 0:
				Leg_error.write("   N = %3d   |     %8.5f     |  %7.5f  |\n" % (N_values[l][i], I_values[l][i],Error_Legendre[l][i]))		
			else:		
				Leg_error.write("        |   N = %3d   |     %8.5f     |  %7.5f  |\n" % (N_values[l][i], I_values[l][i],Error_Legendre[l][i]))

	Leg_error.close()

#Plot runtimes for Legendre

yn = raw_input("Do you want to plot runtimes? y/n \n")
if yn == "y":
	

	Runtime_Files = ["Gauss_Legendre_Runtimes", "Gauss_Laguerre_Runtimes"]
	

	Legends = []
	plt.figure()
	plt.title("Runtimes for Gaussian Legendre and Gaussian Laguerre \n as function of number of data points, N")
	plt.xlabel("$log_{10}(N)$")
	plt.ylabel("$log_{10}(t_R)$")	

	Colours = ["b","g","ro","mo"]
	for i in range(len(Runtime_Files)):
	
		N, Runtimes = readarrays("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/%s.txt" % Runtime_Files[i])[0]
	
		N = array(N); Runtimes = array(Runtimes)
		lN = log10(N)
		lR = log10(Runtimes)

		lininterpol = polyfit(lN,lR,1)

		N_array = linspace(0,2*max(lN),2000)
		linpolynomial = lininterpol[0]*N_array + lininterpol[1]


		plt.plot(N_array,linpolynomial,"%s" % Colours[i],lN,lR,"%s" % Colours[i+2],)

		Legends.append("%.2f$x$ + %.2f" % (lininterpol[0], lininterpol[1]))		
		Legends.append("%s, calculated" % Runtime_Files[i])
		
	
	plt.legend(Legends)
	plt.savefig("Runtime plots Gaussian Quadrature.png")
	plt.show()

#Tabulating runtimes for Gaussian Quadrature

yn = raw_input("Do you want to tabulate Gaussian runtimes? y/n \n")
if yn == "y":
	
	runtime_table = open("Gaussian_runtimes.txt","w+")

	legN, legRuntimes = readarrays("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/Gauss_Legendre_Runtimes.txt")[0]
	
	lagN, lagRuntimes = readarrays("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/Gauss_Laguerre_Runtimes.txt")[0]

	runtime_table.write("  N | Legendre | Laguerre |\n")
	runtime_table.write("    | Runtimes | Runtimes |\n")
	runtime_table.write("    |          |          |\n")
	for i in range(len(lagN)):
		runtime_table.write(" %d |  %6.2f  |  %7.2f |\n" % (legN[i],legRuntimes[i],lagRuntimes[i])) 


	runtime_table.close()


#Tabulating results from Brute Force Monte Carlo

yn = raw_input("Do you want to tabulate Brute Force Monte Carlo results? y/n \n")
if yn == "y":
	Values = open("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/Results_BFMC.txt")

	lines = Values.readlines()

	Values.close()

	N_values = []
	I_values = []
	V_values = []
	E_values = []

	for i in range(len(lines)):
		line_i = lines[i].split()
		N_values.append(float(line_i[7]))
		I_values.append(float(line_i[10]))
		V_values.append(float(line_i[13]))
		E_values.append(float(line_i[16]))

	print N_values
	print I_values
	print V_values
	
	Runtimes = readarrays("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/BFMC_Runtimes.txt")[0][1]

	N_values = array(N_values)
	I_values = array(I_values)

	BFMC_error = open("BFMC_Errors.txt","w+")
	BFMC_error.write("Results of Integral using Brute Force Monte Carlo \n\n")
	BFMC_error.write("            |  Computed Value  |  Average Error  |  Average Variance  |  Average Runtime  |\n")
	for i in range(len(N_values)):
		BFMC_error.write(" N = %5.0e  |     %8.5f     |    %9.7f    |     %.4e     |     %10.5f    |\n" % (N_values[i], I_values[i],E_values[i],V_values[i], Runtimes[i]))

	BFMC_error.close()
	
#Plotting Errors from Brute Force Monte Carlo

yn = raw_input("Do you want to plot Brute Force Monte Carlo errors? y/n \n")
if yn == "y":
	Values = open("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/Results_BFMC.txt")

	lines = Values.readlines()

	Values.close()

	N_values = []
	I_values = []
	V_values = []

	for i in range(len(lines)):
		line_i = lines[i].split()
		N_values.append(float(line_i[7]))
		I_values.append(float(line_i[10]))
		V_values.append(float(line_i[13]))

	N_values = array(N_values)	
	I_values = array(I_values)
	Error_BFMC = abs(I_values-ones(len(I_values))*Exact)
	
	plt.figure()
	plt.plot(log10(N_values),log10(Error_BFMC),'b',log10(N_values),log10(Error_BFMC),'bo')
	plt.show()


#Tabulating results from Importance Sampling Monte Carlo

yn = raw_input("Do you want to tabulate Importance Sampling Monte Carlo results? y/n \n")
if yn == "y":
	Values = open("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/Results_ISMC.txt")

	lines = Values.readlines()

	Values.close()

	N_values = []
	I_values = []
	V_values = []
	E_values = []

	for i in range(len(lines)):
		line_i = lines[i].split()
		N_values.append(float(line_i[2]))
		I_values.append(float(line_i[5]))
		V_values.append(float(line_i[8]))
		E_values.append(float(line_i[11]))

	print N_values
	print I_values
	print V_values
	
	Runtimes = readarrays("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/ISMC_Runtimes.txt")[0][1]

	N_values = array(N_values)
	I_values = array(I_values)

	ISMC_error = open("ISMC_Errors.txt","w+")
	ISMC_error.write("Results of Integral using Importance Sampling Monte Carlo \n\n")
	ISMC_error.write("            |  Computed Value  |  Average Error  |  Average Variance  |  Average Runtime  |\n")
	for i in range(len(N_values)):
		ISMC_error.write(" N = %5.0e  |     %8.5f     |    %9.7f    |     %.4e     |     %10.5f    |\n" % (N_values[i], I_values[i],E_values[i],V_values[i], Runtimes[i]))

	ISMC_error.close()


#Tabulating runtimes for Importance Sampling Monte Carlo

yn = raw_input("Do you want to tabulate Importance Sampling Monte Carlo runtimes? y/n \n")
if yn == "y":

	runtime_table = open("ISMC_Runtimes.txt","w+")

	npN, npRuntimes = readarrays("build-Project_3-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/ISMC_Runtimes.txt")[0]
	
	pN, pRuntimes = readarrays("Project_3_MC/ISMC_Runtimes_para.txt")[0]

	runtime_table.write("     N     | No parallelization | Parallelization | pR / npR |\n")
	runtime_table.write("           |   Runtimes, npR    |  Runtimes, pR   |          |\n")
	runtime_table.write("           |                    |                 |          |\n")
	for i in range(2,len(npN)):
		runtime_table.write(" %9d |        %6.2f      |     %7.2f     |   %4.2f   |\n" % (npN[i],npRuntimes[i],pRuntimes[i],pRuntimes[i]/npRuntimes[i])) 


	runtime_table.close()


