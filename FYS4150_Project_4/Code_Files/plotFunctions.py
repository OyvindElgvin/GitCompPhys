from numpy import *
from matplotlib.pyplot import *
from numba import jit

def plotMostLikelyState(lst,temp,n):

    T = temp
    h = int(n/len(lst[0][1][0])) #Step size
    n_values = linspace(h,n+h,int(n/h))
    fig1, axs1 = subplots(4,1,sharex=True,gridspec_kw={'hspace': 0})
    fig1.text(0.04, 0.5, "$\\langle E\\rangle$", va='center', rotation='vertical')
    fig2, axs2 = subplots(4,1,sharex=True,gridspec_kw={'hspace': 0})
    fig2.text(0.04, 0.5, "$\\langle |M|\\rangle$", va='center', rotation='vertical')
    state = 1
    plt = 0
    for i in lst:
        initial = 0
        if state == 1:
            initial = "Ordered"
        else:
            initial = "Random"
        index = 0
        for j in i[1]:
            E_mean = j
            axs1[plt].plot(log10(n_values),E_mean,label="T="+str(T[index][0])+", "+initial)
            #axs1[plt].set_ylabel("$\\langle E\\rangle$")
            axs1[plt].legend()
            plt += 1
            index += 1

        plt -= 2
        index = 0
        for j in i[3]:
            M_abs = j
            axs2[plt].plot(log10(n_values),M_abs,label="T="+str(T[index][0])+", "+initial)
            #axs2[plt].set_ylabel("$\\langle |M|\\rangle$")
            axs2[plt].legend()
            plt += 1
            index += 1

        state += 1
    axs1[-1].set_xlabel("$log_{10}$(# of MC-cycles)")
    axs2[-1].set_xlabel("$log_{10}$(# of MC-cycles)")
    fig1.savefig("../Figures/Most_Likely_State_E_mean_L_20.pdf")
    fig2.savefig("../Figures/Most_Likely_State_M_abs_L_20.pdf")
    show()

def plotAcceptedConfigurations(lst,temp,n):

    T = temp
    h = int(n/len(lst[0][6][0])) #Step size
    n_values = linspace(h,n+h,int(n/h))
    fig1, axs1 = subplots(4,1,sharex=True,gridspec_kw={'hspace': 0})
    fig1.text(0.04, 0.5, "Accepted configurations", va='center', rotation='vertical')
    state = 1
    plt = 0
    for i in lst:
        initial = 0
        if state == 1:
            initial = "Ordered"
        else:
            initial = "Random"
        index = 0
        for j in i[7]:
            Accepted_configs = j
            axs1[plt].plot(log10(n_values),Accepted_configs,label="T="+str(T[index][0])+", "+initial)
            #axs1[plt].set_ylabel("Acc. conf.")
            axs1[plt].legend()
            plt += 1
            index += 1

        state += 1
    axs1[-1].set_xlabel("# of MC-cycles")
    fig1.savefig("../Figures/Number_of_Accepted_Configs_L_20.pdf")
    show()

def plotProbabilityDistribution(energies,temp,variances,mean_energies,n):
    T = temp
    fig, axs = subplots(2,1,sharex=True,gridspec_kw={'hspace': 0})
    fig.text(0.04, 0.5, "P(E)", va='center', rotation='vertical')
    for i in range(len(energies)):
        variance = variances[i][0]
        E_mean = mean_energies[i][0]
        col = int((max(energies[i])-min(energies[i]))/4)
        N, bins, patches = axs[i].hist(energies[i],col,normed=True,label="T = "+str(T[i][0])+"\n$\\sigma_{E}$ = "+str(sqrt(variance)))
        axs[i].legend()
        axs[i].axvline(x=E_mean-sqrt(variance),linewidth=1,linestyle="--",color='r')
        axs[i].axvline(x=E_mean+sqrt(variance),linewidth=1,linestyle="--", color='r')

        bin_width = bins[1] - bins[0]
        run = 0
        integral = 0
        for h in N:
            if (bins[run] >= E_mean-sqrt(variance) and bins[run] <= E_mean+sqrt(variance)):
                integral += bin_width*h
            run += 1
        print("Data set within first order STD (T = "+str(T[i][0])+"): "+str(integral*100)+"%")

    axs[-1].set_xlabel("Energy, E")
    fig.savefig("../Figures/Probability_Distribution_N_"+str(n)+"_L_20.pdf")
    show()

def plotPhaseTransition(lst,temp,n):
    L = [40,60,80,100]

    T = temp
    E = zeros(len(lst[0][0]))
    M = zeros(len(lst[0][0]))
    Cv = zeros(len(lst[0][0]))
    X = zeros(len(lst[0][0]))

    size = 0
    for i in lst: #Loop over all different L calculations
        print(size)
        index = 0
        for j in i[1]: # Loop over Energies per temperature in given LxL system
            E[index] = j[-1]
            index += 1
        figure(1)
        plot(T,E,label="L = "+str(L[size]))
        index = 0
        for j in i[3]: # Loop over Magnetizations per temperature
            M[index] = j[-1]
            index += 1
        figure(2)
        plot(T,M,label="L = "+str(L[size]))
        index = 0
        for j in i[4]: # Loop over Heat Capacity per temperature
            Cv[index] = j[-1]
            index += 1
        figure(3)
        plot(T,Cv,label="L = "+str(L[size]))
        index = 0
        for j in i[6]: # Loop over Susceptibility per temperature
            X[index] = j[-1]
            index += 1
        figure(4)
        plot(T,X,label="L = "+str(L[size]))
        size += 1
    figure(1)
    legend()
    xlabel("Temperature (T)")
    ylabel("$\\langle E\\rangle (T)$")
    savefig("../Figures/E_of_T_N_"+str(n)+".pdf")
    figure(2)
    legend()
    xlabel("Temperature (T)")
    ylabel("$\\langle |M|\\rangle (T)$")
    savefig("../Figures/M_of_T_N_"+str(n)+".pdf")
    figure(3)
    legend()
    xlabel("Temperature (T)")
    ylabel("$C_{V} (T)$")
    savefig("../Figures/Cv_of_T_N_"+str(n)+".pdf")
    figure(4)
    legend()
    xlabel("Temperature (T)")
    ylabel("$\\chi_{abs} (T)$")
    savefig("../Figures/X_of_T_N_"+str(n)+".pdf")
    show()
