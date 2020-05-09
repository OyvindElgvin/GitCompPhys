# The Ising Model and some of its thermodynamic properties.

This report studies the thermodynamic properties of the two dimensional Ising Model and demonstrates how a second order phase transition can be described using this binary model. The numerical approximation to the critical temperature of the phase transition is found to be close to the exact analytical solution, T = 2.269, in (Onsager 1994). The thermodynamic properties such as oscillating energy in time, the configuration of equilibrium, and temperature dependence are shown in various plots.

The report applies the Metropolis algorithm and shows how it can be used together with the Boltzmann distribution to determine whether or not to make a change of a thermodynamic variable (energy). This, combined with the random number generator **rand()** in *C++*, turns out to give satisfying results and simulate the intriguing phase transition from ferromagnetic to non-magnetic.

All results may be reproduced by running the executables under Results and plotting with **plot_data.py** and **T_plot**.