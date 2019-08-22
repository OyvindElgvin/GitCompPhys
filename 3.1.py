
import numpy as np
import matplotlib.pyplot as plt

"""
1
"""

x = np.sqrt(2)

def f(x):
    return np.arctan(x)

def O(h,x):
    return np.abs(1 / (x**2 + 1) - (f(x+h) - f(x)) / h)


h1 = np.linspace(1e-1,1e-14,10000)

o_liste = O(h1,x)

listeindex = np.argmin(o_liste)
#finner indexen med minste verdi

print("minste feil er: ", o_liste[listeindex],"med en h = ", h1[listeindex])


"""
1x-193-157-176-229:GitCompPhys oyvindengebretsen$ python 3.1.py
minste feil er:  1.57149536345e-06 med en h =  1.000100011e-05
"""


"""
2
"""

h2 = np.logspace(-12, -1, 10000)
# logspace gir lista [1e-1, 1e-2, 1e-3, ...]

# Topunktsmetoden: Forward Euler
def df(x,h):
    return (f(x+h) - f(x)) / h

# Trepunktsmetoden:
def ddf(x,h):
    return (f(x+h) - f(x-h)) / (2*h)


plt.title("Verdien")
plt.plot(h2, df(x,h2))
plt.plot(h2, ddf(x,h2))
plt.axhline(y=1.0/3)
plt.show()

plt.title("Feilen")
plt.loglog(h2, np.abs(df(x,h2)-1.0/3))
plt.loglog(h2, np.abs(ddf(x,h2)-1.0/3))
plt.axhline(y=1.0/3)
plt.show()


"""
3
"""
