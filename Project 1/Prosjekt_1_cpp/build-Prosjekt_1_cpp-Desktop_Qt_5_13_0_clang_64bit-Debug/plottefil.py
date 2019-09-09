
from matplotlib.pyplot import *
import numpy as np

ns= [10, 1e2, 1e3]#, 1e4, 1e5, 1e6, 1e7]

hs = []
errors = []

for n in ns:
    with open("n={}".format(int(n)), "r") as f:
        f.readline()
        x = []
        general = []
        spes = []
        exact = []
        error_general = []
        error_spes = []

        for line in f:
            words = line.split()
            # x
            xword1 = float(words[0][:-1])
            x.append(xword1)

            # general
            yword1 = float(words[1][:-1])
            general.append(yword1)

            # spes
            yword2 = float(words[2][:-1])
            spes.append(yword2)

            # exact
            yword3 = float(words[3][:-1])
            exact.append(yword3)

            # error_general
            yword4 = float(words[4][:-1])
            error_general.append(yword4)

            # error_spes
            try:
                yword5 = float(words[5][:-1])
                error_spes.append(yword5)
            except:
                error_spes.append(np.nan)

        # h
        h = 1/(1+n)
        hs.append(h)
        #print(n, h)


        errors.append(min(error_general[1:-1]))
        #print(errors)


        #plot(x, general,label="Generalized")
        #plot(x, spes,label="Specialized")
        #plot(x, exact,label="Exact")
        #title("Generalized solution with n = %g" %(n))
        #legend()
        #xlabel("x = i*h")
        #ylabel("v")
        #show()

        plot(x, general, label="n = %g" %(n))
        title("Comparing generalized Thomas for various n")
        xlabel("x")
        ylabel("v")
        legend()
        f.close()


#plot(np.log10(np.array(hs)),errors)
#title("log10 Relative Error")
#xlabel("log10(h)")
#ylabel("y")

#plot(x, general, label="n = %g" %n)



show()
