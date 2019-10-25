
import numpy as np
import matplotlib.pyplot as plt
from plotly import offline as py
py.init_notebook_mode()

t = np.linspace(0, 20, 500)
plt.plot(t, np.sin(t))

py.iplot_mpl(plt.gcf())




print("konge")
