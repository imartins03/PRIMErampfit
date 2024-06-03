
import numpy as np
import matplotlib.pyplot as plt
def f(x):
    return x**2

x = np.linspace(-10, 10, 100)
y = f(x)
plt.plot(x, y)
plt.show()
