import numpy as np
from matplotlib import pyplot as plt


def curve(x):
    return 5*x**3

#mask

xx = np.linspace(-10, 10, 100)
y = curve(xx)

saturation_mask = y > 5
y[saturation_mask] = np.nan

x=xx
# x = np.asarray([i*np.ones_like(xx) for i in xx])

idx = np.isfinite(x) & np.isfinite(y)
fit = np.polyfit(xx[idx], y[idx], 1)
fit_fn = np.poly1d(fit)

plt.plot(xx,curve(xx))
plt.plot(xx, y)
plt.plot(xx,fit_fn(xx))
plt.show()




    # x = np.arange(len(y))  #???
    # y = np.asarray(y)
    # y = y.reshape((x.shape[0],-1))