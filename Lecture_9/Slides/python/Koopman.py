import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec
from numpy.random import uniform, normal

params = {'text.usetex' : True,
          'font.size' : 8,
          'font.family' : 'lmodern',
          'text.latex.unicode' : True}
plt.rcParams['text.latex.preamble']=[r'\usepackage{lmodern}']
plt.rcParams['figure.facecolor'] = '1'
plt.rcParams.update(params)

w = 5.33

if __name__ == '__main__':

    lambda_, mu_ = 0.1, 0.9

    def dynamical_system(x):

        # --> Initialize variables.
        y = np.zeros_like(x)

        # --> Initialize variables.
        y[0] = lambda_ * x[0]
        y[1] = mu_ * x[1] + (lambda_**2 - mu_) * x[0]**2

        return y

    x = np.linspace(-2, 2, 20)
    x, y = np.meshgrid(x, x)

    a = np.zeros_like(x)
    b = np.zeros_like(y)

    a[:], b[:] = dynamical_system([x[:], y[:]])

    fig = plt.figure(figsize=(w/2, w/2))
    ax = fig.gca()

    ax.plot(a, b, '.')

    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)

    plt.show()

    plt.show()
