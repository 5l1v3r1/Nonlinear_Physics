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

def discrete_system(x, params):

    # --> Unpack parameters.
    r, a, b, c = params

    # --> Initialize variables.
    xn = np.zeros_like(x)

    # --> Prey's evolution.
    xn[0] = r*x[0]*(1.-x[0]) - a*x[0]*x[1]

    # --> Predator's evolution.
    xn[1] = c*x[1] + b*x[0]*x[1]

    return xn

if __name__ == '__main__':

    # --> Setup the parameters.
    r, a, b, c = 3.9, 2.3, 3.3, 0.1
    params = [r, a, b, c]

    # -->
    ngen = 10000
    x = np.zeros((ngen, 2))
    x[0, :] = 0.5, 0.1

    for i in xrange(1, ngen):
        x[i, :] = discrete_system(x[i-1, :], params)


    fig = plt.figure(figsize=(w/3, w/3))
    ax = fig.gca()

    ax.plot(x[-ngen/2:, 0], x[-ngen/2:, 1], ',', color='royalblue')

    ax.set_xlim(0, 1)
    ax.set_xlabel(r'Preys')

    ax.set_ylim(0, 1.1*x[:, 1].max())
    ax.set_ylabel(r'Predators')

    plt.savefig('../imgs/lotka_volterra.pdf', bbox_inches='tight', dpi=300)
    plt.show()
