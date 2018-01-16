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

from pynamical import cobweb_plot, simulate, phase_diagram, bifurcation_plot
from numba import jit

@jit(nopython=True)
def sine_map(x, r):
    return r*np.sin(np.pi*x)

if __name__ == '__main__':

    pops = simulate(model=sine_map, num_gens=100, rate_min=0.6, rate_max=1, num_rates=1000, num_discard=100, jit=True)
    pops = pops.as_matrix()
    mu = np.linspace(0.6, 1, 1000)

    fig = plt.figure(figsize=(w, w/3))
    ax = fig.gca()

    for i in xrange(mu.size):
        ax.plot(mu[i]*np.ones_like(pops[:, i]), pops[:, i], ',', color='royalblue')


    # --> x-axis.
    ax.set_xlim(0.6, 1)
    ax.set_xlabel(r'$\mu$')

    # --> y-axis.
    ax.set_ylim(0, 1)
    ax.set_ylabel(r'$x^*$')

    plt.savefig('../imgs/sine_map_bifurcation_diagram.pdf', bbox_inches='tight', dpi=300)





    fig = plt.figure(figsize=(w/3, w/3))
    ax = fig.gca()

    # --> Logistic map.
    x = np.linspace(0, 1)
    f = 3.5*x*(1-x)
    ax.plot(x, f, color='royalblue', label=r'Log.')

    # --> Sine map.
    x = np.linspace(0, 1)
    f = 0.9*np.sin(np.pi*x)
    ax.plot(x, f, color='orange', label=r'Sin.')

    # --> Rossler.
    x = np.loadtxt('rossler_map.txt')
    ax.plot(x[:-1], x[1:], ls='', marker='.', ms=1, color='seagreen', label=r'R\"ossler')

    ax.set_xlim(0, 1)
    ax.set_xlabel(r'$x$')

    ax.set_ylim(0, 1)
    ax.set_ylabel(r'$f(x)$')

    ax.set_aspect('equal')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.33), ncol=2)

    plt.savefig('../imgs/unimodal_maps.pdf', bbox_inches='tight', dpi=300)
    plt.show()
