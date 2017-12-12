import numpy as np
import matplotlib.pyplot as plt

params = {'text.usetex' : True,
          'font.size' : 8,
          'font.family' : 'lmodern',
          'text.latex.unicode' : True}
plt.rcParams['text.latex.preamble']=[r'\usepackage{lmodern}']
plt.rcParams.update(params)
colors = [ 'dimgrey', 'royalblue', 'crimson', 'seagreen', 'y' ]

from scipy.integrate import trapz, odeint

w = 5.33

if __name__ == '__main__':

    # --> Define the logistic map.
    def logistic_map(x, r):
        x = r * x * (1-x)
        return x

    # --> Population death.
    # R = [0.5, 1.5, 2.5, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9]
    # for i, r in enumerate(R):
    #     n = 100
    #     x0 = 0.1
    #     x = [x0]
    #
    #     for iteration in xrange(n):
    #         x.append(logistic_map(x[-1], r))
    #     x = np.asarray(x)
    #
    #     fig = plt.figure(figsize=(w, w/4))
    #     ax = fig.gca()
    #     ax.plot(x, '.-', color='royalblue')
    #     ax.set_ylim(-0.02, 1.02)
    #     ax.set_ylabel(r'$x_n$')
    #     ax.set_xlabel(r'$n$')
    #     ax.set_title(r'$\mu = %0.2f$' %r)
    #
    #     fig.savefig('../imgs/logistic_map_%i.pdf' %i, bbox_inches='tight', dpi=300)
    #     plt.close()

    # --> Bifurcation diagram.
    R = np.linspace(2.5, 3.9, 150)
    fig = plt.figure(figsize=(w, w/3))
    ax = fig.gca()

    for i, r in enumerate(R):
        n = 1000
        x0 = 0.1
        x = [x0]
        print r
        for iteration in xrange(n):
            x.append(logistic_map(x[-1], r))
        x = np.asarray(x)[900:]

        ax.plot(np.ones_like(x)*r, x, '.', color='royalblue', ms=1)
        ax.set_ylim(-0.02, 1.02)
        ax.set_ylabel(r'$x_n$')
        ax.set_xlabel(r'$\mu$')

        plt.savefig('../imgs/logistic_map_bifurcation.pdf', bbox_inches='tight', dpi=300)
    plt.show()
