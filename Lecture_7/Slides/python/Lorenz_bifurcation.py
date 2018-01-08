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

def rossler_system(x, t, params):

    # --> Unpack parameters.
    a, b, c = params

    # --> Initialize variables.
    dx = np.zeros_like(x)

    # --> x-equation.
    dx[0] = -x[1] - x[2]
    # --> y-equation.
    dx[1] = x[0] + a*x[1]
    # --> z-equation.
    dx[2] = b + x[2]*(x[0]-c)

    return dx

if __name__ == '__main__':

    # --> Miscellaneous.
    from scipy.signal import find_peaks_cwt
    from scipy.integrate import odeint
    from detect_peak import *


    C = np.linspace(1, 10, 1001)

    fig = plt.figure(figsize=(w, w/3))
    ax = fig.gca()

    for c in C:
        print c

        if c < 2:
            x0 = np.array([0., 1., 0.])

        params = [0.2, 0.2, c]
        t = np.linspace(0, 1000, 10000)
        x = odeint(rossler_system, x0, t, args=(params,))
        x0 = x[-1, :]
        x = x[-5000:, 0]
        peakind = detect_peaks(x)
        if peakind.size == 0:
            z = x[-1]
        else:
            z = x[peakind]
        ax.plot(c*np.ones_like(z), z, ',', color='royalblue')

    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$x_{\max}$')

    plt.savefig('../imgs/rossler_bifurcation_diagram.pdf', bbox_inches='tight', dpi=300)





    C = np.linspace(2.5, 4.5, 1001)

    fig = plt.figure(figsize=(w, w/3))
    ax = fig.gca()

    for c in C:
        print c

        if c < 2:
            x0 = np.array([0., 1., 0.])

        params = [0.2, 0.2, c]
        t = np.linspace(0, 1000, 10000)
        x = odeint(rossler_system, x0, t, args=(params,))
        x0 = x[-1, :]
        x = x[-5000:, 0]
        peakind = detect_peaks(x)
        if peakind.size == 0:
            z = x[-1]
        else:
            z = x[peakind]
        ax.plot(c*np.ones_like(z), z, ',', color='royalblue')

    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$x_{\max}$')

    plt.savefig('../imgs/rossler_bifurcation_diagram_zoom.pdf', bbox_inches='tight', dpi=300)


    C = np.linspace(4.05, 4.25, 1001)

    fig = plt.figure(figsize=(w, w/3))
    ax = fig.gca()

    for c in C:
        print c

        if c < 2:
            x0 = np.array([0., 1., 0.])

        params = [0.2, 0.2, c]
        t = np.linspace(0, 1000, 10000)
        x = odeint(rossler_system, x0, t, args=(params,))
        x0 = x[-1, :]
        x = x[-5000:, 0]
        peakind = detect_peaks(x)
        if peakind.size == 0:
            z = x[-1]
        else:
            z = x[peakind]
        ax.plot(c*np.ones_like(z), z, ',', color='royalblue')

    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$x_{\max}$')
    ax.set_ylim(4.5, 5.5)

    plt.savefig('../imgs/rossler_bifurcation_diagram_zoom_bis.pdf', bbox_inches='tight', dpi=300)


    #####

    params = [0.2, 0.2, 5]
    t = np.linspace(0, 10000, 1000000)
    x0 = np.array([0., 1., 0.])
    x = odeint(rossler_system, x0, t, args=(params,))
    t = t[-500000:]
    x = x[-500000:, 0]
    peakind = detect_peaks(x)
    if peakind.size == 0:
        z = x[-1]
    else:
        z = x[peakind]

    fig = plt.figure(figsize=(w, w/3))
    ax = fig.gca()
    ax.plot(t, x, color='royalblue')
    ax.plot(t[peakind], x[peakind], 'o', color='orange', mfc='None')
    ax.set_xlim(9500, 10000)
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$x$')

    plt.savefig('../imgs/rossler_successive_maxima.pdf', bbox_inches='tight', dpi=300)







    from sklearn.preprocessing import PolynomialFeatures
    A = np.array([np.ones_like(z[:-1]), z**2, z**4, z**6, z**8]).T
    from scipy.linalg import lstsq
    coef = lstsq(A, z[1:])[0]
    print coef

    fig = plt.figure(figsize=(w/3, w/3))
    ax = fig.gca()

    ax.plot(z[:-1], z[1:], '.', color='royalblue')

    xx = np.linspace(0, 12)
    ax.plot(xx, xx, color='lightgray')

    ax.set_xlim(0, 12)
    ax.set_ylim(0, 12)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x_{k}$')
    ax.set_ylabel(r'$x_{k+1}$')

    plt.savefig('../imgs/rossler_lorenz_map.pdf', bbox_inches='tight', dpi=300)

    fig = plt.figure(figsize=(w/3, w/3))
    ax = fig.gca()

    ax.plot(z[:-1], z[1:], '.', color='royalblue')

    xx = np.linspace(0, 12)
    ax.plot(xx, xx, color='lightgray')
    A = library.transform(xx.reshape(-1, 1))
    ax.plot(xx, A.dot(coef), color='gray', ls='--')

    ax.set_xlim(0, 12)
    ax.set_ylim(0, 12)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x_{k}$')
    ax.set_ylabel(r'$x_{k+1}$')

    plt.savefig('../imgs/rossler_lorenz_map_bis.pdf', bbox_inches='tight', dpi=300)

    np.savetxt('rossler_map.txt', z/12.)




    fig, axes = plt.subplots(1, 4, figsize=(w, w/3), sharex=True, sharey=True)
    C = [2.5, 3.5, 4, 5]

    for i, c in enumerate(C):
        params = [0.2, 0.2, c]
        t = np.linspace(0, 10000, 1000000)
        x0 = np.array([0., 1., 0.])
        x = odeint(rossler_system, x0, t, args=(params,))
        t = t[-10000:]
        x = x[-10000:, :]

        axes[i].plot(x[:, 0], x[:, 1], color='royalblue', lw=1)

        axes[i].set_xlim(-10, 10)
        axes[i].set_ylim(-10, 10)

        axes[i].set_xlabel(r'$x$')
        if i == 0:
            axes[i].set_ylabel(r'$y$')
        axes[i].set_title(r'$c = %.1f$' %c)

    plt.savefig('../imgs/rossler_system_period_doubling.pdf', bbox_inches='tight', dpi=300)
    plt.show()
