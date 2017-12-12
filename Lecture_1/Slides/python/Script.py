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

    # --> Pedulum parameters.
    k, omega_0 = 0.1, 1.

    def pendulum(x, t):

        # --> Initialize variables.
        dx = np.zeros_like(x)

        # --> x-equation.
        dx[0] = x[1]

        # --> y-equation.
        dx[1] = -2*k*x[1] - omega_0**2 * np.sin(x[0])

        return dx

    # --> Phase plane.
    x = np.linspace(-2*np.pi, 2*np.pi, 64)
    y = np.linspace(-np.pi, np.pi, 64)
    x, y = np.meshgrid(x, y)

    u, v = pendulum([x[:], y[:]], 0.)

    # --> Plot the phase plane.
    fig = plt.figure(figsize=(w, w/3))
    ax = fig.gca()
    n = 2
    ax.quiver(x[::n, ::n], y[::n, ::n], u[::n, ::n], v[::n, ::n])

    xtickslabels = [r'$-2\pi$', r'$-\pi$', r'$0$', r'$\pi$', r'$2\pi$']
    xticks = [-2*np.pi, -np.pi, 0, np.pi, 2*np.pi]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtickslabels)
    ax.set_xlim(-2*np.pi, 2*np.pi)
    ax.set_xlabel(r'$x$')

    ytickslabels = [r'$-\pi$', r'$\displaystyle -\frac{\pi}{2}$', r'$0$', r'$\displaystyle \frac{\pi}{2}$', r'$\pi$']
    yticks = [-np.pi, -0.5*np.pi, 0, 0.5*np.pi, np.pi]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytickslabels)
    ax.set_ylim(-np.pi, np.pi)
    ax.set_ylabel(r'$y$', rotation=0)

    plt.savefig('../imgs/pendulum_phase_plane.pdf', bbox_inches='tight', dpi=300)




    # --> Fixed points.
    x_bar = np.arange(-2, 3)*np.pi
    y_bar = np.zeros_like(x_bar)

    ax.plot(x_bar, y_bar, 'ro', clip_on=False)

    plt.savefig('../imgs/pendulum_fixed_points.pdf', bbox_inches='tight', dpi=300)





    #--> Trajectoires.
    t = np.linspace(0, 50, 1000)
    for v_init in [0.5, 1, 1.5, 2, 2.5]:
        x0 = [0, v_init]
        x = odeint(pendulum, x0, t)
        ax.plot(x[:, 0], x[:, 1], color='royalblue')

    plt.savefig('../imgs/pendulum_fixed_trajectories_bis.pdf', bbox_inches='tight', dpi=300)

    plt.show()

    # --> Trajectoires (again).

    v_init = [0.5, 1, 1.5, 2, 2.5]
    for i, v in enumerate(v_init):
        x0 = [0, v]
        x = odeint(pendulum, x0, t)

        fig = plt.figure(figsize=(w, w/3))
        ax = fig.gca()
        ax.plot(t, x[:, 0], label=r'$\theta$')
        ax.plot(t, x[:, 1], label=r'$\dot{\theta}$')
        ax.legend(loc=0)
        ax.set_xlabel(r'$t$')

        plt.savefig('../imgs/pendulum_fixed_trajectories_%i.pdf' %i, bbox_inches='tight', dpi=300)

    plt.show()
