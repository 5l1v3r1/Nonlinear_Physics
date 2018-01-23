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

# --> Import FFT package from scipy.
from scipy.fftpack import fft, ifft
# --> Convolution operator.
from scipy.signal import convolve

class KS:

    def __init__(self, nx=2**9, lx=32*np.pi, t=np.linspace(0, 500, 500)):

        ###########################################
        #####     USER-DEFINED QUANTITIES     #####
        ###########################################

        # --> Streamwise extent of the computational domain.
        self.lx = lx
        # --> Number of grid points in the x-direction.
        self.nx = nx
        # --> Time integration.
        self.t = t

        ######################################
        #####     DERIVED QUANTITIES     #####
        ######################################

        # --> Grid-spacing in the x-direction.
        self.dx = 1.0*lx/nx
        # --> Define the grid in the x-direction for computation purposes.
        x = np.linspace(0, lx, nx+1) - lx/2
        self.x = x[:-1] # NOTE: Lat point is removed for coding (periodicity implied).
        # --> Define the Fourier wavenumbers in the x-direction.
        self.kx = (2*np.pi/lx) * np.concatenate((np.arange(0, nx/2 + 1), np.arange(-nx/2+1, 0)))

    def dynamical_system(self, t, u):

        '''
        This function implements the evaluation of the right-hand side of the
        Kuramoto-Sivashinsky equation in Fourier space. Note that the non-linear
        term is evaluated in physical space and that no anti-aliasing technique
        is used.
        '''

        # --> Compute the linear terms.
        rhs = ( self.kx**2 - self.kx**4 )*u
        # --> Evaluate non-linear term.
        rhs += 1j*self.kx*fft( (ifft(u).real)**2 ) * 0.5

        return rhs

    def solve(self, u0):

        '''
        This function solves the Kuramoto-Sivashinsky equation for the parameters
        used when defining the KS object.
        '''

        # --> Import the Initial Value Problem solver from scipy.
        from scipy.integrate import solve_ivp
        # --> Solve the equation from t=0 to t=t.max() using a 4th order accurate
        #     Runger-Kutta method.
        output = solve_ivp(self.dynamical_system, [0., self.t.max()], u0, t_eval=self.t, method='RK45')
        # --> Save the spatio-temporal evolution of the state vector in Fourier space.
        self.u_hat = output['y']
        # --> Save the spatio-temporal evolution of the state vector in physical space.
        self.u = ifft(output['y'], axis=0).real

        return self

if __name__ == '__main__':


    ######################################################################
    #####                                                            #####
    ##### DIRECT NUMERICAL SIMULATION OF THE KURAMOTO-SIVASHINKY EQ. #####
    #####                                                            #####
    ######################################################################

    # --> Setup the problem.
    ks = KS(nx=48, t=np.linspace(0, 2000, 8000), lx=20*np.pi)

    # --> Random initial condition.
    from numpy.random import random
    u0 = random(size=ks.nx)
    u0 = fft( u0 - u0.mean() ) # NOTE: Make sure the initial condition has zero-mean.

    # --> Integrate in time.
    ks.solve(u0)

    # --> Plot the spatio-temporal dynamics.
    fig = plt.figure(figsize=(w, w/3))
    ax = fig.gca()
    ax.contour(ks.t, ks.x, ks.u, cmap=plt.cm.RdBu)

    ax.set_ylim(-ks.lx/2, ks.lx/2)
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$x$')


    fig = plt.figure(figsize=(w, w/3))
    ax = fig.gca()

    ax.loglog(ks.kx[ks.kx>0], abs(ks.u_hat[ks.kx>0]).mean(axis=1))

    plt.show()
