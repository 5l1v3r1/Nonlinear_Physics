# --> Import some functions for Python 2 / Python 3 compatibility.
from __future__ import print_function, division

# --> Import general purpose python libraries.
import numpy as np
import matplotlib.pyplot as plt

params = {'text.usetex' : True,
          'font.size' : 8,
          'font.family' : 'lmodern',
          'text.latex.unicode' : True}
plt.rcParams['text.latex.preamble']=[r'\usepackage{lmodern}']
plt.rcParams['figure.facecolor'] = '1'
plt.rcParams.update(params)

# --> Import pandas to save/read time-series data from disk.
import pandas as pd

# --> Import FENiCS: Finite Element library.
from dolfin import *

# --> Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# --> Define stress tensor
def sigma(u, p, nu):
    return 2*nu*epsilon(u) - p*Identity(len(u))

# --> Define the Pinball class.
class Pinball:





    def __init__(self, re=100):

        # --> Reynolds number.
        self.re = re

        # --> Load the Pinball mesh.
        mesh = Mesh('pinball_mesh.xml')
        self.mesh = mesh




    def dns(self,
            iostep=200,
            dt=0.01,
            T=100,
            initial_condition=None,
            figname='pinball',
            avg=False,
            control=None
            ):

        """
        This function implements a finite-element solver for the unsteady
        Navier-Stokes equation (Pinball configuration). The temporal scheme
        implemented uses a second-order accurate backward-differentation for
        the linear terms and second-order accurate Adam-Bashforth for the
        nonlinear ones.

        Inputs:
        ------

        iostep : (default 200) A plot of the instantaneous vorticity field and
                 a pinbal_00000.h5 file are outpost every iostep.

        dt : (default 0.005) Constant time-step for the temporal integration.

        T : (default T) Final time for the temporal integration.

        initial_condition : (default None) Specify the HDF5 file containing the
                            initial condition.

        figname : (default pinball) Standard filename for the figures.

        avg : (default False) Boolean deciding whether the temporal average
              of the flow has to be computed or not.

        control : (default None) Function handle to compute the control law
                  given the sensor measurements.
        """

        # --> Parameters.
        self.dt = dt
        k = Constant(dt)
        nsteps = int(T/dt)

        #--------------------------------
        #-----     Mesh-related     -----
        #--------------------------------

        print( 'Setting-up the boundary conditions...' )

        # --> Inflow boundary.
        inflow = 'near(x[0], -6)'
        # --> Outflow boundary.
        outflow = 'near(x[0], 20)'
        # --> Top and bottom boundaries.
        walls = 'near(x[1], -6) || near(x[1], 6)'
        # --> Front cylinder.
        front_cylinder = 'on_boundary && x[0]>-2. && x[0]<-0.75 && x[1]>-0.6 && x[1]<0.6'
        # --> Bottom cylinder.
        bottom_cylinder = 'on_boundary && x[0]>-0.6 && x[0]<0.6 && x[1]>-1.5 && x[1]<0.'
        # --> Top cylinder.
        top_cylinder = 'on_boundary && x[0]>-0.6 && x[0]<0.6 && x[1]>0. && x[1]<1.5'


        # --> Inflow velocity profile.
        inflow_profile = ('1', '0')

        #--------------------------------------
        #-----     FE Function spaces     -----
        #--------------------------------------

        print( 'Initializing the function spaces...' )

        # --> Define the velocity and pressure spaces.
        V = VectorFunctionSpace(self.mesh, 'CG', 2)
        Q = FunctionSpace(self.mesh, 'CG', 1)

        # --> Define trial and test functions.
        u, v = TrialFunction(V), TestFunction(V)
        p, q = TrialFunction(Q), TestFunction(Q)

        #---------------------------------------
        #-----     Boundary Conditions     -----
        #---------------------------------------

        print( 'Setting-up the boundary conditions bis...' )

        # --> Sets inflow boundary condition.
        bc_inflow = DirichletBC(V, inflow_profile, inflow)
        # --> Sets boundary conditions on the lateral planes.
        bc_walls = DirichletBC(V, inflow_profile, walls)
        # --> Set boundary condition on front cylinder.
        front_control = Expression(('-0.5*omega*x[1]', '0.5*omega*(x[0]+3*sqrt(3)/4)'), omega=0, degree=2)
        bc_front_cyl = DirichletBC(V, front_control, front_cylinder)
        # --> Set boundary condition on bottom cylinder.
        bottom_control = Expression(('-0.5*omega*(x[1]+0.75)', '0.5*omega*x[0]'), omega=0, degree=2)
        bc_bottom_cyl = DirichletBC(V, bottom_control, bottom_cylinder)
        # --> Set boundary condition on top cylinder.
        top_control = Expression(('-0.5*omega*(x[1]-0.75)', '0.5*omega*x[0]'), omega=0, degree=2)
        bc_top_cyl = DirichletBC(V, top_control, top_cylinder)
        # --> Set pressure boundary condition at the outflow.
        bc_outflow = DirichletBC(Q, Constant(0), outflow)

        # --> Store boundary conditions in a python list.
        bcu = [bc_inflow, bc_walls, bc_front_cyl, bc_bottom_cyl, bc_top_cyl]
        bcp = [bc_outflow]

        #---------------------------------------
        #-----     Variational problem     -----
        #---------------------------------------

        # --> Define functions for solutions at previous and current time steps
        u_nn, u_n, u_ = Function(V), Function(V), Function(V)
        p_n, p_ = Function(Q), Function(Q)

        # --> Load initial condition if needed.
        if initial_condition is not None:
            print( 'Loading initial condition...' )
            u_n, p_n = self.load_solution(initial_condition, up=[u_n, p_n])
            u_nn.assign(u_n)

        # --> Initialize arrays for time-averaging.
        if avg is True:
            u_avg = Function(V, name='mean_flow')
            u_avg.vector()[:] = u_n.vector()[:]
            p_avg = Function(Q, name='mean_pressure')
            p_avg.vector()[:] = p_n.vector()[:]

        # --> Define variational problem for Poisson equation.
        print( 'Assembling Pressure Poisson equation...' )
        Lp = assemble(dot(nabla_grad(p), nabla_grad(q))*dx)
        [bc.apply(Lp) for bc in bcp]
        fp = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

        # --> Define variational problem for velocity update.
        print( 'Assembling velocity update problem...' )
        A = assemble(inner(u, v)*dx)
        b = inner(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

        # --> Set-up the Pressure Poisson solver.
        print( 'Setting-up Pressure Poisson solver...' )
        pressure_solver = LUSolver(Lp, method='umfpack')

        # --> Set-up the correction step solver.
        print( 'Setting-up correction step solver...' )
        update_solver = LUSolver(A, method='umfpack')

        #---------------------------------
        #-----     Time-stepping     -----
        #---------------------------------

        # --> Arrays to store drag and lift.
        drag_time_series, lift_time_series = [], []
        sensors_time_series = []
        if control is not None:
            omega_time_series = []

        t = 0
        output_counter = 1
        for istep in range(1, nsteps+1):

            print('-----------------------------------')
            print('-----     ITERATION %05i     -----' % istep)
            print('-----------------------------------\n')
            # --> Update current time
            t += dt
            print('     --> Time : %0.2f' %t)

            # --> Step 1: Tentative velocity step
            print('     --> Solving the Helmholtz problem.')
            if istep == 1:
                # --> Define the Helholtz operator for 1st order Euler.
                H, f = self._helmholtz_operator_euler(u, v, u_n, p_n, bcu)
                # --> Define the LU solver.
                helhmholtz_solver = LUSolver(H, method='umfpack')
            elif istep == 2:
                # --> Define the Helmholtz operator for BDF2/AB2
                H, f = self._helmholtz_operator_bdf2_ab2(u, v, u_n, u_nn, p_n, bcu)
                # -- Define the LU solver.
                helhmholtz_solver = LUSolver(H, method='umfpack')

            b1 = assemble(f)
            [bc.apply(b1) for bc in bcu]
            helhmholtz_solver.solve(u_.vector(), b1)

            # Step 2: Pressure correction step
            print('     --> Solving the pressure Poisson equation.')
            b2 = assemble(fp)
            [bc.apply(b2) for bc in bcp]
            pressure_solver.solve(p_.vector(), b2)

            # Step 3: Velocity correction step
            print('     --> Divergence-free correction step.\n')
            b3 = assemble(b)
            update_solver.solve(u_.vector(), b3)

            # Update previous solution
            u_nn.assign(u_n)
            u_n.assign(u_)
            p_n.assign(p_)

            # --> Compute drag and lift.
            drag, lift = self.compute_pressure_force(self.mesh, p_n)
            drag_time_series.append(drag)
            lift_time_series.append(lift)

            # --> Sensors measurements.
            sensors = self.sensors_measurements(self.mesh, u_n)
            sensors_time_series.append(sensors)

            # --> Apply control.
            if control is not None:
                omega = control(sensors, t)
                omega_time_series.append(omega)
                front_control.omega = omega[0]
                bottom_control.omega = omega[1]
                top_control.omega = omega[2]

            print('     --> Total lift coefficient : %0.4f' % lift)
            print('     --> Total drag coefficient : %0.4f' % drag)
            print('     --> Angular velocities     : %0.4f, %0.4f, %0.4f \n\n'
                        %(front_control.omega,
                        bottom_control.omega,
                        top_control.omega)
                        )

            # --> Plot solution.
            if np.mod(istep, iostep) == 0:

                # --> Save file.
                filename = figname + '_' + '%05i' %output_counter + '.h5'
                self.outpost_solution(self.mesh, u_, p_, filename)

                if figname is not None:
                    # --> Compute the vorticity field.
                    vorticity = project(curl(u_))

                    # --> Plot the solution.
                    savename = figname + '_' + '%05i' %output_counter
                    self.visualization(self.mesh, vorticity, savename=savename)
                    output_counter += 1
                    plt.close()


            # --> Time-average.
            if avg is True:
                u_avg.vector()[:] += u_n.vector()[:]
                p_avg.vector()[:] += p_n.vector()[:]

        # --> Outpost the mean flow field.
        if avg is True:
            u_avg.vector()[:] /= nsteps
            p_avg.vector()[:] /= nsteps

            # --> Save file.
            filename = 'avg_pinball.h5'
            self.outpost_solution(self.mesh, u_avg, p_avg, filename)

        output = {'lift': lift_time_series,
                  'drag': drag_time_series,
                  'sensors': sensors_time_series,
                  'omega': omega_time_series,
                  }

        output = pd.DataFrame.from_dict(output)

        return output

    def _helmholtz_operator_euler(self, u, v, u_n, p_n, bcu):
        nu = Constant(1/self.re)
        k = Constant(self.dt)

        U = 0.5 * (u_n + u)
        n = FacetNormal(self.mesh)

        problem = inner((u - u_n) / k, v)*dx \
               + inner(dot(u_n, nabla_grad(u_n)), v)*dx \
               + inner(sigma(U, p_n, nu), epsilon(v))*dx \
               + inner(p_n*n, v)*ds - inner(nu*nabla_grad(U)*n, v)*ds

        H = assemble(lhs(problem))
        [bc.apply(H) for bc in bcu]
        f = rhs(problem)

        return H, f

    def _helmholtz_operator_bdf2_ab2(self, u, v, u_n, u_nn, p_n, bcu):
        nu = Constant(1/self.re)
        k = Constant(self.dt)
        n = FacetNormal(self.mesh)

        problem = inner((3*u - 4*u_n + u_nn) / (2*k), v)*dx \
            + inner(dot(1.5*u_n - 0.5*u_nn, nabla_grad(1.5*u_n - 0.5*u_nn)), v)*dx \
            + inner(sigma(u, p_n, nu), epsilon(v))*dx \
            + inner(p_n*n, v)*ds - inner(nu*nabla_grad(u)*n, v)*ds

        H = assemble(lhs(problem))
        [bc.apply(H) for bc in bcu]
        f = rhs(problem)

        return H, f








    def newton(self, savename='steady_state'):

        """
        This function implements a finite-element Newton solver in FENiCS in
        order to compute the stationary equilibrium solutions of the Navier-Stokes
        equations (Pinball configuration). It relies on the basic umfpack
        Newton solver shipped by default with any standard installation of
        FENiCS.
        """

        #--------------------------------
        #-----     Mesh-related     -----
        #--------------------------------

        # --> Define the boundaries of the domain..
        inflow = 'near(x[0], -6)'
        outflow = 'near(x[0], 20)'
        walls = 'near(x[1], -6) || near(x[1], 6)'
        cylinders = 'on_boundary && x[0]>-3. && x[0]<3. && x[1]>-3. && x[1]<3.'

        # --> Inflow velocity profile.
        inflow_profile = ('1', '0')

        #--------------------------------------
        #-----     FE Function spaces     -----
        #--------------------------------------

        # --> Set the Taylor-Hood finite elements.
        P2 = VectorElement('Lagrange', self.mesh.ufl_cell(), 2)
        P1 = FiniteElement('Lagrange', self.mesh.ufl_cell(), 1)

        # --> Define the corresponding function space.
        W = FunctionSpace(self.mesh, MixedElement(P2, P1))

        # --> Define the test functions.
        vq = TestFunction(W)
        (v, q) = split(vq)

        # --> Define the trian function.
        dup = TrialFunction(W)

        #---------------------------------------
        #-----     Variational problem     -----
        #---------------------------------------
        # --> Parameter.
        nu = Constant(1/self.re)

        # --> Define the function.
        up = Function(W)
        (u, p) = split(up)

        # --> Variational problem.
        F = nu*inner(grad(u), grad(v))*dx + inner(grad(u)*u, v)*dx \
            - div(v)*p*dx + div(u)*q*dx

        # --> Compute the Jacobian.
        dF = derivative(F, up, dup)

        # --> Sets the boundary conditions.
        bc_inflow = DirichletBC(W.sub(0), inflow_profile, inflow)
        bc_walls = DirichletBC(W.sub(0), inflow_profile, walls)
        bc_cylinders = DirichletBC(W.sub(0), Constant((0, 0)), cylinders)
        bc_outflow = DirichletBC(W.sub(1), Constant(0), outflow)
        bc = [bc_inflow, bc_walls, bc_cylinders, bc_outflow]

        # --> Define the final variational problem.
        problem = NonlinearVariationalProblem(F, up, bc, dF)

        #------------------------------------
        #-----     Nonlinear solver     -----
        #------------------------------------

        # --> Create the nonlinear solver.
        solver = NonlinearVariationalSolver(problem)

        # --> Solve the problem.
        solver.solve()

        # --> Save the solution to HDF5 file.
        if savename is not None:
            # --> Extract solution.
            (u, p) = up.split()

            filename = savename + '.h5'
            self.outpost_solution(self.mesh, u, p, filename)

            # --> Compute the vorticity field.
            vorticity = project(curl(u))

            # --> Plot the solution.
            self.visualization(self.mesh, vorticity, savename=savename)

            return

        else:
            return u, p





    @staticmethod
    def visualization(mesh, vorticity, psi=None, savename='figure', var_range=[-4, 4]):

        # --> Import various utilities from matplotlib.
        from matplotlib.tri import Triangulation
        from matplotlib import ticker
        from matplotlib.patches import Circle
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        # --> Get mesh information and build-up the Delaunay triangulation.
        x, y = mesh.coordinates()[:, 0], mesh.coordinates()[:, 1]
        tri = Triangulation(x, y)

        # --> Mask the unwanted triangles.
        radius = 0.5
        xmid = tri.x[tri.triangles].mean(axis=1)
        ymid = tri.y[tri.triangles].mean(axis=1)

        xc1, yc1 = 0, 0.75
        xc2, yc2 = 0, -0.75
        xc3, yc3 = -1.5*np.sqrt(0.75), 0
        mask = np.where((xmid-xc1)**2 + (ymid-yc1)**2 < radius**2, 1, 0)
        mask += np.where((xmid-xc2)**2 + (ymid-yc2)**2 < radius**2, 1, 0)
        mask += np.where((xmid-xc3)**2 + (ymid-yc3)**2 < radius**2, 1, 0)
        tri.set_mask(mask)

        # --> Get the vorticity field at the mesh nodes.
        vorticity = vorticity.compute_vertex_values(mesh)

        #-----------------------------------
        #-----     PLOT THE FIGURE     -----
        #-----------------------------------

        # --> Create matplotlib figure.
        fig = plt.figure()
        ax = fig.gca()

        # --> Plot the vorticity field.
        h = ax.tripcolor(tri, vorticity, shading='gouraud', cmap=plt.cm.RdBu)
        h.set_clim(var_range[0], var_range[1])

        # --> Set the axis.
        ax.set_aspect('equal')
        ax.set_xlim(-5, 20)
        ax.set_ylim(-4, 4)

        # --> Setup the colorbar.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.1)
        cb = plt.colorbar(h, cax=cax)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()

        # --> Highlight the cylinders.
        xc, yc, r = 0., 0.75, 0.5
        circle = Circle((xc, yc), r, facecolor='None', edgecolor='k')
        ax.add_patch(circle)

        xc, yc, r = 0., -0.75, 0.5
        circle = Circle((xc, yc), r, facecolor='None', edgecolor='k')
        ax.add_patch(circle)

        xc, yc, r = -1.5*np.sqrt(0.75), 0., 0.5
        circle = Circle((xc, yc), r, facecolor='None', edgecolor='k')
        ax.add_patch(circle)

        # --> Titles, labels, ...
        ax.set_xlabel(r'$x$')
        ax.locator_params(axis='y', nbins=5)
        ax.set_ylabel(r'$y$', rotation=0)

        # --> Save the figure.
        fig.savefig(savename+'.png', bbox_inches='tight', dpi=300)

        return





    @staticmethod
    def compute_streamfunction(mesh, u):
        pass





    @staticmethod
    def compute_pressure_force(mesh, p):

        # --> Define the integration domain.
        class cylinder(SubDomain):
            def inside(self, x, on_boundary):
                return x[0]>-3. and x[0]<3. and x[1]>-3. and x[1]<3.

        # --> ?
        markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 1)
        cylinder().mark(markers, 1)

        # --> Define functionals for drag and lift.
        ds = Measure('ds', subdomain_data=markers)
        n = FacetNormal(mesh)
        D = -p*n[0]*ds(1)
        L = p*n[1]*ds(1)

        # --> Assemble functionals over sub-domain.
        drag = assemble(D)
        lift = assemble(L)

        return drag, lift





    @staticmethod
    def outpost_solution(mesh, u, p, filename):

        # --> Open the HDF5 file.
        #f = HDF5File(mpi_comm_world(), filename, 'w')
        f = HDF5File(MPI.comm_world, filename, 'w')
        # --> Outpost mesh.
        f.write(mesh, 'mesh')
        # --> Outpost velocity field.
        f.write(u, 'velocity')
        # --> Outpost pressure field (if needed)
        if p is not None:
            f.write(p, 'pressure')

        # --> Close the file.
        f.close()

        return





    @staticmethod
    def load_solution(filename, up=None):

        # --> Open HDF5 file.
        f = HDF5File(MPI.comm_world, filename, 'r')

        if up is None:
            mesh = Mesh()
            f.read(mesh, 'mesh', False)
            # --> Load velocity field.
            V = VectorFunctionSpace(mesh, 'CG', 2)
            u = Function(V)
            f.read(u, 'velocity')

            # --> Load pressure field.
            Q = FunctionSpace(mesh, 'CG', 1)
            p = Function(Q)
            f.read(p, 'pressure')

            # --> Close HDF5 file.
            f.close()

            return mesh, u, p
        else:
            f.read(up[0], 'velocity')
            f.read(up[1], 'pressure')
            f.close()

            return up[0], up[1]




    @staticmethod
    def sensors_measurements(mesh, u):

        # --> Create the rack of sensors.
        xs = np.ones((5,))*5.5
        ys = np.zeros_like(xs)
        for i in range(ys.size):
            ys[i] = 0.75*(2-i)

        # # --> Evaluate velocity measurements.
        s = np.zeros((xs.size, 2))
        for i in range(xs.size):
            s[i] = u(xs[i], ys[i])

        return s





if __name__ == "__main__":

    # --> Initialize the pinball object.
    pinball = Pinball(re=120)

    # --> Newton solve.
    pinball.newton()

    # --> Define control.
    def control(sensors, t):

        # --> Initialize omega.
        omega = np.zeros(3)

        # --> Control.
        if t < 1:
            omega[:] = 1.0

        return omega

    # --> Run the (unforced) simulation.
    output = pinball.dns(T=100, iostep=100, initial_condition='steady_state.h5', control=control)

    plt.show()
