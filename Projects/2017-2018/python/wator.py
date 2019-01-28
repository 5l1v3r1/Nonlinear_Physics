import random

class fish(object):
    gperiod = 3 # Periode de gestation.

    def __init__(self):
        self.gclock = 0

    def evolve(self):
        empty = [neighbor for neighbor in self.cell.neighbors if neighbor.content is None]

        if empty:                               # if there's an empty cell nearby:
            if self.gclock >= self.gperiod:     #   if the gestation clock is past due:
                self.cell.place(type(self)())   #     place a new animal at the current position
                self.gclock = 0                 #     reset the gestation clock
            else:
                self.cell.content = None        #   otherwise, empty the current location

            dest = random.choice(empty)         #   choose random empty neighbor and move there
            dest.place(self)

        self.gclock = self.gclock + 1


class shark(fish):
    gperiod = 10 # Periode de gestation.
    speriod  = 3 # Survie sans nourriture.

    def __init__(self):
        self.gclock = 0
        self.sclock = 0

    def evolve(self):
        food = [neighbor for neighbor in self.cell.neighbors if type(neighbor.content) is fish]

        if food:                          # if there is fish nearby:
            prey = random.choice(food)    #   choose a fish at random
            prey.content = None           #   remove it from the grid
            self.sclock = 0               #   reset the starvation clock

        super(shark,self).evolve()        # move and spawn using fish.evolve()

        if self.sclock >= self.speriod:   # if the starvation clock is expired
            self.cell.content = None      #   die and remove ourselves from the grid
        else:
            self.sclock = self.sclock + 1


class cell(object):
    def __init__(self,grid):
        self.content = None

    def place(self,animal):
        self.content = animal
        animal.cell = self


class grid(object):
    def __init__(self,M,N_fish,N_shark):
        """Create an MxM Wator grid and populate it wth N_fish fish and N_shark sharks."""

        self.M = M
        self.cells = [[cell(self) for j in range(M)] for i in range(M)]

        for i in range(M):
            for j in range(M):
                # square-grid neighbors, periodic boundary conditions
                self.cells[i][j].neighbors = (self.cells[i][(j+1) % M],
                                              self.cells[i][(j-1) % M],
                                              self.cells[(i+1) % M][j],
                                              self.cells[(i-1) % M][j])

        self.allfish   = [fish()  for i in range(N_fish)]
        self.allsharks = [shark() for i in range(N_shark)]

        # keep a flattened list of empty cells
        empty = [gridcell for row in self.cells for gridcell in row]

        for animal in self.allfish + self.allsharks:
            # remove a random cell from `empty`, place animal there
            randomcell = empty.pop(random.randrange(len(empty)))
            randomcell.place(animal)

    def evolve(self):
        """Evolve the Wator grid through one time step, return the tuple (N_fish,N_shark)."""

        for animal in self.allfish + self.allsharks:
            animal.evolve()

        # gather all the fish and sharks from the grid
        self.allfish   = [gridcell.content for row in self.cells
                          for gridcell in row if type(gridcell.content) is fish]
        self.allsharks = [gridcell.content for row in self.cells
                          for gridcell in row if type(gridcell.content) is shark]

        return len(self.allfish), len(self.allsharks)

    def spatial_distribution(self):
        x = np.zeros((self.M, self.M))
        for i in xrange(self.M):
            for j in xrange(self.M):
                if type(self.cells[i][j].content) == fish:
                    x[i, j] = 1
                elif type(self.cells[i][j].content) == shark:
                    x[i, j] = -1
        return x


if __name__ == "__main__":

    # --> Import standard python library.
    import numpy as np
    import matplotlib.pyplot as plt

    # --> Define the grid size and the number of generations.
    M, steps = 128, 1000
    mygrid = grid(M=M, N_fish=int(0.2*M**2), N_shark=int(0.05*M**2))

    populations = np.zeros((steps,2))
    counter = 0
    for step in xrange(steps):
        # --> Evolve the populations.
        populations[step,:] = mygrid.evolve()
        # --> Plot the spatial distribution.
        if np.mod(step, 10) == 0:
            x = mygrid.spatial_distribution()
            fig = plt.figure(figsize=(3, 3))
            ax = fig.gca()
            ax.imshow(x, cmap=plt.cm.bwr)
            plt.savefig('wator_world_%i.png' %counter, bbox_inches='tight', dpi=300)
            plt.close()
            counter += 1

    # --> Plot time-series of the shark and fish.
    fig = plt.figure(figsize=(6, 2))
    ax = fig.gca()

    ax.plot(populations[:,0], label='fish')
    ax.plot(populations[:,1], label='shark')

    ax.set_ylabel('Population')
    ax.set_xlabel('Generation')

    ax.legend(loc=0)

    plt.savefig('populations_time_series.pdf', bbox_inches='tight', dpi=300)

    # --> Save time series.
    np.savetxt('populations.dat', populations)

    plt.show()
