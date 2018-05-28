import numpy as np

class params:

    def __init__(self):
        self.Lx, self.Ly = 10, 10
        self.Nx, self.Ny = 20, 20
        self.dx, self.dy = 1.0/self.Nx, 1.0/self.Ny

    def mesh_init(self):
        Lx, Ly = self.Lx, self.Ly
        Nx, Ny = self.Nx, self.Ny
        self.x = np.linspace(-Lx, Lx, 2*Lx*Nx + 1)
        self.y = np.linspace(-Ly, Ly, 2*Ly*Ny + 1)

m = params()
m.mesh_init()
print(m.x)
