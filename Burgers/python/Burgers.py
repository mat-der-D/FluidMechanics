import numpy as np

class Burgers:

    def __init__(self):
        self.Lx = 10.0
        self.Nx = 100
        self.dx = self.Lx / self.Nx
        self.xx = np.linspace(0.0, self.Lx, self.Nx, endpoint = False)
        self.dt = 0.0001
        self.t_end = 10000
        self.nu = 1.0
        self.t_out = 100
        self.A_one = self.mat_CN(self.dx, self.Nx, self.dt, self.nu)

    def mat_CN(self, dx, Nx, dt, nu):
        lmb = nu * dt / (dx ** 2)
        E = np.eye(Nx)
        A = - lmb * np.roll(E, -1, 1) - lmb * np.roll(E, +1, 1) \
            + (1 + 2*lmb) * E
        return np.linalg.inv(A)
        
    def f(self, u, t):
        '''非線形項'''
        du = (np.roll(u, -1) - np.roll(u, +1)) / (self.dx ** 2)
        return - u * du

    def u_init(self):
        return np.sin(2 * np.pi * self.xx / self.Lx) + 0.5

    def Next_u(self, u, t, A_one):
        u_mid = u + self.dt * self.f(u, t)
        u_new = np.dot(A_one, u_mid)
        return u_new

    def u_write(self, f_out, u):
        for n in range(self.Nx):
            f_out.write(str(self.xx[n]) + " " + str(u[n]) + "\n")
        f_out.write("\n\n")
    
    
b = Burgers()
A_one = b.A_one
u = b.u_init()

fb = open('Burgers.dat', 'w')
b.u_write(fb, u)

for it in range(b.t_end +1):
    t = it * b.dt
    u = b.Next_u(u, t, A_one)
    if it % b.t_out == 0:
        b.u_write(fb, u)

fb.close()
