import numpy as np

class advection:
    def __init__(self):
        self.Lx = 50
        self.Nx = 500
        self.dt = 0.01
        self.Nt = 50000
        self.dx = self.Lx / self.Nx
        self.x = np.linspace(0, self.Lx, self.Nx, endpoint = False)
        self.out = 500
        E = np.eye(self.Nx)
        lmb = self.dt/(4*self.dx)
        A_coef = lmb * np.roll(E, +1, axis = 1) \
                 - E \
                 - lmb *np.roll(E, -1, axis = 1)
        B_coef = - lmb * np.roll(E, +1, axis = 1) \
                 - E \
                 + lmb * np.roll(E, -1, axis = 1)
        self.A_inv_B = np.dot(np.linalg.inv(A_coef), B_coef)

    def Next_u(self, u):
        return np.dot(self.A_inv_B, u)
    
adv = advection()
u = np.sin(4*np.pi*adv.x / adv.Lx)

f = open("adv.dat","w")
for ix in range(adv.Nx):
    f.write(str(adv.x[ix]) + " " + str(u[ix]) + "\n")

f.write("\n\n")

for it in range(1, adv.Nt + 1):
    u = adv.Next_u(u)
    if it % adv.out == 0:
        for ix in range(adv.Nx):
            f.write(str(adv.x[ix]) + " " + str(u[ix]) + "\n")

        f.write("\n\n")

f.close()
