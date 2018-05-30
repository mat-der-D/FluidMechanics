import numpy as np

class advection:
    def __init__(self):
        self.Lx = 50
        self.Nx = 500
        self.dt = 0.01
        self.Nt = 10000
        self.dx = self.Lx / self.Nx
        self.x = np.linspace(0, self.Lx, self.Nx, endpoint = False)
        self.out = 100

    def Next_u(self, u):
        u_p = np.roll(u, -1)
        u_m = np.roll(u, +1)
        du = (u_p - u_m)/(2 * self.dx) * self.dt
        return u + du

adv = advection()
u = np.sin(4*np.pi*adv.x / adv.Lx)

f = open("adv.dat","w")
for ix in range(adv.Nx):
    f.write(str(adv.x[ix]) + " " + str(u[ix]) + "\n")

f.write("\n\n")

for it in range(adv.Nt + 1):
    u = adv.Next_u(u)
    if it % adv.out == 0:
        for ix in range(adv.Nx):
            f.write(str(adv.x[ix]) + " " + str(u[ix]) + "\n")

        f.write("\n\n")

f.close()

