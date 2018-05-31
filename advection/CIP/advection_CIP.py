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

    def Next_u_v(self, u, v):
        u_up = np.roll(u, -1)
        v_up = np.roll(v, -1)
        D = self.dx
        a = (v + v_up)/(D**2) + 2*(u - u_up)/(D**3)
        b = 3*(u_up - u)/(D**2) - (2*v + v_up)/D
        xi = self.dt
        u_next = a*xi**3 + b*xi**2 + v*xi + u
        v_next = 3*a*xi**2 + 2*b*xi + v
        return u_next, v_next

    def derivative_x(self, u):
        u_p = np.roll(u, -1)
        u_m = np.roll(u, +1)
        return (u_p - u_m)/(2*self.dx)
    
adv = advection()
u = np.sin(4*np.pi*adv.x / adv.Lx)
v = adv.derivative_x(u)

f = open("adv.dat","w")
for ix in range(adv.Nx):
    f.write(str(adv.x[ix]) + " " + str(u[ix]) + "\n")

f.write("\n\n")

for it in range(adv.Nt + 1):
    u, v = adv.Next_u_v(u,v)
    if it % adv.out == 0:
        for ix in range(adv.Nx):
            f.write(str(adv.x[ix]) + " " + str(u[ix]) + "\n")

        f.write("\n\n")

f.close()

