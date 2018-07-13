import numpy as np

class diffusion2D_ADI:

    def __init__(self):
        self.Lx, self.Lz = 50.0, 50.0
        self.Nx, self.Nz = 100, 100
        self.dx, self.dz = self.Lx/self.Nx, self.Lz/self.Nz
        self.D = 1.0
        self.dt = 0.001
        self.Nt = 10000
        self.out = 100

        self.x = np.linspace(0., self.Lx, self.Nx, endpoint = False)
        self.z = np.linspace(0., self.Lz, self.Nz, endpoint = False)
        self.zz, self.xx = np.meshgrid(self.z, self.x)

        E_x, E_z = np.eye(self.Nx), np.eye(self.Nz)
        lmb_x = self.D*self.dt / (2*self.dx**2)
        lmb_z = self.D*self.dt / (2*self.dz**2)

        A_x = self.init_mat(lmb_x, self.Nx)
        A_z = self.init_mat(lmb_z, self.Nz)
        B_x = self.init_mat(-lmb_x, self.Nx)
        B_z = self.init_mat(-lmb_z, self.Nz)

        self.A_1 = np.dot(B_x, np.linalg.inv(A_x))
        self.A_2 = np.dot(B_z, np.linalg.inv(A_z))

    def init_mat(self, x, n_size):
        E = np.eye(n_size)
        A = -x*np.roll(E, +1, 0) - x*np.roll(E, -1, 0) \
            + (1 + 2*x)*E
        return A

    def init_val(self):
        u = np.cos(4*np.pi*self.xx/self.Lx)*np.sin(4*np.pi*self.zz/self.Lz)
        return u

    def Next_u(self, u):
        return np.dot(self.A_1, np.dot(u, self.A_2))

    def output_u(self, u, f):
        for i in range(self.Nx):
            for k in range(self.Nz):
                f.write(str(self.x[i]) + " " +\
                        str(self.z[k]) + " " +\
                        str(u[i,k]) + "\n")
            f.write("\n")
        f.write("\n")

df = diffusion2D_ADI()
u = df.init_val()

f = open("diffusion.dat", "w")

df.output_u(u, f)

for it in range(1, df.Nt + 1):
    t = it * df.dt
    u = df.Next_u(u)
    if it % df.out == 0:
        df.output_u(u, f)
        print("it=" + str(it))

f.close()
