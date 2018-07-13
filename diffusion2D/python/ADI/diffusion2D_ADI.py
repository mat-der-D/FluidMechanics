import numpy as np

class diffusion2D_ADI:

    def __init__(self):
        self.Lx, self.Lz = 15.0, 10.0
        self.Nx, self.Nz = 100, 100
        self.dx, self.dz = self.Lx/self.Nx, self.Lz/self.Nz
        self.D = 0.1
        self.dt = 0.005
        self.Nt = 20000
        self.out = 200
        self.alpha = self.D*self.dt/(2*self.dx**2)
        self.beta = self.D*self.dt/(2*self.dz**2)
        x = np.linspace(0.0, self.Lx, self.Nx, endpoint=False)
        z = np.linspace(0.0, self.Lz, self.Nz, endpoint=False)
        self.zz, self.xx = np.meshgrid(z, x)
        
    def F(self, c, u):
        # 普通に使う場合は, x に関する時間発展.
        # z に関する時間発展に使う場合は, あらかじめ u を転置
        u_plus  = np.roll(u, +1, 1)
        u_minus = np.roll(u, -1, 1)
        return c*u_plus + (1 - 2*c)*u + c*u_minus

    def init_mat(self, c, N):
        E = np.eye(N)
        A = (1 + 2*c)*E - c*np.roll(E, 1, 1) - c*np.roll(E, -1, 1)
        Ainv = np.linalg.inv(A)
        return Ainv

    def Next_u(self, u, A_x, A_z):
        u_mid = np.dot(A_x, self.F(self.beta, u))
        u_next_T = np.dot(A_z, self.F(self.alpha, u_mid.T))
        return u_next_T.T

    def init_u(self):
        u = np.sin(2*np.pi*self.xx/self.Lx) \
            + np.cos(4*np.pi*self.zz/self.Lz)
        return u

    def exact_u(self, t):
        u = np.sin(2*np.pi*self.xx/self.Lx) \
            * np.exp(- 4*np.pi**2*self.D/(self.Lx**2)*t)\
            + np.cos(4*np.pi*self.zz/self.Lz) \
            * np.exp(- 16*np.pi**2*self.D/(self.Lz**2)*t)
        return u
    
    def output_u(self, u, f):
        for ix in range(self.Nx):
            for iz in range(self.Nz):
                f.write(str(self.xx[ix, iz]) + " " + \
                        str(self.zz[ix, iz]) + " " + \
                        str(u[ix, iz]) + "\n")
            f.write("\n")
        f.write("\n")

    def norm(self, v):
        v_flat = v.flatten()
        return np.sqrt(np.dot(v_flat, v_flat)/(self.Nx*self.Nz))

df = diffusion2D_ADI()
Ax = df.init_mat(df.alpha, df.Nx)
Az = df.init_mat(df.beta, df.Nz)
u_init = df.init_u()
u = u_init

f = open("diffusion.dat", "w")
ff = open("exact.dat", "w")
g = open("error.dat", "w")
h = open("relative_error.dat", "w")


exact = df.exact_u(0.0)

df.output_u(u, f)
df.output_u(exact, ff)

for it in range(1, df.Nt + 1):
    t = it * df.dt
    u = df.Next_u(u, Ax, Az)
    if it % df.out == 0:
        df.output_u(u, f)
        # exact = u_init*np.exp(- 4*np.pi**2*df.D/(df.Lx**2)*t)
        exact = df.exact_u(t)
        # error = (np.abs(exact - u)).max()
        error = df.norm(exact - u)
        relative_error = error*np.exp(+ 4*np.pi**2*df.D/(df.Lx**2)*t)
        df.output_u(exact, ff)
        g.write(str(t) + " " + str(error) + "\n")
        h.write(str(t) + " " + str(relative_error) + "\n")
    if it % (10*df.out) == 0:
        print(it)
        
f.close()
g.close()
h.close()
