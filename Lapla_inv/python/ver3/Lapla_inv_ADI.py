import numpy as np

class Lapla_inv_ADI:

    def __init__(self):
        self.Lx = 10.
        self.Nx = 100
        self.Lz = 10.
        self.Nz = 100
        self.dx, self.dz = self.Lx/self.Nx, self.Lz/self.Nz
        self.x = np.linspace(0., self.Lx, self.Nx, endpoint=False)
        self.z = np.linspace(0., self.Lz, self.Nz+1, endpoint=True)
        self.zz, self.xx = np.meshgrid(self.z, self.x)
        
        self.dt = 1./( 4*(self.dx**(-2) + self.dz**(-2)) )
        
        E_x, E_z = np.eye(self.Nx), np.eye(self.Nz)
        lmb_x = self.dt / (2*self.dx**2)
        lmb_z = self.dt / (2*self.dz**2)

        A_x = self.init_mat(lmb_x, self.Nx, False)
        A_z = self.init_mat(lmb_z, self.Nz+1, True)
        B_x = self.init_mat(-lmb_x, self.Nx, False)
        B_z = self.init_mat(-lmb_z, self.Nz+1, True)

        self.A_1 = np.dot(B_x, np.linalg.inv(A_x))
        self.A_2 = np.dot(B_z, np.linalg.inv(A_z))
        self.A_l = B_x
        self.A_r = np.linalg.inv(A_z)

    def init_mat(self, x, n_size, bdry):
        E = np.eye(n_size)
        A = -x*np.roll(E, +1, 0) - x*np.roll(E, -1, 0) \
            + (1 + 2*x)*E
        if bdry:
            A[n_size -1, 0] = 0.
            A[0, n_size -1] = 0.
        
        return A

    def Next_u(self, u, omega):
        new = np.dot(self.A_1, np.dot(u, self.A_2)) \
              + (self.dt/2)* \
              (np.dot(self.A_l, np.dot(omega, self.A_r)) \
               + omega)
        new[:,0] = 0.
        new[:,self.Nz] = 0.
        return new

    def u_Lapla_u(self, u):
        u_xx = self.derivative_xx(u)
        u_zz = self.derivative_zz(u)
        return u_xx + u_zz

    def derivative_xx(self, f):
        ddf = (np.roll(f, -1, 0) + np.roll(f, +1, 0) - 2*f) \
              / ( self.dx**2 )
        return ddf

    def derivative_zz(self, f):
        ddf = (np.roll(f, -1, 1) + np.roll(f, +1, 1) - 2*f) \
              / ( self.dz**2 )
        ddf[:,0] = ddf[:,1]
        ddf[:,self.Nz] = ddf[:,self.Nz-1]
        return ddf

    def norm(self, u):
        return np.sqrt( (u**2).sum()/(self.Nx*(self.Nz + 1)) )
    
li = Lapla_inv_ADI()
    
psi = np.sin(2*np.pi*li.zz/li.Lz)
omega = - li.u_Lapla_u(psi)
psi *= 0.

f = open("error.d", "w")

for i in range(10000):
    psi = li.Next_u(psi, omega)
    err = li.norm(omega + li.u_Lapla_u(psi))
    f.write(str(i) + " " + str(err) + "\n")

f.close()
