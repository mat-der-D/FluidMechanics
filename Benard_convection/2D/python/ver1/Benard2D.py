import sys
import numpy as np
from scipy.ndimage.interpolation import shift

class Benard_2D:

    def __init__(self):
        self.Lx = 2*np.sqrt(2) #= 2*pi/(pi/sqrt(2))
        self.Nx = 30
        self.Lz = 1.
        self.Nz = 30
        self.dx, self.dz = self.Lx/self.Nx, self.Lz/self.Nz        
        self.x = np.linspace(0., self.Lx, self.Nx, endpoint=False)
        self.z = np.linspace(0., self.Lz, self.Nz+1, endpoint=True)
        self.zz, self.xx = np.meshgrid(self.z, self.x)
        self.dt = 0.001
        self.Nt = 5000
        self.tout = 100
        
        self.T_U = 0.
        self.T_L = 1.
        self.Ra = 27*np.pi**4/4 * 3
        self.Pr = 0.01

        # For diffusion eq.
        lmb_x = self.dt / (2*self.dx**2)
        lmb_z = self.dt / (2*self.dz**2)

        # used for T
        self.A_1, self.A_2, self.A_3 = \
            self.init_diffusion(lmb_x, lmb_z, self.T_L, self.T_U)

        # used for omega
        self.B_1, self.B_2, self.B_3 = \
            self.init_diffusion(self.Pr*lmb_x, self.Pr*lmb_z, \
                                0., 0.)

        # For Lapla_inv
        E_x, E_z = np.eye(self.Nx), np.eye(self.Nz+1)
        
        self.mu_x = 1./( 2*(1 + (self.dx/self.dz)**2) )
        self.mu_z = 1./( 2*(1 + (self.dz/self.dx)**2) )
        self.mu_omega = 1./( 2*(self.dx**(-2) + self.dz**(-2)) )

        self.A_Li_x = self.mu_x*(np.roll(E_x, +1, 0) \
                                 + np.roll(E_x, -1, 0))
        self.A_Li_z = self.mu_z*(np.roll(E_z, +1, 0) \
                                 + np.roll(E_z, -1, 0))
        self.A_Li_z[:,0] = 0.
        self.A_Li_z[:,self.Nz] = 0.
        
    def init_diffusion(self, lmb_x, lmb_z, u_L, u_U):
        A_x = self.init_mat_p(lmb_x, self.Nx)
        A_z = self.init_mat_f(lmb_z, self.Nz-1, u_L, u_U)
        B_x = self.init_mat_p(-lmb_x, self.Nx)
        B_z = self.init_mat_f(-lmb_z, self.Nz-1, u_L, u_U)
        U_BC = np.zeros((self.Nx, self.Nz-1))
        U_BC[:,0] = u_L
        U_BC[:,self.Nz-2] = u_U
        
        A_1 = np.dot(B_x, np.linalg.inv(A_x))
        A_2 = np.dot(B_z, np.linalg.inv(A_z))
        A_3 = lmb_z \
              * np.dot(A_1 + np.eye(self.Nx), \
                       np.dot(U_BC, np.linalg.inv(A_z)))
        return A_1, A_2, A_3

        
    def init_mat_p(self, x, n_size):
        E = np.eye(n_size)
        A = -x*np.roll(E, +1, 0) -x*np.roll(E, -1, 0) \
            + (1 + 2*x)*E
        return A

    def init_mat_f(self, x, n_size, u_L, u_U):
        E = np.eye(n_size)
        N = np.roll(E, +1, 0) + np.roll(E, -1, 0)
        N[n_size-1, 0], N[0, n_size-1] = 0., 0.
        A = -x*N + (1 + 2*x)*E
        return A

    def Next_psi_omega_T(self, psi, omega, T):
        psi_mid, omega_mid, T_mid = \
                self.nonlin_Euler(psi, omega, T)
        omega_new = self.Next_diff_u(omega_mid)
        T_new = self.Next_diff_T(T_mid)
        psi_new = self.psi_omega(omega_new, psi_mid)
        return psi_new, omega_new, T_new
    
    def Next_diff_u(self, u):
        u_new = np.zeros((self.Nx, self.Nz+1))
        u_new[:,1:self.Nz] = \
            np.dot(self.B_1, np.dot(u[:,1:self.Nz], self.B_2)) \
            + self.B_3
        return u_new

    def Next_diff_T(self, T):
        T_new = np.zeros((self.Nx, self.Nz+1))
        T_new[:,0] = self.T_L
        T_new[:,self.Nz] = self.T_U
        T_new[:,1:self.Nz] = \
                np.dot(self.A_1, np.dot(T[:,1:self.Nz], self.A_2)) \
                + self.A_3
        return T_new
    
    def derivative_x(self, f):
        return (np.roll(f, -1, 0) - np.roll(f, +1, 0))/(2*self.dx)

    def derivative_z(self, f):
        df = (np.roll(f, -1, 1) - np.roll(f, +1, 1))/(2*self.dz)
        df[:,0] = (f[:,1] - f[:,0])/(self.dz)
        df[:,self.Nz] = (f[:,self.Nz] - f[:,self.Nz-1])/(self.dz)
        return df

    def derivative_xx(self, f):
        ddf = (np.roll(f, -1, 0) + np.roll(f, +1, 0) - 2*f) \
              / ( self.dx**2 )
        return ddf

    def derivative_zz(self, f):
        '''境界値は 0'''
        ddf = (np.roll(f, -1, 1) + np.roll(f, +1, 1) - 2*f) \
              / ( self.dz**2 )
        ddf[:,0] = 0.
        ddf[:,self.Nz] = 0.
        return ddf
    
    def Jacobi(self, a, b):
        a_x = self.derivative_x(a)
        a_z = self.derivative_z(a)
        b_x = self.derivative_x(b)
        b_z = self.derivative_z(b)
        return a_z*b_x - a_x*b_z

    def nonlin_omega(self, psi, omega, T):
        T_x = self.derivative_x(T)
        out = - self.Ra*self.Pr*T_x - self.Jacobi(psi, omega)
        return out

    def nonlin_Euler(self, psi, omega, T):
        dlt_omega = self.dt*self.nonlin_omega(psi, omega, T)
        dlt_T = - self.dt*self.Jacobi(psi, T)
        omega_new = omega + dlt_omega
        T_new = T + dlt_T
        psi_new = self.psi_omega(omega_new, psi)
        return psi_new, omega_new, T_new

    def u_Lapla_u(self, u):
        '''境界値は 0'''
        u_xx = self.derivative_xx(u)
        u_zz = self.derivative_zz(u)
        return u_xx + u_zz

    def norm(self, u):
        return np.sqrt( (u**2).sum()/(self.Nx*(self.Nz + 1)) )
    
    def psi_omega(self, omega, psi_old):
        '''\nabla^2 psi = omega を解く'''
        eps = 10.**(-3) # 誤差評価
        psi = psi_old
        istop = 10000 #エラー終了
        error=1.
        i = 0
        while error>eps:
            for icount in range(10):
                psi = np.dot(self.A_Li_x, psi) \
                      + np.dot(psi, self.A_Li_z) \
                      - self.mu_omega * omega
                psi[:, 0] = 0.
                psi[:, self.Nz] = 0.
            error = self.norm(omega - self.u_Lapla_u(psi)) \
                    / self.norm(omega)
            i += 1
            if i >= istop:
                print("Error occurred!! error=" + str(error))
                g = open("E_Laplapsi.d", "w")
                gg = open("E_psi.d", "w")
                h = open("E_omega.d", "w")

                self.output_u(psi, gg)
                self.output_u(self.u_Lapla_u(psi), g)
                self.output_u(omega, h)
                g.close()
                h.close()
                sys.exit()
        return psi

    def output_u(self, u, ff):
        for i in range(self.Nx):
            for k in range(self.Nz+1):
                ff.write(str(self.x[i]) + " " + \
                         str(self.z[k]) + " " + \
                         str(u[i, k])   + "\n")
            ff.write("\n")
        ff.write("\n")

    def output_T(self, T, ff):
        for i in range(2*self.Nx+1):
            ii = i % self.Nx
            for k in range(self.Nz+1):
                ff.write(str(i * self.dx) + " " + \
                         str(self.z[k])   + " " + \
                         str(T[ii, k])    + "\n")
            ff.write("\n")
        ff.write("\n")

    def output_velocity(self, psi, ff):
        u = self.derivative_z(psi)
        w = - self.derivative_x(psi)
        for i in range(self.Nx):
            for k in range(self.Nz+1):
                ff.write(str(self.x[i]) + " " + \
                         str(self.z[k]) + " " + \
                         str(u[i,k]) + " " + \
                         str(w[i,k]) + "\n")
            # ff.write("\n")
        ff.write("\n\n")
        
    def initialize_psi_omega_T(self):
        psi = np.sin(np.pi*self.zz/self.Lz) \
              * np.exp( np.cos(2*np.pi*self.xx/self.Lx) )
        omega = self.u_Lapla_u(psi)
        T = - self.zz + 1 + np.sin(np.pi*self.zz/self.Lz)
        return psi, omega, T

bc = Benard_2D()

psi, omega, T = bc.initialize_psi_omega_T()

ff = open("vel.d", "w")
gg = open("Temp.d", "w")
hh = open("norm_omega.d", "w")

# bc.output_u(psi, ff)
bc.output_velocity(psi, ff)
bc.output_T(T, gg)
hh.write(str(0.) + " " + str(bc.norm(omega)) + "\n")
# bc.output_u(omega, ff)

for it in range(1, bc.Nt + 1):
    psi, omega, T = bc.Next_psi_omega_T(psi, omega, T)
    if it % bc.tout == 0:
        t = it * bc.dt
        print("it =" + str(it) + "\n")
        bc.output_velocity(psi, ff)
        bc.output_T(T, gg)
        hh.write(str(t) + " " + str(bc.norm(omega)) + "\n")

ff.close()
gg.close()
hh.close()
