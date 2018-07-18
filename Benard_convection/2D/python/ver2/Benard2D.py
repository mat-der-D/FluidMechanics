import sys
import numpy as np

class Settings:

    def __init__(self):
        self.Lx = 2*np.sqrt(2)
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

        E_x = np.eye(self.Nx)
        E_z = np.eye(self.Nz+1)
        
        self.A_der_x \
            = (np.roll(E_x, +1, 0) - np.roll(E_x, -1, 0)) \
              / (2*self.dx)
        self.A_der_z = np.zeros((self.Nz+1, self.Nz+1))
        self.A_der_z[:, 1:self.Nz] \
            = ( (np.roll(E_z, +1, 1) - np.roll(E_z, -1, 1)) \
              / (2*self.dz) )[:, 1:self.Nz]
        self.A_der_z[0, 0:3] \
            = np.array([-3.,4.,-1.])/(2*self.dz)
        self.A_der_z[self.Nz, self.Nz-2:self.Nz+1] \
            = np.array([1.,-4.,3.])/(2*self.dz)
                
        self.A_der_xx \
            = (np.roll(E_x, +1, 0) - 2*E_x \
               + np.roll(E_x, -1, 0)) / (self.dx**2)
        self.A_der_zz = np.zeros((self.Nz+1, self.Nz+1))
        self.A_der_zz[:,1:self.Nz] \
            = ( (np.roll(E_z, +1, 0) - 2*E_z \
                 + np.roll(E_z, -1, 0)) \
                / (self.dz**2) )[:, 1:self.Nz]
        self.A_der_zz[0, 0:4] \
            = np.array([2., -5., 4., -1.]) / (self.dz**2)
        self.A_der_zz[self.Nz, self.Nz-3:self.Nz+1] \
            = np.array([-1., 4., -5., 2.]) / (self.dz**2)

    def derivative_x(self, f):
        return np.dot(self.A_der_x, f)
    # def derivative_x(self, f):
    #     df = (np.roll(f, -1, 0) - np.roll(f, +1, 0))/(2*self.dx)
    #     return df
    
    def derivative_z(self, f):
        return np.dot(f, self.A_der_z)
    # def derivative_z(self, f):
    #     df = (np.roll(f, -1, 1) - np.roll(f, -1, 1))/(2*self.dz)
    #     df[:,0] = (f[:,1] - f[:,0])/(self.dz)
    #     df[:,self.Nz] = (f[:,self.Nz] - f[:,self.Nz-1])/(self.dz)
    #     return df
    
    def derivative_xx(self, f):
        return np.dot(self.A_der_xx, f)

    def derivative_zz(self, f):
        return np.dot(f, self.A_der_zz)

    def Lapla(self, f):
        f_xx = self.derivative_xx(f)
        f_zz = self.derivative_zz(f)
        return f_xx + f_zz

    def BC(self, f, BC_L, BC_U):
        f_out = f
        f_out[:,0] = BC_L
        f_out[:,self.Nz] = BC_U
        return f_out

    def norm(self, f):
        return np.sqrt( (f**2).sum() / np.size(f) )
    
    def output_u(self, u, ff):
        for i in range(self.Nx):
            for k in range(self.Nz+1):
                ff.write(str(self.x[i]) + " " + \
                         str(self.z[k]) + " " + \
                         str(u[i, k]) + "\n")
            ff.write("\n")
        ff.write("\n")

    def output_vel(self, psi, ff):
        u = + self.derivative_z(psi)
        w = - self.derivative_x(psi)
        for i in range(self.Nx):
            for k in range(self.Nz+1):
                ff.write(str(self.x[i]) + " " + \
                         str(self.z[k]) + " " + \
                         str(u[i, k]) + " " + \
                         str(w[i, k]) + "\n")
        ff.write("\n\n")

    def output_T(self, T, ff):
        for i in range(self.Nx):
            for k in range(self.Nz+1):
                ff.write(str(self.x[i]) + " " + \
                         str(self.z[k]) + " " + \
                         str(T[i, k]) + "\n")
        ff.write("\n\n")
                
class Lapla_inv:

    def __init__(self):
        self.Set = Settings()
        
        self.mu_x = 1./( 2*(1 + (self.Set.dx/self.Set.dz)**2) )
        self.mu_z = 1./( 2*(1 + (self.Set.dz/self.Set.dx)**2) )
        self.mu_omega = 1./( 2*(self.Set.dx**(-2) \
                                + self.Set.dz**(-2)) )
        E_x = np.eye(self.Set.Nx)
        E_z = np.eye(self.Set.Nz+1)
        self.A_Li_x = self.mu_x*(np.roll(E_x, +1, 0) \
                                 + np.roll(E_x, -1, 0))
        self.A_Li_z = self.mu_z*(np.roll(E_z, +1, 0) \
                                 + np.roll(E_z, -1, 0))
        self.A_Li_z[:,0] = 0.
        self.A_Li_z[:,self.Set.Nz] = 0.

    def norm(self, f):
        return np.sqrt( (f**2).sum() / np.size(f) )

    def LaplaInv(self, omega, psi_old):
        '''\nabla^2 psi = omega を解く'''
        eps = 10.**(-3)
        psi = psi_old
        istop = 10000 # Error halt
        error = 1.
        i = 0
        while error > eps:
            for icount in range(10):
                psi = np.dot(self.A_Li_x, psi) \
                      + np.dot(psi, self.A_Li_z) \
                      - self.mu_omega * omega
                psi[:, 0] = 0.
                psi[:, self.Set.Nz] = 0.
            error = self.norm(omega - self.Set.Lapla(psi)) \
                    / self.norm(omega)
            i += 1
            if i >= istop:
                print("Program is halted: " + \
                      "in LaplaInv; " + \
                      "'error' failed to be smaller " +\
                      "than the given threshold.")
                f1 = open("E_Laplapsi.d", "w")
                f2 = open("E_psi.d", "w")
                f3 = open("E_omega.d", "w")

                self.Set.output_u(self.Set.Lapla(psi), f1)
                self.Set.output_u(psi, f2)
                self.Set.output_u(omega, f3)
                
                f1.close()
                f2.close()
                f3.close()
                sys.exit()
        return psi
                
class Diffusion:

    def __init__(self):
        self.Set = Settings()
        self.lmb_x = self.Set.dt / (2*self.Set.dx**2)
        self.lmb_z = self.Set.dt / (2*self.Set.dz**2)

    def init_diff_mat(self, D, BC_L, BC_U):
        A_x = self.init_mat_p(self.Set.Nx, D*self.lmb_x)
        A_z = self.init_mat_f(self.Set.Nz-1, D*self.lmb_z)
        B_x = self.init_mat_p(self.Set.Nx, - D*self.lmb_x)
        B_z = self.init_mat_f(self.Set.Nz-1, - D*self.lmb_z)
        U_BC = np.zeros((self.Set.Nx, self.Set.Nz-1))
        U_BC[:,0] = BC_L
        U_BC[:,self.Set.Nz-2] = BC_U

        A_1 = np.dot(B_x, np.linalg.inv(A_x))
        A_2 = np.dot(B_z, np.linalg.inv(A_z))
        A_3 = D*self.lmb_z \
              * np.dot(A_1 + np.eye(self.Set.Nx), \
                       np.dot(U_BC, np.linalg.inv(A_z)))
        return A_1, A_2, A_3
        
    def init_mat_p(self, n_size, x):
        E = np.eye(n_size)
        A = - x*np.roll(E, -1, 0) + (1 + 2*x)*E \
            - x*np.roll(E, +1, 0)
        return A
        
    def init_mat_f(self, n_size, x):
        A = self.init_mat_p(n_size, x)
        A[n_size-1,0] = 0.
        A[0,n_size-1] = 0.
        return A
        
    def evolve_diff(self, u, A_1, A_2, A_3):
        return np.dot(A_1, np.dot(u, A_2)) + A_3
    
        
class Benard_2D:

    def __init__(self):
        self.Set = Settings()
        self.Li = Lapla_inv()
        self.Df = Diffusion()

        # Diffusion process of omega
        self.A_omega_1, self.A_omega_2, self.A_omega_3 \
            = self.Df.init_diff_mat(self.Set.Pr, 0., 0.)

        # Diffusion process of T
        self.A_T_1, self.A_T_2, self.A_T_3 \
            = self.Df.init_diff_mat(1., self.Set.T_L, self.Set.T_U)

    def TT_trans_T_TTold(self, T, TT_old):
        return self.Li.LaplaInv(self.Set.derivative_x(T), TT_old)
        
    def OMEGA_trans_omega_TT(self, omega, TT):
        # hoge = open("test.d", "a")

        # self.Set.output_u(self.Set.Lapla(omega), hoge)
        
        # hoge.close()
        
        # OMEGA = (1. - self.Set.Pr)*self.Set.Lapla(omega) \
        #         + self.Set.Ra*self.Set.Pr*self.Set.derivative_x(T)
        OMEGA = (1. - self.Set.Pr)*omega \
                + self.Set.Ra*self.Set.Pr*TT
        return OMEGA

    def omega_trans_OMEGA_TT(self, OMEGA, TT):
        # Lap_omega = (OMEGA - \
        #              self.Set.Ra*self.Set.Pr \
        #              *self.Set.derivative_x(T)) \
        #              / (1. - self.Set.Pr)
        omega = (OMEGA - \
                 self.Set.Ra*self.Set.Pr*TT) \
                 / (1. - self.Set.Pr)
        return omega
        
    def Next_psi_omega_T_TT(self, psi, omega, T, TT):
        psi_mid, omega_mid, T_mid = \
                self.nonlin_Euler(psi, omega, T)
        TT_mid = self.TT_trans_T_TTold(T_mid, TT)
        OMEGA_mid = self.OMEGA_trans_omega_TT(omega_mid, TT_mid)
        OMEGA_new = self.Next_diff_u(OMEGA_mid)
        T_new = self.Next_diff_T(T_mid)
        TT_new = self.TT_trans_T_TTold(T_new, TT_mid)
        omega_new \
            = self.omega_trans_OMEGA_TT(OMEGA_new, TT_new)
        psi_new = self.Li.LaplaInv(omega_new, psi_mid)
        return psi_new, omega_new, T_new

    def Next_diff_u(self, u):
        u_new = np.zeros((self.Set.Nx, self.Set.Nz+1))
        u_new[:, 1:self.Set.Nz] \
            = self.Df.evolve_diff(\
                    u[:, 1:self.Set.Nz], \
                    self.A_omega_1, self.A_omega_2, self.A_omega_3)
        return u_new

    def Next_diff_T(self, T):
        T_new = np.zeros((self.Set.Nx, self.Set.Nz+1))
        T_new[:,0] = self.Set.T_L
        T_new[:,self.Set.Nz] = self.Set.T_U
        T_new[:,1:self.Set.Nz] \
            = self.Df.evolve_diff(\
                    T[:, 1:self.Set.Nz], \
                    self.A_T_1, self.A_T_2, self.A_T_3)
        return T_new
    
    def Jacobi(self, a, b):
        a_x = self.Set.derivative_x(a)
        a_z = self.Set.derivative_z(a)
        b_x = self.Set.derivative_x(b)
        b_z = self.Set.derivative_z(b)
        
        return a_z*b_x - a_x*b_z

    def nonlin_Euler(self, psi, omega, T):
        dlt_omega = - self.Set.dt * self.Jacobi(psi, omega)
        dlt_T = - self.Set.dt * self.Jacobi(psi, T)
        omega_new = omega + dlt_omega
        T_new = T + dlt_T
        psi_new = self.Li.LaplaInv(omega_new, psi)
        return psi_new, omega_new, T_new
        
    def initialize_psi_omega_T_TT(self):
        psi = np.sin(np.pi*self.Set.zz/self.Set.Lz) \
              * np.exp( np.cos(2*np.pi*self.Set.xx/self.Set.Lx) )
        psi = self.Set.BC(psi, 0., 0.)
        omega = self.Set.Lapla(psi)
        T = - self.Set.zz + 1 \
            + np.sin(np.pi*self.Set.zz/self.Set.Lz)
        T = self.Set.BC(T, self.Set.T_L, self.Set.T_U)
        TT = np.zeros((self.Set.Nx, self.Set.Nz+1))
        return psi, omega, T, TT

Set = Settings()
bc = Benard_2D()

psi, omega, T, TT = bc.initialize_psi_omega_T_TT()

lap_omega = Set.Lapla(omega)

ff = open("vel.d", "w")
gg = open("Temp.d", "w")
hh = open("norm_omega.d", "w")
kk = open("omega.d", "w")
kkkk = open("lap_omega.d", "w")

# bc.Set.output_u(psi, ff)
bc.Set.output_vel(psi, ff)
bc.Set.output_u(T, gg)
bc.Set.output_u(omega, kk)
bc.Set.output_u(lap_omega, kkkk)
hh.write(str(0.) + " " + str(bc.Li.norm(omega)) + "\n")
# bc.Set.output_u(omega, ff)

# sys.exit()

for it in range(1, bc.Set.Nt + 1):
    psi, omega, T = bc.Next_psi_omega_T_TT(psi, omega, T, TT)
    if it % bc.Set.tout == 0:
        t = it * bc.Set.dt
        print("it =" + str(it))
        bc.Set.output_vel(psi, ff)
        bc.Set.output_u(T, gg)
        bc.Set.output_u(omega, kk)
        hh.write(str(t) + " " + str(bc.Li.norm(omega)) + "\n")

ff.close()
gg.close()
hh.close()
