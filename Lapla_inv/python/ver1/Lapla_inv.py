#######################################
# 空間2次元の Poisson 方程式
# \nabla^2 u = f
# を Jacobi 法で解く.
#######################################


import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

class Lapla_inv:

    def __init__(self):
        # 空間パラメーター
        self.Lx, self.Lz = 50., 50.
        self.Nx, self.Nz = 500, 500
        self.dx, self.dz = self.Lx/self.Nx, self.Lz/self.Nz
        self.x = np.linspace(0., self.Lx, self.Nx, endpoint=False)
        self.z = np.linspace(0., self.Lz, self.Nz+1, endpoint=True)
        self.zz, self.xx = np.meshgrid(self.z, self.x)
        # 時間パラメーター
        self.Nt = 100
        self.dt = 1/( 4*(self.dx**(-2) + self.dz**(-2)) )
        # 時間前進の行列
        self.lmb_x = self.dt/(self.dx**2)
        self.lmb_z = self.dt/(self.dz**2)
        E_x, E_z = np.eye(self.Nx), np.eye(self.Nz+1)
        self.A_x = (np.roll(E_x, +1, 0) + np.roll(E_x, -1, 0))/\
                   (2*(1 + (self.dx/self.dz)**2))
        self.A_z = (np.roll(E_z, +1, 0) + np.roll(E_z, -1, 0))/\
                   (2*(1 + (self.dz/self.dx)**2))
        self.lmb = 1/( 2*(self.dx**(-2) + self.dz**(-2)) )

    def init_f(self):
        # f = np.sin(4*np.pi*self.xx/self.Lx) *\
        #     np.sin(2*np.pi*self.zz/self.Lz)
        f = np.sin(4*np.pi*self.zz/self.Lz)
        return f

    def BC(self, u):
        u_out = u
        u_out[:,0] = 0.
        u_out[:,self.Nz] = 0.
        return u_out
    
    def u_next_u(self, u, f):
        u_next = - self.lmb * f \
                 + np.dot(self.A_x, u) \
                 + np.dot(u, (self.A_z).T)
        return self.BC(u_next)

    def u_Lapla_u(self, u):
        u_xx = (np.roll(u, +1, 1) - 2*u + np.roll(u, -1, 1))/(self.dx**2)
        u_zz = (np.roll(u, +1, 0) - 2*u + np.roll(u, -1, 0))/(self.dz**2)
        return u_xx + self.BC(u_zz)
    

lp = Lapla_inv()
f_rhs = lp.init_f()


u = np.zeros((lp.Nx, lp.Nz+1))
for i in range(1000):
    u = lp.u_next_u(u, f_rhs)

fig = plt.figure()
ax = Axes3D(fig)
# ax.plot_wireframe(lp.xx, lp.zz, f_rhs - lp.u_Lapla_u(u))
# ax.plot_wireframe(lp.xx, lp.zz, f_rhs)
# ax.plot_wireframe(lp.xx, lp.zz, lp.u_Lapla_u(u))
ax.plot_wireframe(lp.xx, lp.zz, u)
print(u[5,:])

plt.show()
