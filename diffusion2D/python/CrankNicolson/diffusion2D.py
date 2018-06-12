#########################################
# 2次元拡散方程式
#   \p_t u = D \nabla^2 u
# の数値積分
#    - Crank-Nicolson 法
#########################################
import numpy as np

class diffusion:

    def __init__(self):
        self.Lx, self.Ly = 10.0, 10.0
        self.Nx, self.Ny = 50, 50
        self.NxNy = self.Nx * self.Ny
        self.dx, self.dy = self.Lx/self.Nx, self.Ly/self.Ny
        self.dt = 0.01
        self.D = 4.0
        self.out = 10
        self.tend = 100

    def matrix_unit(self,n,i,j):
        A = np.zeros((n,n))
        A[i,j] = 1.0
        return A
        
    def init_mesh(self):
        x = np.linspace(0.0, self.Lx, self.Nx, endpoint=False)
        y = np.linspace(0.0, self.Ly, self.Ny, endpoint=False)
        xx, yy = np.meshgrid(x, y)        
        return xx, yy

    def increment_decrement_x(self):
        E_mini = np.eye(self.Nx)
        A_inc_mini = np.roll(E_mini, +1, axis=1)
        A_dec_mini = np.roll(E_mini, -1, axis=1)
        A_inc = np.zeros((self.NxNy, self.NxNy))
        A_dec = np.zeros((self.NxNy, self.NxNy))
        print("ok!")
        for iy in range(self.Ny):
            iyy = iy*self.Nx
            A_inc[iyy : iyy+self.Nx, iyy : iyy+self.Nx] = A_inc_mini
            A_dec[iyy : iyy+self.Nx, iyy : iyy+self.Nx] = A_dec_mini
        print("great!")
        return A_inc, A_dec

    def increment_decrement_y(self):
        E = np.eye(self.NxNy)
        A_inc = np.roll(E, + self.Nx, axis=1)
        A_dec = np.roll(E, - self.Nx, axis=1)
        print("cool!")
        return A_inc, A_dec
    
    def advance_matrix(self):
        # A (u_new) = B (u_old)
        A_x_inc, A_x_dec = self.increment_decrement_x()
        A_y_inc, A_y_dec = self.increment_decrement_y()
        lmb = self.D*self.dt / (4*self.dx**2)
        E = np.eye(self.NxNy)
        print("wonderful!")
        A = (1 + 4*lmb)*E - lmb*(A_x_inc + A_x_dec + A_y_inc + A_y_dec)
        B = (1 - 4*lmb)*E + lmb*(A_x_inc + A_x_dec + A_y_inc + A_y_dec)
        print("amazing!")
        A_inv = np.linalg.inv(A)
        print("crazy!")
        return np.dot(A_inv, B)

    def Next_u(self, u, A):
        return np.dot(A, u)

    def u_init(self):
        u_0 = np.zeros((self.Nx, self.Ny))
        u_0[self.Nx//2,self.Ny//2] = 1.0
        return u_0


f = open("diff.dat","w")

dif = diffusion()
A = dif.advance_matrix()
u = dif.u_init()
xx, yy = dif.init_mesh()

for iy in range(dif.Ny):
    for ix in range(dif.Nx):
        f.write(str(xx[ix,iy]) + " " + \
                str(yy[ix,iy]) + " " + \
                str(u[ix,iy]) + "\n")
    f.write("\n")
f.write("\n")
print("written!")
                
u_line = np.reshape(u, dif.NxNy)

for it in range(1,dif.tend+1):
    u_line = dif.Next_u(u_line, A)
    if it % dif.out == 0:
        u = np.reshape(u_line, (dif.Ny, dif.Nx))
        for iy in range(dif.Ny):
            for ix in range(dif.Nx):
                f.write(str(xx[ix,iy]) + " " + \
                        str(yy[ix,iy]) + " " + \
                        str(u[ix,iy]) + "\n")
            f.write("\n")
        f.write("\n")
        print("written!")

f.close()
