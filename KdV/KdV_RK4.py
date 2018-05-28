import numpy as np
import sys


def derivative_x(u, dx): # u は np.array を想定
    u_plus  = np.roll(u, -1)
    u_minus = np.roll(u, +1)
    du = (u_plus - u_minus)/(2*dx)
    return du

def derivative_xxx(u, dx): # u は np.array を想定
    u_pp = np.roll(u, -2)
    u_p  = np.roll(u, -1)
    u_m  = np.roll(u, +1)
    u_mm = np.roll(u, +2)
    dddu = (u_pp - 2*u_p + 2*u_m - u_mm)/(2*dx**3)
    return dddu

def F(u, dx):
    u_x = derivative_x(u, dx)
    u_xxx = derivative_xxx(u, dx)
    # if(u[0] != u[0]):
    #     sys.exit()
    return - 6*u*u_x - u_xxx

def RK4(u, dt):
    k1 = F(u, dx)
    k2 = F(u + k1*dt/2, dx)
    k3 = F(u + k2*dt/2, dx)
    k4 = F(u + k3*dt, dx)
    u_new = u + (k1 + 2*k2 + 2*k3 + k4)*dt/6
    return u_new

L = 30.0
N = 300
dx = L/N
dt = 0.001
M = 20000
out = 100
# count = 0



x = np.linspace(0, L, N, endpoint = False)
u = np.sin(2*np.pi*x/L) # 初期条件

# f = open("outdata/u0.dat","w")
# for ix in range(0, N):
#     f.write(str(x[ix]) + " " + str(u[ix]) + "\n")
# f.close()

f_anim = open("outdata/u_anim.dat","w")
for ix in range(0, N):
    f_anim.write(str(x[ix]) + " " + str(u[ix]) + "\n")
f_anim.write("\n\n")


for i in range(1, M+1):
    u = RK4(u, dt)
    if (i % out == 0):
        # count += 1
        # f = open("outdata/u" + str(count) + ".dat","w")
        # for ix in range(0, N):
        #     f.write(str(x[ix]) + " " + str(u[ix]) + "\n")
        # f.close()

        for ix in range(0,N):
            f_anim.write(str(x[ix]) + " " + str(u[ix]) + "\n")
        f_anim.write("\n\n")

f_anim.close()
