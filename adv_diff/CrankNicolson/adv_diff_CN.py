##########################################
# 一次元移流拡散方程式
#   \p_t u + c \p_x u = \nu \p_x^2 u
#   --- 周期境界条件
#   --- Crank-Nicolson 法により
#       移流項・拡散項を差分化
##########################################

import numpy as np

class adv_diff:

    def __init__(self):
        # 空間のパラメーター
        self.Lx = 50.0 # 領域サイズ
        self.Nx = 500 # 格子点数
        self.dx = self.Lx / self.Nx # メッシュ幅
        # 時間のパラメーター
        self.Nt = 100000 # 時間ステップ数
        self.dt = 0.01 # 時間刻み幅
        self.out = 100 # 出力ステップ周期
        # 方程式のパラメーター
        self.c = 1.0 # 移流速度
        self.nu = 0.01 # 拡散速度
        # 空間座標
        self.x = np.linspace(0, self.Lx, self.Nx, endpoint=False)

    def coefficient(self):
        lmb = self.c*self.dt / (4*self.dx)
        mu  = self.nu*self.dt / (4*self.dx**2)
        E = np.eye(self.Nx)
        E_p = np.roll(E, +1, axis=1)
        E_m = np.roll(E, -1, axis=1)
        A_coeff = \
            (lmb - mu)*E_p + (1 + 2*mu)*E + (- lmb - mu)*E_m
        B_coeff = \
            (- lmb + mu)*E_p + (1 - 2*mu)*E + (lmb + mu)*E_m
        return np.dot(np.linalg.inv(A_coeff), B_coeff)

    def Next_u(self, u, A):
        return np.dot(A, u)

    def initial_u(self):
        # return np.sin(2*np.pi*self.x / self.Lx)
        return - 1 + self.x / 25

ad = adv_diff()
u = ad.initial_u()
A = ad.coefficient()

f = open("adv_diff.dat","w")

for ix in range(0, ad.Nx):
    f.write(str(ad.x[ix]) + " " + str(u[ix]) + "\n")
f.write("\n\n")
    
for it in range(1, ad.Nt + 1):
    u = ad.Next_u(u, A)
    if it % ad.out == 0:
        for ix in range(0, ad.Nx):
            f.write(str(ad.x[ix]) + " " + str(u[ix]) + "\n")
        f.write("\n\n")
