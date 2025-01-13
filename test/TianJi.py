import numpy as np
from itertools import permutations
from functools import partial
import matplotlib.pyplot as plt
import math


alpha = 1/np.array([9.80, 9.84, 9.89, 9.93, 9.95, 9.98])
beta = 1/np.array([9.81, 9.89, 9.91, 9.93, 9.94, 9.96])
v = np.array([1, 2, 3, 4, 5, 6])
# alpha = np.array([0, 1, 2])
# beta = np.array([1, 2, 3])
# v = np.array([1, 2, 3])
var = 1
n = len(v)
min_saddle_point_residual = float('inf')
ans_x = np.array([0] * int(math.factorial(n)))
ans_y = np.array([0] * int(math.factorial(n)))

# Note: Player 1 want to maximize xAy, and Player 2 want to minimize xAy

perms = np.array(list(permutations(range(len(v)))))


def val(x, y, u):
    if x > y:
        return u
    elif x == y:
        return 0
    else:
        return -u


def utility(x, y, variant):
    a = np.array([alpha[id] for id in x])
    b = np.array([beta[id] for id in y])
    util = np.array([val(a[i], b[i], v[i]) for i in range(len(x))])
    if variant == 1:
        return np.sum(util)
    else:
        return val(np.sum(util), 0, 1)


def calculate_matrix(perms, var):
    ans = np.zeros((len(perms),len(perms)))
    for i in range(len(perms)):
        for j in range(len(perms)):
            ans[i, j] = utility(perms[i], perms[j], var)
    return ans


def xA(x, A):
    return np.dot(x, A).flatten()


def Ay(A, y):
    yT = y.reshape(-1, 1)
    return np.dot(A, yT).flatten()


def optimal_player2(x, A): # player2 plays optimal under x
    return np.min(xA(x, A))


def optimal_player1(A, y): # player1 plays optimal under y
    return np.max(Ay(A, y))


def saddle_point_residual(x, A, y):
    return optimal_player1(A, y) - optimal_player2(x, A)


def attribute(arr):
    arr_non_negative = np.maximum(arr, 0)
    arr_sum = np.sum(arr_non_negative)
    if arr_sum == 0:
        return np.full_like(arr, 1 / len(arr), dtype=float)
    else:
        return arr_non_negative / arr_sum


A = calculate_matrix(perms, var)

# print(A)


class FP:
    def __init__(self, n):
        n = math.factorial(n)
        self.n = n
        self.count = 1
        self.x = np.array([1/n] * n)
        self.y = np.array([1/n] * n)
        self.sumx = np.array([1/n] * n)
        self.sumy = np.array([1/n] * n)

    def eval(self):
        spr = saddle_point_residual(self.x, A, self.y)
        global min_saddle_point_residual, ans_x, ans_y
        if spr < min_saddle_point_residual:
            min_saddle_point_residual = spr
            ans_x = self.x
            ans_y = self.y
        return spr

    def update(self):
        avgy = self.sumy/self.count
        Ay_vec = Ay(A, avgy)
        self.x = np.zeros_like(Ay_vec, dtype=float)
        max_index = np.argmax(Ay_vec)
        self.x[max_index] = 1
        self.sumx += self.x

        self.count += 1

        avgx = self.sumx/self.count
        xA_vec = xA(avgx, A)
        self.y = np.zeros_like(xA_vec, dtype=float)
        min_index = np.argmin(xA_vec)
        self.y[min_index] = 1
        self.sumy += self.y


class RM:
    def __init__(self, n):
        n = math.factorial(n)
        self.n = n
        self.x = np.array([1/n] * n)
        self.y = np.array([1/n] * n)
        self.rx = np.array([0] * n)
        self.ry = np.array([0] * n)

    def eval(self):
        spr = saddle_point_residual(self.x, A, self.y)
        global min_saddle_point_residual, ans_x, ans_y
        if spr < min_saddle_point_residual:
            min_saddle_point_residual = spr
            ans_x = self.x
            ans_y = self.y
        return spr

    def update(self):
        gx = Ay(A, self.y)
        lossx = np.sum(gx * self.x)
        self.rx = self.rx + gx - lossx
        self.x = attribute(self.rx)

        gy = xA(self.x, A)
        lossy = np.sum(gy * self.y)
        self.ry = self.ry + lossy - gy
        self.y = attribute(self.ry)


class RMP:
    def __init__(self, n):
        n = math.factorial(n)
        self.n = n
        self.x = np.array([1/n] * n)
        self.y = np.array([1/n] * n)
        self.rx = np.array([0] * n)
        self.ry = np.array([0] * n)

    def eval(self):
        spr = saddle_point_residual(self.x, A, self.y)
        global min_saddle_point_residual, ans_x, ans_y
        if spr < min_saddle_point_residual:
            min_saddle_point_residual = spr
            ans_x = self.x
            ans_y = self.y
        return spr

    def update(self):
        gx = Ay(A, self.y)
        lossx = np.sum(gx * self.x)
        self.rx = self.rx + gx - lossx
        self.rx = np.maximum(self.rx, 0)
        self.x = attribute(self.rx)

        gy = xA(self.x, A)
        lossy = np.sum(gy * self.y)
        self.ry = self.ry + lossy - gy
        self.ry = np.maximum(self.ry, 0)
        self.y = attribute(self.ry)


fp = FP(n)
rm = RM(n)
rmp = RMP(n)
x = [1]
y_fp = [fp.eval()]
y_rm = [rm.eval()]
y_rmp = [rmp.eval()]
for i in range(2, 1000):
    fp.update()
    rm.update()
    rmp.update()
    x.append(i)
    y_fp.append(fp.eval())
    y_rm.append(rm.eval())
    y_rmp.append(rmp.eval())
    if i % 50 == 0:
        # plt.figure(figsize=(10, 6))
        # plt.plot(x, y_fp, linestyle='-', color='b', linewidth=2, markersize=8)
        # plt.plot(x, y_rm, linestyle='-', color='r', linewidth=2, markersize=8)
        # plt.plot(x, y_rmp, linestyle='-', color='g', linewidth=2, markersize=8)
        # plt.yscale('log')
        # plt.title('Log Scale Line Plot', fontsize=16)
        # plt.xlabel('X-axis', fontsize=12)
        # plt.ylabel('Y-axis (Log Scale)', fontsize=12)
        # plt.grid(False)
        # plt.xticks(fontsize=10)
        # plt.yticks(fontsize=10)
        # plt.tight_layout()
        # plt.show()
        print("Min Saddle Point Residual: ", min_saddle_point_residual)
        print("Game Value for player 1: ", np.sum(ans_x * Ay(A, ans_y)))
    if min_saddle_point_residual < 0.001:
        break

# print(x)
# print(y_fp)
# print(y_rm)
# print(y_rmp)

plt.figure(figsize=(10, 6))
plt.plot(x, y_fp, linestyle='-', color='b', linewidth=2, markersize=8)
plt.plot(x, y_rm, linestyle='-', color='r', linewidth=2, markersize=8)
plt.plot(x, y_rmp, linestyle='-', color='g', linewidth=2, markersize=8)
plt.yscale('log')
plt.title('Log Scale Line Plot', fontsize=16)
plt.xlabel('X-axis', fontsize=12)
plt.ylabel('Y-axis (Log Scale)', fontsize=12)
plt.grid(False)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
output_dir = "/data/zec022/"
plt.savefig(os.path.join(output_dir, 'final_convergence_plot.png'), dpi=300, bbox_inches='tight')
print(min_saddle_point_residual)
print(ans_x)
print(ans_y)
print("Game Value for player 1: ", np.sum(ans_x * Ay(A, ans_y)))


