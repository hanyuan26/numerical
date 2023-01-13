import sympy as sy
import numpy as np
import matplotlib.pyplot as plt


x, t = sy.symbols('x, t')
f = sy.sin(2 * sy.pi * x) ** 2
l = 0. * t
r = 0. * t
def parapolic_forth(f, xl, xr, yb, yt, M, N):


    h = (xr - xl) / M
    k = (yt - yb) / N
    D = 1
    sigma = D * k / h ** 2

    m = M -1
    n = N -1
    w = np.zeros((N+1, m))
    x_list = [(xl + h * i) for i in range(M+1)]
    t_list = [(yb + k * j) for j in range(N+1)]

    A0 = [1 - 2 * sigma] * m
    A1 = [sigma] * (m - 1)
    A = np.diag(A0) + np.diag(A1, 1) + np.diag(A1, -1)
    l_side = [0] * (N+1)
    r_side = [0] * (N+1)
    for i in range(1, m):
        w[0, i] = f.evalf(subs={'x': x_list[i]})
    for i in range(0, N):
        tmp1 = [0.] * (m-2)
        tmp1.insert(0, l_side[i])
        tmp1.append(r_side[i])
        tmp1 = np.array(tmp1)
        s = sigma * tmp1

        w[i+1, :] = np.matmul(A, w[i, :]) + s
        a = np.reshape(l_side, (N+1, 1))
    w_new0 = np.concatenate((np.reshape(l_side, (N+1, 1)), w), axis=1)
    w_new1 = np.concatenate((w_new0, np.reshape(r_side, (N+1, 1))), axis=1)

    fig = plt.figure()
    ax = plt.gca(projection='3d')
    xx, tt = np.meshgrid(x_list, t_list)
    ax.plot_surface(xx, tt, w_new1, cmap='jet')
    plt.show()


def case0():
    a = parapolic_forth(f, 0, 1, 0, 1, 10, 250)
    print(a)

case0()

