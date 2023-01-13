import sympy as sy
import matplotlib.pyplot as plt
import numpy as np
import time


from systems_of_equations2 import Sor


x = sy.symbols('x')

def fitting(a, b, kind=1):
    n = len(a)
    a = np.array(a)
    b = np.array(b)
    m = kind+1
    """初始化矩阵A"""
    A = np.ones((n, m), dtype=float)
    for i in range(n):
        for j in range(1, m):
            A[i][j] = a[i] ** j
    """计算对称矩阵"""
    AA = np.zeros((m, m))
    sum = 0
    for i in range(m):
        for j in range(m):
            for k in range(n):
                sum += A[k][i] * A[k][j]
            AA[i][j] = sum
            sum = 0
    """计算等号后面那个东西"""
    """记住这里sum的处理方式"""
    A_b = np.zeros(m)
    sum = 0
    for i in range(m):
        for k in range(n):
            sum += A[k][i] * b[k]
        A_b[i] = sum
        sum = 0

    x0 = np.zeros(m)
    c = Sor(1.2, AA, x0, A_b, 20)

    """计算误差"""
    r = np.zeros(n)
    for i in range(n):
        sum = 0
        for j in range(m):
            sum += A[i][j] * c[j]
            r[i] = b[i] - sum
    SE = 0
    for i in range(n):
        SE += r[i] ** 2

    RMSE = np.sqrt(SE/n)



    return c

def test0():
    a = [1., -1., 1.]
    b = [2., 1., 3.]
    c = fitting(a, b)

    f = c[0] + c[1] * x
    x_list = np.linspace(-1.5, 3, 20)
    y_list = np.zeros(len(x_list))
    for i in range(len(x_list)):
        y_list[i] = f.evalf(subs={'x': x_list[i]})

    plt.scatter(a, b)
    plt.plot(x_list, y_list, '-.')
    plt.show()


if __name__ == "__main__":
    test0()
