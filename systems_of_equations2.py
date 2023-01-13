import sympy as sy
import matplotlib.pyplot as plt
import numpy as np
import time


u, v = sy.symbols('u, v')

def gaussian_elimination(e0, b):
    e = np.array(e0)
    n = len(e)
    """想好有几个循环过程"""
    for i in range(n-1):
        for j in range(i+1, n):
            a = e[j][i] / e[i][i]
            for k in range(i, n):
                e[j][k] = e[j][k] - e[i][k] * a
            b[j] -= b[i] * a
    x = [0] * n
    for i in range(n-1, -1, -1):
        x0 = 0
        for j in range(n-1, i-1, -1):
            x0 += x[j] * e[i][j]
        x[i] = (b[i]- x0) / e[i][i]


    return x

"""用于检查矩阵是否是严格对角占优矩阵"""
def test_diagonally(A):
    n = len(A)
    for i in range(n):
        if abs(A[i, i]) < abs(np.sum(A[i, :])):
            print("this matrix is not strictly diagonally domination. ")
            break

def lu_factorization(e0, b):
    pass

def jacobi(A, x0, b, k):

    n = len(A)
    A = np.array(A)
    test_diagonally(A)
    x0 = np.array(x0)
    b = np.array(b)

    x_list = np.zeros([k+1, n], dtype=float)
    L = np.zeros([n, n])
    U = np.zeros([n, n])
    LU = np.zeros([n, n])
    D_inver = np.zeros(n)
    D = np.zeros(n)


    x_list[0] = x0

    if n!=len(A[0]):
        print('A is not a square matrix')
    else:
        for i in range(n):
            for j in range(n):
                if i == j:
                    D[i] = A[i][j]
                    D_inver[i] = 1. / A[i][j]
                elif i < j:
                    U[i][j] = A[i][j]
                elif i > j:
                    L[i][j] = A[i][j]
                LU[i][j] = L[i][j] + U[i][j]

        for m in range(1, k+1):
            for i in range(0, n):
                temp = np.zeros(n)
                for j in range(0, n):
                    temp[i] += LU[i][j] * x_list[m-1][j]
                x_list[m][i] = D_inver[i] * (b[i] - temp[i])
    return x_list[k]

def gauss_seidel(A, x0, b, k):
    n = len(A)
    A = np.array(A)
    test_diagonally(A)
    x0 = np.array(x0)
    b = np.array(b)

    x_list = np.zeros([k + 1, n], dtype=float)
    L = np.zeros([n, n])
    U = np.zeros([n, n])

    D_inver = np.zeros(n)
    D = np.zeros(n)

    x_list[0] = x0

    if n != len(A[0]):
        print('A is not a square matrix')
    else:
        for i in range(n):
            for j in range(n):
                if i == j:
                    D[i] = A[i][j]
                    D_inver[i] = 1. / A[i][j]
                elif i < j:
                    U[i][j] = A[i][j]
                elif i > j:
                    L[i][j] = A[i][j]


        for m in range(1, k + 1):
            for i in range(0, n):
                temp = np.zeros(n)
                for j in range(0, n):
                    temp[i] += U[i][j] * x_list[m - 1][j] + L[i][j] * x_list[m][j]

                    x_list[m][i] = D_inver[i] * (b[i] - temp[i])
                    """这里有点巧妙，得多看看"""
    return x_list[k]

def Sor(omiga, A, x0, b, k):
    n = len(A)
    A = np.array(A)
    test_diagonally(A)
    x0 = np.array(x0)
    b = np.array(b)

    x_list = np.zeros([k + 1, n], dtype=float)
    x_list[0] = x0
    temp_x_list = np.zeros([k + 1, n], dtype=float)
    temp_x_list[0] = x0
    D = np.zeros(n)
    ng_omiga_UD = np.zeros([n, n])
    omiga_L = np.zeros([n, n])
    omiga_b = np.zeros(n)

    if n != len(A[0]):
        print('A is not a square matrix')
    else:
        """将A拆解成三个矩阵。不过这不是单纯的拆解，为了后面方便计算，omiga也被乘进去了"""
        for i in range(n):
            omiga_b[i] = omiga * b[i]
            for j in range(n):
                if i == j:
                    D[i] = A[i][j]
                elif i < j:
                    ng_omiga_UD[i][j] = -omiga * A[i][j]
                elif i > j:
                    omiga_L[i][j] = omiga * A[i][j]
        """计算omiga*L+D 和-U+(1-omiga)*D"""
        for i in range(n):
            for j in range(n):
                if i == j:
                    omiga_L[i][j] = D[i]
                    ng_omiga_UD[i][j] = (1. - omiga) * D[i]
        """求出这个的逆矩阵，就相当于解方程组了"""
        temp_omiga_L = np.zeros([n, n])
        inver_omiga_LD = np.zeros([n, n], dtype=float)
        for m in range(n):
            temp_e = np.zeros(n)
            temp_e[m] = 1.
            for i in range(n):
                x0 = 0.
                for j in range(n):
                    x0 += temp_omiga_L[m][j] * omiga_L[i][j]
                temp_omiga_L[m][i] = (temp_e[i] - x0) / omiga_L[i][i]
        """需要转置一下才能成为真正的逆矩阵"""
        for i in range(n):
            for j in range(n):
                inver_omiga_LD[i][j] = temp_omiga_L[j][i]
        """计算xk"""
        """分成两步，先把xk和前面的乘起来，在加上omiga_b；下面的iijj循环是为了跟外面的逆矩阵相乘"""
        for m in range(1, k+1):
            temp = np.zeros(n)
            for i in range(n):
                for j in range(n):
                    temp[i] += ng_omiga_UD[i][j] * x_list[m - 1][j]
                temp[i] += omiga_b[i]
            for ii in range(n):
                for jj in range(n):
                    x_list[m][ii] += inver_omiga_LD[ii][jj] * temp[jj]
    return x_list[k]

def QR(A):
    m = len(A)
    n = len(A[0])
    A = np.array(A)

    Q = np.zeros((m, n), dtype=float)
    R = np.zeros((n, n))
    Y = np.zeros((m, n), dtype=float)

    Y[:, 0] = A[:, 0]
    sum0 = 0
    for i in range(m):
        sum0 += Y[i][0] ** 2
    Q[:, 0] = Y[:, 0] / np.sqrt(sum0)
    R[0][0] = np.sqrt(sum0)
    """怎么想到将计算R和Y的部分拆开来，如何将R和Y的迭代联系到一起是一个问题"""
    """首先R和Y的迭代列是一致的，那么就可以先把i提取出来作为列的迭代。
    Y的计算要等到R计算完才能算，所以分成两个。j的取值就是要计算的该列共有多少行。k则是需要多少个加法。"""
    for i in range(1, n):
        Y[:, i] = A[:, i]
        sum0 = 0
        for j in range(i):
            for k in range(m):
                R[j, i] += Q[k, j] * A[k, i]

        for j in range(m):
            for k in range(i):
                Y[j, i] -= Q[j, k] * R[k, i]
            sum0 += Y[j, i] ** 2
        R[i][i] = np.sqrt(sum0)
        Q[:, i] = Y[:, i] / R[i][i]
    return Q, R

"""为了提高QR的程序局部性"""
def QR1(A):
    m = len(A)
    n = len(A[0])
    A = np.array(A, dtype='float64')

    A_t = np.zeros((n, m), dtype='float64')
    Q = np.zeros((n, m), dtype='float64')
    R = np.zeros((n, n), dtype='float64')
    Y = np.zeros((n, m), dtype='float64')

    for i in range(n):
        for j in range(m):
            A_t[i, j] = A[j, i]
    Y[0, :] = A_t[0, :]
    sum0 = 0
    for i in range(m):
        sum0 += Y[0, i] ** 2
    R[0, 0] = np.sqrt(sum0)
    Q[0, :] = Y[0, :] / R[0, 0]
    for i in range(1, n):
        Y[i, :] = A_t[i, :]
        sum0 = 0
        for j in range(i):
            for k in range(m):
                R[i, j] += Q[j, k] * A_t[i, k]
        for j in range(m):
            for k in range(i):
                Y[i, j] -= Q[k, j] * R[i, k]
            sum0 += Y[i, j] ** 2
        R[i, i] = np.sqrt(sum0)
        Q[i, :] = Y[i, :] / R[i, i]

    return Q, R


"""改进的格拉姆-施密特正交方法"""
"""求一个R就要求一个Y，不能再向原来那样全部求完R，再求Y了"""
def QR_(A):
    m = len(A)
    n = len(A[0])
    A = np.array(A, dtype='float64')

    A_t = np.zeros((n, m), dtype='float64')
    Q = np.zeros((n, m), dtype='float64')
    R = np.zeros((n, n), dtype='float64')
    Y = np.zeros((n, m), dtype='float64')

    for i in range(n):
        for j in range(m):
            A_t[i, j] = A[j, i]
    Y[0, :] = A_t[0, :]
    sum0 = 0
    for i in range(m):
        sum0 += Y[0, i] ** 2
    R[0, 0] = np.sqrt(sum0)
    Q[0, :] = Y[0, :] / R[0, 0]
    for i in range(1, n):
        Y[i, :] = A_t[i, :]
        sum0 = 0
        for j in range(i):
            for k in range(m):
                R[i, j] += Q[j, k] * Y[i, k]
            Y[i, :] -= Q[j, :] * R[i, j]
        for j in range(m):
            sum0 += Y[i, j] ** 2
        R[i, i] = np.sqrt(sum0)
        """真该仔细的检查每一步骤"""
        Q[i, :] = Y[i, :] / R[i, i]

    return Q, R


def newton_method(f1, f2, DF):
    pass








"""试试雅克比迭代怎么样， 精确解是2，-1,1"""
def test0():
    A = np.array([[3, 1, -1], [2, 4, 1], [-1, 2, 5]])
    x0 = np.array([0, 0, 0])
    b = np.array([4, 1, 1])
    x = jacobi(A, x0, b, 50)
    print(x)

def test1():
    A = np.array([[3, 1, -1], [2, 4, 1], [-1, 2, 5]])
    x0 = np.array([0, 0, 0])
    b = np.array([4, 1, 1])
    x = gauss_seidel(A, x0, b, 21)
    print(x)

def test2():
    """当k为14时收敛到精确解，但是呢omiga也是可以变得不知道该怎么取值,似乎为1.25就是比较好的了"""
    A = np.array([[3, 1, -1], [2, 4, 1], [-1, 2, 5]])
    x0 = np.array([0, 0, 0])
    b = np.array([4, 1, 1])
    x = Sor(1.25, A, x0, b, 14)
    print(x)

"""试试QR分解"""
def test3():
    A = [[1., -4], [2., 3.], [2., 2.]]
    a = QR(A)
    Q = a[0]
    R = a[1]

"""验证一下普通QR方法会导致的不正交"""
def test4():
    delta = 1.e-20
    A = [[1., -4], [2., 3.], [2., 2.]]
    A1 = [[1., 1., 1.], [delta, 0., 0.], [0., delta, 0.], [0., 0., delta]]
    a = QR(A1)
    Q = a[0]
    R = a[1]

    a1 = QR_(A1)
    Q1 = a1[0]
    R1 = a1[1]

    sum0 = 0.
    sum1 = 0.
    for i in range(len(A)):
        sum0 += round(Q[i, -2], 15) * round(Q[i, -1], 15)
        sum1 += round(Q1[-2, i], 15) * round(Q1[-1, i], 15)
    print(sum0)
    print(sum1)

"""试试程序的局部性"""
"""果然还是步长小的局部性好"""
def test5():
    A = [[1., -4], [2., 3.], [2., 2.]]
    delta = 1.e-20
    A1 = [[1., 1., 1.], [delta, 0., 0.], [0., delta, 0.], [0., 0., delta]]

    start_time1 = time.perf_counter()
    a = QR(A1)
    end_time1 = time.perf_counter()

    start_time2 = time.perf_counter()
    a1 = QR1(A1)
    end_time2 = time.perf_counter()

    print('time1=', end_time1-start_time1)
    print('time2=', end_time2 - start_time2)







if __name__ == "__main__":
    test5()