import sympy as sy
import matplotlib.pyplot as plt
import numpy as np
import time


from systems_of_equations2 import Sor
from systems_of_equations2 import gauss_seidel

x = sy.symbols('x')

"""
这是一个拉格朗日插值。
输入：x,y对应好的两组列向量。
输出：插值多项式
"""
def lagrange_interpolation(a, b):
    n = len(a)
    l, l1 = 0., 0.
    for i in range(n):
        xi = 1.
        xk = 1.
        for j in range(n):

            if i != j:
                xi *= x - a[j]
                xk *= a[i] - a[j]

        l += b[i] * (xi / xk)
        l1 = sy.factor(l)
    return l1

"""这个专门是给高斯积分时候计算系数用的，要不对于上面那个，我觉得可以把输入改一下b可以作为一个默认值。"""
def lagrange_interpolation_(x_list):
    n = len(x_list)
    l_list = []
    for i in range(n):
        xi = 1.
        xk = 1.
        for j in range(n):

            if i != j:
                xi *= x - x_list[j]
                xk *= x_list[i] - x_list[j]

        l_list.append(xi / xk)
    return l_list
"""做一个《有限元法》p430的例子。"""
def fem_case1():
    a = [0, 1, 3]
    b = [1, 1, 5]
    p = lagrange_interpolation(a, b)
    print(p)

def newton_divided_differences(a,b):
    n = len(a)
    e1 = [[0.] * n] * n
    e1[0] = b
    e = np.array(e1)
    for i in range(1, n):
        for j in range(0, n - i):
            e0 = e[i-1, j+1] - e[i-1, j]
            a0 = a[j+i] - a[j]
            e[i, j] = e0 / a0

    x0 = []
    for i in range(n):
        x00 = e[i, 0]
        x0.append(x00)

    y = x0[0]
    for i in range(1, n):
        yi = x0[i]
        for j in range(i):
            yi *= x - a[j]
        y += yi

    y1  = sy.factor(y)
    return y1

def chebyshev_interpolation(a, n):
    """用这个东西只能选择插值点啊，难不成选完之后我还得按照这些点去做实验找数据不成？"""
    x_list = []
    for i in range(1, n+1):
        cs = sy.cos(sy.pi*(2*i-1)/(2*n))
        x0 = 0.5*(a[0]+a[-1])+0.5*(a[-1]-a[0])*cs
        xi = float(x0)
        x_list.append(xi)
    return x_list

def bezier_curves(a0):
    """t的取值只能是【0,1】,
    第一个控制点和第一个点的斜率就是端点的斜率，第2个控制点和最后一点的斜率是右端点的斜率"""
    t =sy.symbols('t')
    n = len(a0)
    a = np.array(a0)
    Sx = []
    Sy = []
    for i in range(n):
        bx = 3. * (a[i, 2] - a[i, 0])
        cx = 3. * (a[i, 4] - a[i, 2]) - bx
        dx = a[i, 6] - a[i, 0] - bx - cx

        by = 3. * (a[i, 3] - a[i, 1])
        cy = 3. * (a[i, 5] - a[i, 3]) - by
        dy = a[i, 7] - a[i, 1] - by - cy

        x = a[i, 0] + bx * t + cx * t ** 2. + dx * t ** 3.
        y = a[i, 1] + by * t + cy * t ** 2. + dy * t ** 3.
        Sx.append(x)
        Sy.append(y)

    return Sx, Sy

"""或许可以给他试试函数修饰符，把这个类尘封起来吧"""
class CubicSplines():
    def __init__(self, a, b):
        self.n = len(a)
        self.a = np.array(a)
        self.b = np.array(b)
        self.A = np.zeros((self.n, self.n))
        self.delta_b = np.zeros(self.n)

    def cal(self):
        self.delta_list = self.delta()
        self.Delta_list = self.Delta()
        self.matrix()

    """计算delta"""

    def delta(self):
        delta_list = np.zeros(self.n - 1)
        for i in range(1, self.n):
            delta_list[i - 1] = self.a[i] - self.a[i - 1]
        return delta_list

    def Delta(self):
        Delta_list = np.zeros(self.n - 1)
        for i in range(1, self.n):
            Delta_list[i - 1] = self.b[i] - self.b[i - 1]
        return Delta_list

    """计算矩阵和向量"""

    def matrix(self):
        for i in range(1, self.n - 1):
            for j in range(i - 1, i + 2):
                if i == j:
                    self.A[i][j] = 2 * (self.delta_list[i - 1] + self.delta_list[i])
                elif i < j:
                    self.A[i][j] = self.delta_list[i]
                elif i > j:
                    self.A[i][j] = self.delta_list[i - 1]
        for i in range(1, self.n - 1):
            self.delta_b[i] = 3 * (
                        self.Delta_list[i] / self.delta_list[i] - self.Delta_list[i - 1] / self.delta_list[i - 1])

    def cal_plot(self):
        """使用迭代法解线性方程组 求解系数"""
        x0 = np.zeros(self.n)
        c = Sor(1.2, self.A, x0, self.delta_b, 20)
        d_list = np.zeros(self.n - 1)
        b_list = np.zeros(self.n - 1)
        for i in range(1, self.n):
            d_list[i - 1] = (c[i] - c[i - 1]) / (3. * self.delta_list[i - 1])
            b_list[i - 1] = (self.Delta_list[i - 1] / self.delta_list[i - 1]) - (self.delta_list[i - 1] / 3) * (
                        2 * c[i - 1] + c[i])
        """求出每段方程"""
        S = []
        for i in range(self.n - 1):
            s = self.b[i] + b_list[i] * (x - self.a[i]) + c[i] * (x - self.a[i]) ** 2 + d_list[i] * (x - self.a[i]) ** 3
            S.append(s)
        """绘图"""
        m = 10
        plt.scatter(self.a, self.b)
        for i in range(self.n - 1):
            y_list = np.zeros(m)
            x_list = np.linspace(self.a[i], self.a[i + 1], m)
            for j in range(len(x_list)):
                y_list[j] = S[i].evalf(subs={'x': x_list[j]})
            plt.plot(x_list, y_list, '-.c')
        plt.show()

    def cubic_splines(self):
        """求解矩阵"""
        self.A[0][0] = 1.
        self.A[self.n - 1][self.n - 1] = 1.

    def curvature_adjust(self, v1, v2):
        """设置左右端点的2阶导数曲率"""
        self.A[0][0] = 2.
        self.A[self.n - 1][self.n - 1] = 2.
        self.delta_b[0] = v1
        self.delta_b[self.n - 1] = v2

    def clamp(self, v1, vn):
        """设置左右端点的1阶导数曲率"""
        self.A[0][0] = 2 * self.delta_list[0]
        self.A[0][1] = self.delta_list[0]
        self.A[-1][-1] = 2 * self.delta_list[-1]
        self.A[-1][-2] = self.delta_list[-1]

        self.delta_b[0] = 3. * (self.Delta_list[0] / self.delta_list[0] - v1)
        self.delta_b[-1] = 3. * (vn - self.Delta_list[-1] / self.delta_list[-1])

    def parabolically_terminated(self):
        """第1个和最后一个曲线设置为抛物线"""
        self.A[0, 0:2] = np.array([1., -1.])
        self.A[-1][-2:] = np.array([-1, 1.])

    def not_a_kont(self):
        """第1,2和倒数两个曲线分别设置为同一条曲线"""
        self.A[0][0] = self.delta_list[1]
        self.A[0][1] = -(self.delta_list[0] + self.delta_list[1])
        self.A[0][2] = self.delta_list[0]
        self.A[self.n - 1][self.n - 1] = self.delta_list[-2]
        self.A[self.n - 1][self.n - 2] = -(self.delta_list[-2] + self.delta_list[-1])
        self.A[self.n - 1][self.n - 3] = self.delta_list[-1]

def cubic_splines(a, b, kind='a', *v):
    n = len(a)
    a = np.array(a)
    b = np.array(b)

    """计算delta"""
    delta_list = np.zeros(n - 1)
    Delta_list = np.zeros(n - 1)
    for i in range(1, n):
        delta_list[i - 1] = a[i] - a[i - 1]
        Delta_list[i - 1] = b[i] - b[i - 1]
    """计算矩阵和向量"""
    A = np.zeros((n, n))
    delta_b = np.zeros(n)
    for i in range(1, n - 1):
        for j in range(i - 1, i + 2):
            if i == j:
                A[i][j] = 2 * (delta_list[i - 1] + delta_list[i])
            elif i < j:
                A[i][j] = delta_list[i]
            elif i > j:
                A[i][j] = delta_list[i - 1]
    for i in range(1, n - 1):
        delta_b[i] = 3 * (Delta_list[i] / delta_list[i] - Delta_list[i - 1] / delta_list[i - 1])

    """根据不同的要求改变矩阵A和向量delta_b"""
    """a自然样条, b设置左右端点的2阶导数曲率, c设置左右端点的1阶导数曲率, d第1个和最后一个曲线设置为抛物线, e第1,2和倒数两个曲线分别设置为同一条曲线"""
    if kind == 'a':
        A[0][0] = 1.
        A[-1][-1] = 1.
    elif kind == 'b':
        v1 = v[0]
        v2 = v[1]
        A[0][0] = 2.
        A[-1][-1] = 2.
        delta_b[0] = v1
        delta_b[-1] = v2
    elif kind == 'c':
        v1 = v[0]
        vn = v[1]
        A[0][0] = 2 * delta_list[0]
        A[0][1] = delta_list[0]
        A[-1][-1] = 2 * delta_list[-1]
        A[-1][-2] = delta_list[-1]

        delta_b[0] = 3. * (Delta_list[0] / delta_list[0] - v1)
        delta_b[-1] = 3. * (vn - Delta_list[-1] / delta_list[-1])
    elif kind == 'd':
        A[0, 0:2] = np.array([1., -1.])
        A[-1][-2:] = np.array([-1, 1.])
    elif kind == 'e':
        A[0][0] = delta_list[1]
        A[0][1] = -(delta_list[0] + delta_list[1])
        A[0][2] = delta_list[0]
        A[n - 1][n - 1] = delta_list[-2]
        A[n - 1][n - 2] = -(delta_list[-2] + delta_list[-1])
        A[n - 1][n - 3] = delta_list[-1]
    """求解方程"""
    x0 = np.zeros(n)
    c = Sor(1.2, A, x0, delta_b, 20)
    """计算多项式系数"""
    d_list = np.zeros(n - 1)
    b_list = np.zeros(n - 1)
    for i in range(1, n):
        d_list[i - 1] = (c[i] - c[i - 1]) / (3. * delta_list[i - 1])
        b_list[i - 1] = (Delta_list[i - 1] / delta_list[i - 1]) - (delta_list[i - 1] / 3) * (
                2 * c[i - 1] + c[i])
    """求出每段方程"""
    S = []
    for i in range(n - 1):
        s = b[i] + b_list[i] * (x - a[i]) + c[i] * (x - a[i]) ** 2 + d_list[i] * (x - a[i]) ** 3
        S.append(s)
    """绘图"""
    return S


def case0():
    """绘制插值函数的图像, 某些年的人口数量。下面的取值（0.95,1.05）不过是为了满足5%的置信水平而已"""
    """没什么好看的，似乎在持续上升"""
    a = [1960, 1970, 1990, 2000]
    b = [3039585530, 3707475887, 5281653820, 6079603571]
    y = newton_divided_differences(a, b)
    min_x = 0.95*min(a) if min(a) > 0 else 1.05*min(a)
    max_x = 1.05 * max(a) if max(a) > 0 else 0.95 * max(a)
    x_list = np.linspace(min_x, max_x, 100)
    y_list = []
    for i in range(len(x_list)):
        y0 = y.evalf(subs={'x':x_list[i]})
        y_list.append(y0)

    plt.plot(x_list, y_list, '-.')
    y11 = str(y)
    plt.legend([y11], loc='upper center')
    plt.show()

def case1():
    """观察一下龙格现象"""
    """似乎插值点越多，龙格现象越明显"""
    f = 1. / (1. + 12 * x ** 2.)
    x1_list = np.linspace(-1., 1., 14)
    x2_list = np.linspace(-1., 1., 24)
    x0_list = np.linspace(-1., 1., 100)
    y1_list = []
    y2_list = []
    f_list = []
    for i in range(len(x1_list)):
        y1 = f.evalf(subs={'x': x1_list[i]})
        y1_list.append(y1)
    for i in range(len(x2_list)):
        y2 = f.evalf(subs={'x': x2_list[i]})
        y2_list.append(y2)
    for i in range(len(x0_list)):
        f0 = f.evalf(subs={'x': x0_list[i]})
        f_list.append(f0)

    y1 = newton_divided_differences(x1_list, y1_list)
    y2 = newton_divided_differences(x2_list, y2_list)

    y11_list = []
    y12_list = []
    x_list = x0_list.copy()
    for i in range(len(x_list)):
        y11 = y1.evalf(subs={'x': x_list[i]})
        y11_list.append(y11)
        y12 = y2.evalf(subs={'x': x_list[i]})
        y12_list.append(y12)

    plt.plot(x_list, y11_list, '-.')
    plt.plot(x_list, y12_list, '-.')
    plt.plot(x0_list, f_list, ':', 'c')
    plt.legend(['15', '25', 'exact'], loc='best')
    plt.show()

def case2():
    """比较一下lagrange插值和newton插值的速度"""
    """似乎牛顿插值更快"""
    a = [0, 2, 3]
    b = [1, 2, 4]
    t1_s = time.perf_counter()
    y1 = lagrange_interpolation(a, b)
    t1_e = time.perf_counter()
    t1 = t1_e - t1_s
    t2_s = time.perf_counter()
    y2 = lagrange_interpolation(a, b)
    t2_e = time.perf_counter()
    t2 = t2_e - t2_s
    print('t1=%.5f' % t1)
    print('t2=%.5f' % t2)

def case3():
    """观察一下均匀采样和切比雪夫采样的插值函数之间的差异，对于sinx"""
    """不太明显啊，n取20,30,40.直到40对于均匀采样有明显的龙格现象。还是改成比较它们的误差比较明显吧.确实是这样。"""
    a = [0., 10.]
    n = 4
    f = sy.sin(x)
    x1_list = np.linspace(a[0], a[1], n)
    y1_list = []
    for i in range(len(x1_list)):
        y1 = f.evalf(subs={'x':x1_list[i]})
        y1_list.append(y1)
    x2_list = chebyshev_interpolation(a, n)
    y2_list = []
    for i in range(len(x2_list)):
        y2 = f.evalf(subs={'x':x2_list[i]})
        y2_list.append(y2)

    f1 = newton_divided_differences(x1_list, y1_list)
    f2 = newton_divided_differences(x2_list, y2_list)

    x_list = np.linspace(a[0], a[1], 100)
    f1_list = []
    f2_list = []
    f_list = []
    for i in range(len(x_list)):
        f1i = f1.evalf(subs={'x':x_list[i]})
        f1_list.append(f1i)
        f2i = f2.evalf(subs={'x':x_list[i]})
        f2_list.append(f2i)
        fi = f.evalf(subs={'x':x_list[i]})
        f_list.append(fi)

    plt.plot(x_list, f1_list, ":")
    plt.plot(x_list, f2_list, '-.')
    plt.plot(x_list, f_list, "-")
    plt.legend(['uniform','chebyshev', 'exact'], loc='best')
    plt.show()

    e_f1 = []
    e_f2 = []
    for i in range(len(x_list)):
        ef1 = f_list[i] - f1_list[i]
        ef2 = f_list[i] - f2_list[i]
        e_f1.append(ef1)
        e_f2.append(ef2)

    plt.plot(x_list, e_f1)
    plt.plot(x_list, e_f2)
    plt.legend(['uniform', 'chebyshev'], loc='best')
    plt.show()

"""与三次样条相关"""
def case4():
    """这图像一点也不像抛物线啊"""
    a = [0, 1, 2]
    b = [3, -2, 1]
    c = CubicSplines(a, b)
    c.cal()
    c.parabolically_terminated()
    c.cal_plot()

def case45():
    """这图像一点也不像抛物线啊"""
    a = [0, 1, 2]
    b = [3, -2, 1]
    c = CubicSplines(a, b)
    c.cal()
    c.not_a_kont()
    c.cal_plot()

def case3_4_7():
    """使用clamp三次样条曲线绘制cos在0，pi/2上的图像"""
    """对于端点的一次导数，当然是选择和cos相同的值才能保证样条曲线的精确性啊.似乎没什么差别。"""
    x_list = np.linspace(0., 2., 50)
    m = len(x_list)
    y_list = np.zeros(m)
    for j in range(m):
        y_list[j] = np.cos(x_list[j])
    plt.plot(x_list, y_list, '-.')

    """S中将样条曲线左右端点取为了cos对应的一阶导数"""
    a = np.linspace(0., np.pi/2, 5)
    n = len(a)
    b = np.zeros(n)
    for i in range(n):
        b[i] = np.cos(a[i])
    S = cubic_splines(a, b, 'c', 0., -1.)
    """绘制三次样条曲线，这一段代码有必要重复使用"""
    """m是对每一段样条曲线上取m个点来离散的绘制样条曲线，散点图是画出每段样条曲线的端点。"""
    m = 10
    plt.scatter(a, b)
    for i in range(n - 2):
        y_list = np.zeros(m)
        x_list = np.linspace(a[i], a[i + 1], m)
        for j in range(len(x_list)):
            y_list[j] = S[i].evalf(subs={'x': x_list[j]})
        plt.plot(x_list, y_list, '-.c')
    x_last = np.linspace(a[-2], 2., m)
    y_last = np.zeros(m)
    for i in range(len(x_last)):
        y_last[i] = S[-1].evalf(subs={'x': x_last[i]})
    plt.plot(x_last, y_last, '-.c')
    plt.show()

def case3_4_8():
    """使用clamp三次样条曲线绘制cos在0，pi/2上的图像"""
    """对于端点的一次导数，当然是选择和cos相同的值才能保证样条曲线的精确性啊.似乎没什么差别。"""
    x_list = np.linspace(0., 2., 50)
    m = len(x_list)
    y_list = np.zeros(m)
    for j in range(m):
        y_list[j] = np.sin(x_list[j])
    plt.plot(x_list, y_list, '-.')

    """S中将样条曲线左右端点取为了cos对应的一阶导数"""
    a = np.linspace(0., np.pi/2, 5)
    n = len(a)
    b = np.zeros(n)
    for i in range(n):
        b[i] = np.sin(a[i])
    S = cubic_splines(a, b, 'c', 1., 0.)
    """绘制三次样条曲线，这一段代码有必要重复使用"""
    """m是对每一段样条曲线上取m个点来离散的绘制样条曲线，散点图是画出每段样条曲线的端点。"""
    m = 10
    plt.scatter(a, b)
    for i in range(n - 2):
        y_list = np.zeros(m)
        x_list = np.linspace(a[i], a[i + 1], m)
        for j in range(len(x_list)):
            y_list[j] = S[i].evalf(subs={'x': x_list[j]})
        plt.plot(x_list, y_list, '-.c')
    x_last = np.linspace(a[-2], 2., m)
    y_last = np.zeros(m)
    for i in range(len(x_last)):
        y_last[i] = S[-1].evalf(subs={'x': x_last[i]})
    plt.plot(x_last, y_last, '-.c')
    plt.show()
    
def case3_4_9():
    """使用clamp三次样条曲线绘制ln在1,3上的图像"""
    """它的导数是个递减函数，而误差对于同一函数不同的x应该是x越小，误差越大"""
    x_list = np.linspace(1., 3., 50)
    m = len(x_list)
    y_list = np.zeros(m)
    for j in range(m):
        y_list[j] = np.log(x_list[j])
    plt.plot(x_list, y_list, '-.')
    plt.show()
    """S中将样条曲线左右端点取为了ln对应的一阶导数"""
    a = np.linspace(1., 3., 15)
    n = len(a)
    b = np.zeros(n)
    for i in range(n):
        b[i] = np.log(a[i])
    S = cubic_splines(a, b, 'c', 1., 0.2)
    """绘制三次样条曲线，这一段代码有必要重复使用"""
    """m是对每一段样条曲线上取m个点来离散的绘制样条曲线，散点图是画出每段样条曲线的端点。"""
    m = 10
    plt.scatter(a, b)
    for i in range(n - 1):
        y_list = np.zeros(m)
        x_list = np.linspace(a[i], a[i + 1], m)
        for j in range(len(x_list)):
            y_list[j] = S[i].evalf(subs={'x': x_list[j]})
        plt.plot(x_list, y_list, '-.c')
    plt.show()

def case3_4_10():
    """看看在1,3上要去n=?个点可以使最大误差小于error"""
    """怎么才能确定最大误差？解决了这个问题就可以吧代码普遍化了"""
    error = 0.5*10**-7

    for n in range(5, 20):
        mul = 1.
        for j in range(n):
            mul *= 1+(2/n)*j
            mul /= j+1
        if mul <= error:
            print("i=%d" % n)
            print("mul=%.8f" % mul)
            break

def case3_4_11():
    """使用自然三次样条和clamp三次样条，clamp样条端点斜率采用源数据的两点斜率"""
    """clamp预测的误差更小一点"""
    a = [1960, 1970, 1990, 2000]
    n = len(a)
    b = [3039585530, 3707475885, 5281653820,  6079603571]
    plt.scatter(a, b)

    S1 = cubic_splines(a, b, 'a')
    S2 = cubic_splines(a, b, 'c', 66789035.7, 79794975.1)
    m = 8
    for i in range(n - 1):
        y_list = np.zeros(m)
        x_list = np.linspace(a[i], a[i + 1], m)
        for j in range(len(x_list)):
            y_list[j] = S1[i].evalf(subs={'x': x_list[j]})
        plt.plot(x_list, y_list, '-.c')
    plt.show()
    for i in range(n - 1):
        y_list = np.zeros(m)
        x_list = np.linspace(a[i], a[i + 1], m)
        for j in range(len(x_list)):
            y_list[j] = S2[i].evalf(subs={'x': x_list[j]})
        plt.plot(x_list, y_list, '-.')
    y1_1980 = S1[1].evalf(subs={'x': 1980})
    print('y1_1980=%.f' % (y1_1980-4452584592.))
    y2_1980 = S2[1].evalf(subs={'x': 1980})
    print('y2_1980=%.f' % (y2_1980 - 4452584592.))

def case5():
    """看看二维贝塞尔曲线"""
    """它的端点斜率是由端点和控制点决定的,当左端点和第一个控制点相同，右端点和第二个控制点相同时，贝塞尔曲线是直线"""
    a = [[1., 1.], [1., 3.], [3., 3.], [2., 2.]]
    bc = bezier_curves(a)
    xt = bc[0]
    yt = bc[1]

    t_list = np.linspace(0., 1., 50)
    x_list = []
    y_list = []
    for i in range(len(t_list)):
        xi = xt.evalf(subs={'t' : t_list[i]})
        yi = yt.evalf(subs={'t' : t_list[i]})
        x_list.append(xi)
        y_list.append(yi)

    plt.plot(x_list, y_list, '-')
    plt.plot([a[0][0], a[1][0]], [a[0][1], a[1][1]], '-o')
    plt.plot([a[2][0], a[3][0]], [a[2][1], a[3][1]], '-o', color='c')
    plt.show()
    pass

def case6():
    """绘制多条贝塞尔曲线， 实际当中，点的值应该就是这样输入的"""
    a0 = [[0., 2., 1., 2., 1., 1., 0., 1.],
          [0, 1, 1, 1, 1, 0, 0, 0]]
    t_list = np.linspace(0., 1., 30)
    x_list = []
    y_list = []
    bc = bezier_curves(a0)
    xt = bc[0]
    yt = bc[1]
    for i in range(len(a0)):
        for j in range(len(t_list)):
            xi = xt[i].evalf(subs={'t': t_list[j]})
            yi = yt[i].evalf(subs={'t': t_list[j]})
            x_list.append(xi)
            y_list.append(yi)
    plt.plot(x_list, y_list, '-')
    plt.show()

if __name__ == "__main__":
    fem_case1()

