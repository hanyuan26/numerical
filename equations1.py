import sympy as sy
import matplotlib.pyplot as plt
import numpy as np


x = sy.symbols('x')

def bisection_method(f, a, b, tol):
    tol0 = 0.5*10**-tol
    fa = f.evalf(subs={'x': a})
    fb = f.evalf(subs={'x': b})
    i = 0
    if fa * fb >= 0:
        print('ab之间没有根')
        return 0.
    else:
        while (b-a)/2. > tol0:
            fa = f.evalf(subs={'x': a})
            fb = f.evalf(subs={'x': b})
            c = (a+b)/2.
            fc = f.evalf(subs={'x': c})
            if fc == 0.:
                break
            if fa*fc < 0.:
                b = c
            elif fb*fc < 0.:
                a = c
            i += 1
        print('bm(i)=%d' % i)
        xc = (a+b)/2.
    return xc

def FPI(f, x0, tol):
    x_list = [x0]
    f_list = []
    f0 = f.evalf(subs={'x': x_list[0]})
    x_list.append(f0)
    f_list.append(f0)
    while abs(x_list[-1]-x_list[-2]) > tol:
        fi = f.evalf(subs={'x': x_list[-1]})
        f_list.append(fi)
        x_list.append(fi)
    return x_list, f_list

def FPI0(f, x0, k):
    """想画出那种几何图形，用容差的话它不收敛就没法完成计算。所以还是用迭代步数来控制一下"""
    x_list = [0] * (k+1)
    f_list = [0] * (k+1)
    x_list[0] = x0
    for i in range(k):
        f_list[i] = f.evalf(subs={'x': x_list[i]})
        x_list[i+1] = f_list[i]
    return x_list, f_list

def secant_method(f, x0, x1, tol):
    x_list = []
    x_list.append(x0)
    x_list.append(x1)
    while abs(x_list[-1] - x_list[-2]) >= tol:
        f0 = f.evalf(subs={'x': x_list[-2]})
        f1 = f.evalf(subs={'x': x_list[-1]})
        x = x_list[-1] - (f1*(x_list[-1] - x_list[-2])) / (f1 - f0)
        x_list.append(x)
    return x_list

def false_position(f, a, b, tol):
    i = 0

    while abs(b-a)/2. > tol:

        fa = f.evalf(subs={'x': a})
        fb = f.evalf(subs={'x': b})
        if fa * fb >= 0:
            print('ab之间没有根')
            break
        else:
            c = a - (fa*(a - b)) / (fa - fb)
            fc = f.evalf(subs={'x': c})

            if fc == 0.:
                break
            if fa*fc < 0.:
                b = c
            else:
                a = c
            i += 1
    print('sp(i)=%d' % i)
    xc = (a+b)/2.
    return xc

def test0():
    """使用二分法来求A的n次根"""
    """有时需要根据A中的值调整右端点的值，这时迭代步数会增加，不过幅度并不大。迭代次数似乎仅受到区间大小的影响"""
    """对于x^4-2，我把左右端点分别换成-10,10，采用这个形式的算法竟然说，它们中间不存在根。笑死。
    好吧不是这样用的，应该根据函数图像确定一个较小的范围，在此范围内寻找根.可是电脑又不知道根在哪"""
    A = [2, 3, 50]
    n = [2., 3., 4.]
    for j in n:
        for i in A:
            xc = bisection_method(x**j-i, 1., 10., 8)
            print('%.f的%.f次根xc = %.8f' % (i, j, xc))

def test1():
    """使用不动点迭代计算出x+sy.cos(x)-sy.sin(x), x0=0,r = np.pi / 4.
    每一步的x,gx，绝对误差e，ei+1/ei"""
    """可以看到在后期，ei+1/ei保持一个常值，这个值等于abs(g'(r))<1，说明误差在稳定的下降.这就是所谓的线性收敛
    收敛速度就是该常数值"""
    fpi = FPI0(x+sy.cos(x)-sy.sin(x), 0., 20)
    r = np.pi / 4.

    error = []
    s_error = []
    for i in range(len(fpi[0])):
        error.append(abs(r - fpi[0][i]))
    for i in range(len(fpi[0])-1):
        s = error[i+1]/error[i]
        s_error.append(s)
    s_error.insert(0, 0.1)
    print('xi\tgi\terror\ts_error')
    for i in range(len(fpi[0])):
        print('%.7f\t%.7f\t%.8f\t%.3f' % (fpi[0][i], fpi[1][i], error[i], s_error[i]))
    l = range(len(error))
    plt.plot(l, error, '-o')
    plt.plot(l, s_error)
    plt.legend(['error', 'ei+1/ei'], loc='best')
    plt.show()

def test11():
    """不动点迭代是怎么收敛的结合几何图形"""
    """这个方法不太行，fpi设置的是容差，如果不收敛就画不出图来，所以换成了下面那个test12，迭代多少步就画多少"""

    f = -1.5*x+2.5
    fpi = FPI(f, 0.8, 0.01)

    xs = np.linspace(0., 2., 100)
    fx_list = []
    for i in range(len(xs)):
        fx = f.evalf(subs={'x': xs[i]})
        fx_list.append(fx)

    x_list = fpi[0].copy()
    x_list.pop()
    g_list = fpi[1]
    xs_list = x_list.copy()
    gs_list = g_list.copy()
    for i in range(len(x_list)):
        xs_list.insert(2*i+1, x_list[i])
        gs_list.insert(2*i+1, g_list[i])
    gs_list.insert(0, 0)
    gs_list.pop()

    plt.plot(xs_list, gs_list, '->')
    plt.plot(xs, xs, '-')
    plt.plot(xs, fx_list, '-')
    plt.show()

def test12():
    """不动点迭代是怎么收敛的结合图形，-1.5*x+2.5，-0.5*x+1.5，x0=1，前者发散后者收敛。
    (1+2*x**3)/(1+3*x**2),x0=0.5"""
    """这个样子搞看不出来起点在哪啊，还有这个箭头的方向真不好.斜率小于1的话会收敛，越小越容易收敛"""
    f = -1.5*x+2.5
    fpi = FPI0(f, 0.8, 20)

    x_list = fpi[0].copy()
    x_list.pop()
    g_list = fpi[1]
    g_list.pop()
    xs_list = x_list.copy()
    gs_list = g_list.copy()
    for i in range(len(x_list)):
        xs_list.insert(2*i+1, x_list[i])
        gs_list.insert(2*i+1, g_list[i])
    gs_list.insert(0, 0)
    gs_list.pop()

    min_x = float(min(xs_list))
    max_x = float(max(xs_list))
    xs = np.linspace(min_x, max_x, 100)
    fx_list = []
    for i in range(len(xs)):
        fx = f.evalf(subs={'x': xs[i]})
        fx_list.append(fx)

    plt.plot(xs_list, gs_list, '->')
    plt.plot(xs, xs, '-')
    plt.plot(xs, fx_list, '-')
    plt.show()

def test2():
    """割线方法的几何过程。f = x**3.+x-1. x0=0., x1=1."""
    """最先的两个点形成割线与y=0形成交点，就是新的xi，该x在原函数上取得另一个点，与前一个点再形成割线计算下一个点"""
    f = sy.exp(x) + sy.sin(x) - 4.
    sm = secant_method(f, 1., 2., 1.*10**-10)
    f_list = []
    for i in range(len(sm)):
        f0 = f.evalf(subs={'x': sm[i]})
        f_list.append(f0)

    x_list = np.linspace(min(sm), max(sm), 50)
    y_list = []
    for i in range(len(x_list)):
        y0 = f.evalf(subs={'x': x_list[i]})
        y_list.append(y0)

    print(sm[-1])
    plt.plot(sm[:4], f_list[:4], '-o')
    plt.plot(x_list, y_list)
    plt.plot([min(sm), max(sm)], [0., 0.])
    plt.show()


def test3():
    """不动点迭代"""
    fp = false_position(x**3.+x-1., 0., 1., 10**-2)
    print(fp)

def test13():
    """试试割线方法收敛的速度"""
    """割线方法是介于线性和牛顿方法（二次）之间的方法，可以称作alpha = (1. + 5**0.5) / 2.阶
    但是这个东西不像fpi那样稳定，也许是因为没有精确解"""
    f = x**3 + x - 1.
    sm = secant_method(f, 0., 1., 10**-8)
    r = 0.682327803828
    alpha = (1. + 5**0.5) / 2.

    error = []
    s_error = []
    for i in range(len(sm)):
        error.append(abs(r - sm[i]))
    for i in range(len(sm) - 1):
        s = error[i + 1] / error[i] ** alpha
        s_error.append(s)
    s_error.insert(0, 0.1)
    print('xi\terror\ts_error')
    for i in range(len(sm)):
        print('%.7f\t%.8f\t%.3f' % (sm[i], error[i], s_error[i]))

    l = range(len(error))
    plt.plot(l, error, '-o')
    plt.plot(l, s_error)
    plt.legend(['error', '(ei+1/ei)**alpha'], loc='best')
    plt.show()

def test4():
    """比较一下几种算法速度"""
    f = sy.exp(x) + sy.sin(x) - 4.
    f1 = sy.exp(x) + x - 7.

    bm = bisection_method(f1, 1., 2., 5)
    fpi = FPI(f1, 1., 10**-5)
    sm = secant_method(f1, 1., 2., 10**-5)
    fp = false_position(f1, 1., 2., 10**-3)

    print('i(fpi)= %d' % len(fpi[0]))
    print('i(sm)= %d' % len(sm))

if __name__ == "__main__":
    test4()
