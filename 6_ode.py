import sympy as sy
import numpy as np
import matplotlib.pyplot as plt

t, w = sy.symbols('t, w')
"""y' = ty + t**3 的精确解是y = 3exp(t**2 / 2) - t**2 -2"""
"""【0， 1】，y0=1"""

def euler(f, y0, inter, n):
    w_list = [0] * (n+1)
    t_list = [0] * (n+1)
    h = (inter[1] - inter[0]) / n
    w_list[0] = y0
    t_list[0] = inter[0]
    for i in range(n):
        diff_y = f.evalf(subs={'t': t_list[i], 'w': w_list[i]})
        w_list[i+1] = w_list[i] + h * diff_y
        t_list[i+1] = inter[0] + h * (i + 1)
    return t_list, w_list

def explicit_trapezoid(f, y0, inter, n):
    w_list = [0] * (n + 1)
    t_list = [0] * (n + 1)
    h = (inter[1] - inter[0]) / n
    w_list[0] = y0
    t_list[0] = inter[0]
    for i in range(n):
        diff_y1 = f.evalf(subs={'t': t_list[i], 'w': w_list[i]})
        diff_y2 = f.evalf(subs={'t': (t_list[i] + h), 'w': (w_list[i] + h * diff_y1)})
        w_list[i + 1] = w_list[i] + 0.5 * h * (diff_y1 + diff_y2)
        t_list[i + 1] = inter[0] + h * (i + 1)
    return t_list, w_list

def midpoint(f, y0, inter, n):
    w_list = [0] * (n + 1)
    t_list = [0] * (n + 1)
    h = (inter[1] - inter[0]) / n
    w_list[0] = y0
    t_list[0] = inter[0]

    for i in range(n):
        s_1 = f.evalf(subs={'t': t_list[i], 'w': w_list[i]})
        s_2 = f.evalf(subs={'t': (t_list[i] + 0.5 * h), 'w': (w_list[i] + 0.5 * h * s_1)})
        w_list[i + 1] = w_list[i] + h * s_2
        t_list[i + 1] = inter[0] + h * (i + 1)
    return t_list, w_list

def RK2(f, y0, inter, n):
    w_list = [0] * (n + 1)
    t_list = [0] * (n + 1)
    h = (inter[1] - inter[0]) / n
    w_list[0] = y0
    t_list[0] = inter[0]

    for i in range(n):
        s_1 = f.evalf(subs={'t': t_list[i], 'w': w_list[i]})
        s_2 = f.evalf(subs={'t': (t_list[i] + h), 'w': (w_list[i] + h * s_1)})
        w_list[i + 1] = w_list[i] + h * (s_1 + s_2) / 2.
        t_list[i + 1] = inter[0] + h * (i + 1)
    return t_list, w_list

def RK3(f, y0, inter, n):
    w_list = [0] * (n + 1)
    t_list = [0] * (n + 1)
    h = (inter[1] - inter[0]) / n
    w_list[0] = y0
    t_list[0] = inter[0]

    for i in range(n):
        s_1 = f.evalf(subs={'t': t_list[i], 'w': w_list[i]})
        s_2 = f.evalf(subs={'t': (t_list[i] + h), 'w': (w_list[i] + h * s_1)})
        s_3 = f.evalf(subs={'t': (t_list[i] + 0.5 * h), 'w': (w_list[i] + 0.5 * h * 0.5 * (s_1 + s_2))})
        w_list[i + 1] = w_list[i] + h * (s_1 + s_2 + 4. * s_3) / 6.
        t_list[i + 1] = inter[0] + h * (i + 1)
    return t_list, w_list

def RK4(f, y0, inter, n):
    w_list = [0] * (n + 1)
    t_list = [0] * (n + 1)
    h = (inter[1] - inter[0]) / n
    w_list[0] = y0
    t_list[0] = inter[0]

    for i in range(n):
        s_1 = f.evalf(subs={'t': t_list[i], 'w': w_list[i]})
        s_2 = f.evalf(subs={'t': (t_list[i] + 0.5 * h), 'w': (w_list[i] + 0.5 * h * s_1)})
        s_3 = f.evalf(subs={'t': (t_list[i] + 0.5 * h), 'w': (w_list[i] + 0.5 * h * s_2)})
        s_4 = f.evalf(subs={'t': t_list[i] + h, 'w': w_list[i] + h * s_3})
        w_list[i + 1] = w_list[i] + h * (s_1 + 2. * s_2 + 2. * s_3 + s_4) / 6.
        t_list[i + 1] = inter[0] + h * (i + 1)
    return t_list, w_list

def ode23(f, y0, inter, h, T):
    t_list = []
    w2_list = []
    w3_list = []
    t_list.append(inter[0])
    w2_list.append(y0)
    w3_list.append(y0)
    i = 0

    while t_list[-1] < inter[1]:
        s_1 = f.evalf(subs={'t': t_list[i], 'w': w2_list[i]})
        s_2 = f.evalf(subs={'t': (t_list[i] + h), 'w': (w2_list[i] + h * s_1)})
        s_3 = f.evalf(subs={'t': (t_list[i] + 0.5 * h), 'w': (w2_list[i] + 0.5 * h * 0.5 * (s_1 + s_2))})

        w2_list.append(w2_list[i] + h * (s_1 + s_2) / 2.)
        w3_list.append(w3_list[i] + h * (s_1 + s_2 + 4. * s_3) / 6.)

        error = abs(w3_list[i+1] - w2_list[i+1]) / abs(w2_list[i+1])
        error1 = abs(w3_list[i+1] - w2_list[i+1]) / max(abs(w2_list[i+1]), 1.)
        """主要的区别就在这里啊"""
        if error1 <= T:
            t_list.append(t_list[i] + h)
            h *= 0.8 * (T / error) ** (1 / (2 + 1))
            i += 1
        else:
            h *= 0.5
            w2_list.pop()
            w3_list.pop()
            """这里把t超出右端点时候换成右端点"""
    h = inter[1] - t_list[-2]
    s_1 = f.evalf(subs={'t': t_list[i-1], 'w': w2_list[i-1]})
    s_2 = f.evalf(subs={'t': (t_list[i-1] + h), 'w': (w2_list[i-1] + h * s_1)})
    w2_list[-1] = w2_list[i-1] + h * (s_1 + s_2) / 2.
    t_list[-1] = inter[1]
    return t_list, w2_list

def ode45(f, y0, inter, h, T):
    t_list = []
    w4_list = []
    w5_list = []
    t_list.append(inter[0])
    w4_list.append(y0)
    w5_list.append(y0)
    i = 0

    while t_list[-1] < inter[1]:
        s_1 = f.evalf(subs={'t': t_list[i], 'w': w4_list[i]})
        s_2 = f.evalf(subs={'t': (t_list[i] + 0.25*h), 'w': (w4_list[i] + 0.25*h*s_1)})
        s_3 = f.evalf(subs={'t': (t_list[i] + 0.375*h), 'w': (w4_list[i] + h * ((3/32)*s_1 + (9/32)*s_2))})
        s_4 = f.evalf(subs={
            't': (t_list[i] + (12/13) * h),
            'w': (w4_list[i] + h * ((1932/2197) * s_1 + (-7200/2197) * s_2 + (7296/2197)*s_3))})
        s_5 = f.evalf(subs={
            't': (t_list[i] + h),
            'w': (w4_list[i] + h * ((439/216) * s_1 + (-8) * s_2 + (3680/513)*s_3 + (-845/4104) * s_4))})
        s_6 = f.evalf(subs={
            't': (t_list[i] + 0.5 * h),
            'w': (w4_list[i] + h * ((-8/27) * s_1 + 2. * s_2 + (-3544/2565)*s_3 + (1859/4104)*s_4 + (-11/40)*s_5))})

        w4_list.append(w4_list[i] + h * ((25/216)*s_1 + (1408/2565)*s_3 + (2197/4104)*s_4 + (-1/5)*s_5))
        w5_list.append(w5_list[i] + h * ((16/135)*s_1 + (6656/12825)*s_3 + (28561/56430)*s_4 + (-9/50)*s_5 + (2/55)*s_6))

        error = abs(w5_list[i+1] - w4_list[i+1]) / abs(w4_list[i+1])

        if error <= T:
            t_list.append(t_list[i] + h)
            h *= 0.8 * (T / error) ** (1 / (4 + 1))

            i += 1
            w4_list[i] = w5_list[i]

        else:
            h *= 0.5
            w4_list.pop()
            w5_list.pop()
            """这里把t超出右端点时候换成右端点"""
    return t_list, w4_list

def test_0():
    """绘制同一常微分方程的精确解与步长不同（0.1和0.2）时解的精确度对比"""
    """可以看到随着步长的缩短，数值解距离精确解越近。"""
    euler_1 = euler(t * w + t ** 3., 1.0, [0., 1.], 5)
    plt.plot(euler_1[0], euler_1[1], '-o')

    euler_ = euler(t * w + t ** 3., 1.0, [0., 1.], 10)
    plt.plot(euler_[0], euler_[1], '->')

    euler_2 = euler(t * w + t ** 3., 1.0, [0., 1.], 40)
    plt.plot(euler_2[0], euler_2[1], '-.')

    f = 3. * sy.exp(0.5 * t ** 2.) - t ** 2. - 2.
    y_list = []
    for i in euler_[0]:
        y = f.evalf(subs={'t': i})
        y_list.append(y)
    plt.plot(euler_[0], y_list, '-<')
    plt.legend(['0.2', '0.1', '0.05', 'exact'], loc='best')
    plt.show()

def test_1():
    """在斜率场中绘制同一常微分方程在不同初值下的路径"""
    """可以看到当初值不同的时候，y的路径也是不同的。他们从初值出发，沿着斜率场中箭头的方向移动到右端点。"""
    x_list = np.linspace(0., 4., 30)
    y_list = np.linspace(-1., 2., 30)
    x, y = np.meshgrid(x_list, y_list)
    y_copy = []
    for i in range(len(y_list)):
        y_list[i] = 1. * y_list[i] * (1. - y_list[i])
        y_copy.append(y_list[i])
    y_copy1 = np.array(y_copy)

    u_list = np.ones((30, 30))

    v_list = np.tile(y_copy1, (30, 1))
    v_list = np.transpose(v_list)

    euler_ = euler(w * (1. - w), 1.5, [0., 4.], 30)
    plt.plot(euler_[0], euler_[1], '->')
    euler_1 = euler(w * (1. - w), 0.2, [0., 4.], 30)
    plt.plot(euler_1[0], euler_1[1], '-o')
    plt.quiver(x, y, u_list, v_list, color='c', width=0.003)
    plt.show()

def test_2():
    """比较中点方法和欧拉方法对于某无辜问题的精确程度"""
    """通过该例子的对比发现，欧拉方法在一些所谓无辜例子中无法得到很好的解答，但是中点方法却可以。"""
    et = explicit_trapezoid(-4.*t**3.*w**2., 1/10001, [-10., 0.], 1000)
    plt.plot(et[0], et[1], '-')
    euler_1 = euler(-4.*t**3.*w**2., 1/10001, [-10., 0.], 1000)
    plt.plot(euler_1[0], euler_1[1], '-.')

    plt.show()

def test_3(f):
    """对比一下两个二阶方法的精度,再加上四阶的方法"""
    """不愧都是二阶方法，函数图像基本一致，这精度几乎分辨不出来啊"""
    et = explicit_trapezoid(f, 1., [0., 2.], 20)
    plt.plot(et[0], et[1], '->')
    md = midpoint(f, 1., [0., 2.], 20)
    plt.plot(md[0], md[1], '-o')
    rk4 = RK4(f, 1., [0., 2.], 20)
    plt.plot(rk4[0], rk4[1], '-<')
    plt.show()

def test_4(f, y, y0, inter):
    """对比了1,2,4阶方法的步长-误差log-log图像"""
    """四阶方法的斜率是4,说明随着步长减半，误差变为2e-4.1,2阶对应斜率1,2"""
    n_list0 = np.logspace(0., 7., num=8, base=2)
    n_list = [int(5 * i) for i in n_list0]
    exact = y.evalf(subs={'t': 1.})

    h_list = []
    es_error_list = []
    mp_error_list = []
    rk4_error_list = []
    for n in n_list:
        h = (inter[1] - inter[0]) / n
        h_list.append(h)

        es = euler(f, y0, inter, n)
        es_error = float(abs(exact - es[1][-1]))
        es_error_list.append(es_error)

        mp = explicit_trapezoid(f, y0, inter, n)
        mp_error = float(abs(exact - mp[1][-1]))
        mp_error_list.append(mp_error)
        
        rk4 = RK4(f, y0, inter, n)
        rk4_error = float(abs(exact - rk4[1][-1]))
        rk4_error_list.append(rk4_error)

    lgt = range(len(h_list))
    h_list1 = [np.log10(h_list[i]) for i in lgt]
    es_error_list1 = [np.log10(es_error_list[i]) for i in lgt]
    mp_error_list1 = [np.log10(mp_error_list[i]) for i in lgt]
    rk4_error_list1 = [np.log10(rk4_error_list[i]) for i in lgt]

    plt.plot(h_list1, es_error_list1, '-o')
    plt.plot(h_list1, mp_error_list1, '-<')
    plt.plot(h_list1, rk4_error_list1, '->')
    plt.show()

def test_5():
    """比较ode23和精确解"""
    """搞的挺好的"""
    od23 = ode23(t*w + t**3, 1., [0., 1.], 0.5, 0.001)
    plt.plot(od23[0], od23[1], '-<')

    f = 3. * sy.exp(0.5 * t ** 2.) - t ** 2. - 2.
    y_list = []
    for i in od23[0]:
        y = f.evalf(subs={'t': i})
        y_list.append(y)
    plt.plot(od23[0], y_list, '-o')

    plt.show()

def test_6():
    """ode45的效果"""
    """容差调成1e-7这么快就算出来了"""
    od45= ode45(t * w + t ** 3, 1., [0., 1.], 0.5, 0.0000001)
    plt.plot(od45[0], od45[1], '-<')

    f = 3. * sy.exp(0.5 * t ** 2.) - t ** 2. - 2.
    y_list = []
    for i in od45[0]:
        y = f.evalf(subs={'t': i})
        y_list.append(y)
    plt.plot(od45[0], y_list, '-o')

    plt.show()

def test_7():
    """y'=10(1-y), 1. - 0.5 * sy.exp(-10 * t), y0=0.5, [0,100],ode45和精确解比较"""
    """这上下浮动也太恐怖了吧"""
    od45 = ode45(10.*(1 - w), 0.5, [0., 20.], 0.5, 0.0001)
    plt.plot(od45[0][10:], od45[1][10:], '-<')

    f = 1. - 0.5 * sy.exp(-10 * t)
    y_list = []
    for i in od45[0]:
        y = f.evalf(subs={'t': i})
        y_list.append(y)
    plt.plot(od45[0][10:], y_list[10:], '-o')

    plt.show()


if __name__ == "__main__":
    test_5()

