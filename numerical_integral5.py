import sympy as sy
import re

from interpolation3 import lagrange_interpolation, lagrange_interpolation_


class Integral():
    def __init__(self):
        pass
    
    def cal(self, var='x'):
        pass
    
    def cal_error(self, f_c):
        pass
    
    
x, y, z, t = sy.symbols("x, y, z, t")
f1 = sy.ln(x)
arange = [1., 2.]

"""
一个简单的牛顿科特斯计算算法，1阶线性。
输入：一个函数和积分的两个端点
输出：该函数在积分区间上的积分值
"""
def newton_cotes1(f, array):

    f_0 = f.evalf(subs={x:array[0]})
    f_1 = f.evalf(subs={x:array[1]})
    out = (array[1] - array[0]) * (0.5*f_0 + 0.5*f_1)
    return out
def newton_cotes2(f, array):
    x_1 = 0.5*(array[0] + array[1])
    f_0 = f.evalf(subs={x: array[0]})
    f_1 = f.evalf(subs={x: x_1})
    f_2 = f.evalf(subs={x: array[1]})
    out = (array[1] - array[0]) * (f_0 + 4.*f_1 + f_2) / 6.
    return out

def com_newton_cotes1(f,array,a):
    h = (array[1] - array[0]) / a
    sum_f_i = 0
    for i in range(1, a):
        x_i = array[0] + i * h
        sum_f_i += f.evalf(subs={x: x_i})

    f_0 = f.evalf(subs={x: array[0]})
    f_1 = f.evalf(subs={x: array[1]})
    out = 0.5 * h * (f_0 + 2.*sum_f_i + f_1)
    return out
def com_newton_cotes2(f,array,m):
    while m % 2 != 0:
        m = input("m = ")
    h = (array[1] - array[0]) / m / 2.
    sum_f_envy = 0
    sum_f_odd = 0.
    for i in range(1, m+1):
        x_odd = array[0] + (2 * i - 1) * h
        sum_f_odd += f.evalf(subs={x: x_odd})
    for i in range(1, m):
        x_envy = array[0] + (2 * i) * h
        sum_f_envy += f.evalf(subs={x: x_envy})

    f_0 = f.evalf(subs={x: array[0]})
    f_1 = f.evalf(subs={x: array[1]})
    out = h * (f_0 + 4. * sum_f_odd + 2. * sum_f_envy + f_1) / 3.
    return out

"""总觉得使用将函数和积分区间作为类的参数传入不太好。"""
class NewtonCotes1():
    """输入为函数和积分区间"""
    def __init__(self, f, _range):
        self.f = f
        self._range = _range
        self.h = self._range[1] - self._range[0]

    """当想使用复合积分的时候，可以指定需要划分的区间个数"""
    def cal(self, n=1):
        self.n = n
        if self.n == 1:
            f_0 = self.f.evalf(subs={x: self._range[0]})
            f_1 = self.f.evalf(subs={x: self._range[1]})
            self.output = 0.5 * self.h * (f_0 + f_1)                 
        else:
            self.h = (self._range[1] - self._range[0]) / n
            sum_f_i = 0
            for i in range(1, n):
                x_i = self._range[0] + i * self.h
                sum_f_i += self.f.evalf(subs={x: x_i})

            f_0 = self.f.evalf(subs={x: self._range[0]})
            f_1 = self.f.evalf(subs={x: self._range[1]})
            self.output = 0.5 * self.h * (f_0 + 2. * sum_f_i + f_1)          
        return self.output

    def cal_error(self, f_c=1.):
        if self.n == 1:
            self.error = self.h ** 3 * f_c / 12
        else:
            self.error = (self._range[1] - self._range[0]) * self.h ** 2 * f_c / 12
        return self.error


"""我想构造一个更加通用的积分器，排除上面那种里面指定了自变量是x，我想实现对于某函数里面任意自变量的积分"""
class NewtonCotes1_v1():
    def __init__(self, f, _range):
        self.f = f
        self._range = _range
        self.h = self._range[1] - self._range[0]

    """当想使用复合积分的时候，可以指定需要划分的区间个数"""

    def cal(self, var, n=1):
        var = sy.symbols(var)
        self.n = n
        if self.n == 1:
            f_0 = self.f.subs({var: self._range[0]})
            f_1 = self.f.subs({var: self._range[1]})
            self.output = 0.5 * self.h * (f_0 + f_1)
        else:
            self.h = (self._range[1] - self._range[0]) / n
            sum_f_i = 0
            for i in range(1, n):
                var_i = self._range[0] + i * self.h
                sum_f_i += self.f.evalf(subs={var: var_i})

            f_0 = self.f.evalf(subs={var: self._range[0]})
            f_1 = self.f.evalf(subs={var: self._range[1]})
            self.output = 0.5 * self.h * (f_0 + 2. * sum_f_i + f_1)
        return self.output

    """我是没办法自动计算这个fc的导数了， 还是让ueser自己算了填进来吧。要不给个默认值？还是再写个类计算它。"""
    def cal_error(self, f_c=1.):
        if self.n == 1:
            self.error = self.h ** 3 * f_c / 12
        else:
            self.error = (self._range[1] - self._range[0]) * self.h ** 2 * f_c / 12
        return self.error
    
    
class NewtonCotes2(Integral):
    def __init__(self, f, _range):
        super().__init__()
        self.f = f
        self._range = _range
        self.h = 0.5 * (self._range[1] - self._range[0])
        """对于这个n我想起来以前被它的耦合给折腾过。"""
    
    def cal(self, n=1):
        self.n = n
        if self.n == 1:
            x_1 = 0.5 * (self._range[0] + self._range[1])
            f_0 = self.f.evalf(subs={x: self._range[0]})
            f_1 = self.f.evalf(subs={x: x_1})
            f_2 = self.f.evalf(subs={x: self._range[1]})
            self.output = self.h * (f_0 + 4. * f_1 + f_2) / 3.
        else:
            while self.n % 2 != 0:
                self.n = input("n = ")
            h = (self._range[1] - self._range[0]) / self.n / 2.
            sum_f_even = 0
            sum_f_odd = 0.
            for i in range(1, self.n + 1):
                x_odd = self._range[0] + (2 * i - 1) * h
                sum_f_odd += self.f.evalf(subs={x: x_odd})
            for i in range(1, self.n):
                x_even = self._range[0] + (2 * i) * h
                sum_f_even += self.f.evalf(subs={x: x_even})

            f_0 = self.f.evalf(subs={x: self._range[0]})
            f_1 = self.f.evalf(subs={x: self._range[1]})
            self.output = h * (f_0 + 4. * sum_f_odd + 2. * sum_f_even + f_1) / 3.
        return self.output

    def cal_error(self, f_c=1.):
        if self.n == 1:
            self.error = self.h ** 5 * f_c / 90
        else:
            self.error = (self._range[1] - self._range[0]) * self.h ** 4 * f_c / 180
        return self.error


class NewtonCotes2_v1(Integral):
    def __init__(self, f, _range):
        super().__init__()
        self.f = f
        self._range = _range
        self.h = 0.5 * (self._range[1] - self._range[0])
        """对于这个n我想起来以前被它的耦合给折腾过。"""

    def cal(self, var='x', n=1):
        var = sy.symbols(var)
        self.n = n
        if self.n == 1:
            x_1 = 0.5 * (self._range[0] + self._range[1])
            f_0 = self.f.subs({var: self._range[0]})
            f_1 = self.f.subs({var: x_1})
            f_2 = self.f.subs({var: self._range[1]})
            self.output = self.h * (f_0 + 4. * f_1 + f_2) / 3.
        else:
            while self.n % 2 != 0:
                self.n = input("n = ")
            h = (self._range[1] - self._range[0]) / self.n / 2.
            sum_f_even = 0
            sum_f_odd = 0.
            for i in range(1, self.n + 1):
                x_odd = self._range[0] + (2 * i - 1) * h
                sum_f_odd += self.f.subs({var: x_odd})
            for i in range(1, self.n):
                x_even = self._range[0] + (2 * i) * h
                sum_f_even += self.f.subs({var: x_even})

            f_0 = self.f.subs({var: self._range[0]})
            f_1 = self.f.subs({var: self._range[1]})
            self.output = h * (f_0 + 4. * sum_f_odd + 2. * sum_f_even + f_1) / 3.
        return self.output

    def cal_error(self, f_c=1.):
        if self.n == 1:
            self.error = self.h ** 5 * f_c / 90
        else:
            self.error = (self._range[1] - self._range[0]) * self.h ** 4 * f_c / 180
        return self.error
    
"""应用了一个sy对x的求导"""
def gauss_x(n=3):
    f = (x**2 - 1) ** n
    f_list = []
    f_list.append(f)
    for i in range(n):
        f_list.append(sy.diff(f_list[i], x))
    n_ = 1
    for i in range(n):
        n_ *= (i + 1)
    f_ = f_list[-1] / (2**n) / n_
    x_list = sy.solve(f_)
    x_list = [float(i) for i in x_list]
    return x_list

class GaussianQuadrature():
    def __init__(self, f, _range):
        self.f = f
        self._range = _range

    """n指的是高斯点的个数。"""
    def cal(self, n=2):
        x_list = gauss_x(n)
        x_newlist = [0.5 * (self._range[1]-self._range[0]) * i + 0.5 * (self._range[0] + self._range[1]) for i in x_list]
        l_list = lagrange_interpolation_(x_list)

        alpha_list = []
        value = 0
        for i in range(len(x_list)):
            """这个alpha_list是不是应该这样调整我还不是很确定。这个积分法估计还得设置成能修改的。"""
            alpha_list.append(0.5 * (self._range[1]-self._range[0]) * NewtonCotes2(l_list[i], [-1, 1]).cal())
            value += alpha_list[i] * self.f.evalf(subs={x: x_newlist[i]})
        return value


class GaussianQuadrature_v1():
    def __init__(self, f, _range):
        self.f = f
        self._range = _range

    """n指的是高斯点的个数。"""
    def cal(self, var='x', n=2):
        var = sy.symbols(var)
        x_list = gauss_x(n)
        x_newlist = [0.5 * (self._range[1]-self._range[0]) * i + 0.5 * (self._range[0] + self._range[1]) for i in x_list]
        l_list = lagrange_interpolation_(x_list)

        alpha_list = []
        value = 0
        for i in range(len(x_list)):
            """这个alpha_list是不是应该这样调整我还不是很确定。这个积分法估计还得设置成能修改的。"""
            alpha_list.append(0.5 * (self._range[1]-self._range[0]) * NewtonCotes2(l_list[i], [-1, 1]).cal())
            value += alpha_list[i] * self.f.subs({var: x_newlist[i]})
        return value


"""我需要一个积分器，能够在二维、三维积分的时候选择不同的积分方式。
我想利用已有的类。"""
class Integral0(Integral):
    def __init__(self, f, _range):
        super().__init__()
        self.f = f
        self._range = _range

    """可能需要一个函数来检查输入的f是几维的。不过opensees里面维度是自己指定的。"""
    def identify_dimen(self):
        f_new = ''.join(re.split(r'[^A-Za-z]', str(self.f)))
        self.dimen = len("".join(set(f_new)))

    def assign_integral(self, *args):
        pass

    def cal(self):
        pass

    def cal_error(self, f_c):
        pass


def case0():
    f = sy.ln(x)
    array = [1., 2.]
    newtoncotes1 = NewtonCotes1(f, array)
    value = newtoncotes1.cal()
    error = newtoncotes1.cal_error(1)
    print('value=%.4f, error=%.4f' % (value, error))
def case1():
    f = x ** 2 * y ** 2
    newtoncotes2 = NewtonCotes2_v1(f, [-1, 1])
    value = newtoncotes2.cal('x')
    gauss = GaussianQuadrature_v1(value, [-1, 1])
    value1 = gauss.cal('y')
    print(value1)

if __name__ == '__main__':
    case1()
