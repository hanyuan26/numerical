"""没法确定定义域，所以没法在两点和三点公式中选择"""
import sympy as sy
import matplotlib.pyplot as plt
import numpy as np


x = sy.symbols('x')

class Differentiation():
    f = ""
    value = 0.
    error = 0.
    def __init__(self, f, x):
        pass
    def calcualte(self):
        pass
    def choose_h(self):

        pass
    def geterror(self):
        pass

class Diff1(Differentiation):
    def __init__(self, f, x):
        self.f = f
        self.x = x
    def calcualte(self):
        value = (self.f.evalf(subs={'x': self.x + h}) - self.f.evalf(subs={'x': self.x - h})) / (2. * h)
        return value

def two_point_forward_df(f, x, h):
    value = (f.evalf(subs={'x': (x+h)}) - f.evalf(subs={'x': x})) / h
    return value

def three_point_center_df(f, x, h):
    value = (f.evalf(subs={'x': x+h}) - f.evalf(subs={'x': x-h})) / (2. * h)
    return value

def three_point_center_df2(f, x, h):
    value = (f.evalf(subs={'x': x + h}) + f.evalf(subs={'x': x - h}) - 2. * f.evalf(subs={'x': x})) / (2. * h)
    return value

def diff3(f, x, h):
    value = (-f.evalf(subs={'x': x-2*h}) + 2. * f.evalf(subs={'x': x - h}) - 2.*f.evalf(subs={'x': x + h}) + f.evalf(subs={'x': x+2.*h})) / (2.*h**3.)
    return value

def diff4(f, x, h):
    value = (f.evalf(subs={'x': x-2*h}) - 4. * f.evalf(subs={'x': x - h})
             + 6. * f.evalf(subs={'x': x})
             - 4.*f.evalf(subs={'x': x + h}) + f.evalf(subs={'x': x+2.*h})) / h**4.
    return value

def choose_h(f, x):
    diff_3 = diff3(f, x, 0.05)
    h = (3. * 1.e-15 / diff_3) ** (1 / 3)
    return h








