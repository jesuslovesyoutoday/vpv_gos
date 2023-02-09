from sympy import integrate, symbols, diff, pi, evalf, log, cos
import matplotlib.pyplot as plt
import numpy as np
import math


def kompton(e, a):
    return e/(1 + e*(1-cos(a))/0.511e6)
    
x = np.linspace(-np.pi, np.pi, 1000)

a, e, b = symbols('a e b')
r = 2.81e-13
E = 0.6617e6
a1 = 8.13
a2 = 9.46
a3 = a2
#e = E/(0.511e6 + E*(1 - cos(alpha)))
N0 = 10e8

alpha = a2 * math.pi / 180
e = E/(0.511e6 + E*(1 - b))
sigma = 2*math.pi*r**2*((1+e)/(e**2) * ((2 + 2 *e)/(1 + 2*e) - log(1 + 2*e)/e) + log(1 + 2*e)/(2*e) - (1 + 3*e)/((1 + 2*e)**2))
dsigma = -diff(sigma, b)
sigma = sigma.evalf(subs = {b:cos(alpha)})
dsigma = dsigma.evalf(subs = {b:cos(alpha)})
N1 = N0 * dsigma * 3*64/(4*math.pi*r**2)
N = N1 * 1.47e-25 * 3*64/(4*math.pi*r**2)
E = kompton(E, alpha)
print('E = ', E)
print('DSigma = ', dsigma)

e = E/(0.511e6 + E*(1 - b))
sigma = 2*math.pi*r**2*((1+e)/(e**2) * ((2 + 2 *e)/(1 + 2*e) - log(1 + 2*e)/e) + log(1 + 2*e)/(2*e) - (1 + 3*e)/((1 + 2*e)**2))
dsigma = -diff(sigma, b)
sigma = sigma.evalf(subs = {b:cos(alpha)})
dsigma = dsigma.evalf(subs = {b:cos(alpha)})
N1 = N * dsigma * 3*64/(4*math.pi*r**2)
N2 = N1 * 1.47e-25 * 3*64/(4*math.pi*r**2)
E = kompton(E, alpha)

print('E = ', E)
print('N2 = ', N2)
print('Sigma = ', sigma)
print('DSigma = ', dsigma)

