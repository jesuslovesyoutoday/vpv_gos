from sympy import integrate, symbols, diff, cos, log, pi, solve, simplify
import matplotlib.pyplot as plt
import numpy as np
import math

def kompton(e, a):
    return e/(1 + e*(1-math.cos(a))/0.511)

#------------------------------------------------#

x = np.linspace(-np.pi, np.pi, 100)

a, e, b = symbols('a e b')
e = 0.6617e6/(0.511e6 + 0.6617e6*(1 - cos(a)))
r = 2.81e-15

E = 0.6617e6

angles = [0]
energy = [0.6617]
k = 0

EEE = [6e6, 0.6e6, 0.6e5, 0.6e4, 0.6]
#EEE = [6e6, 0.6e6]
Y_ = []
X_ = []

#------------------------------------------------#

for E in EEE:

    e = E/(0.511e6 + E*(1 - b))
    sigma = 2*pi*(r**2) * ((1+e)/(e**2) * ((2 + 2 *e)/(1 + 2*e) - log(1 + 2*e)/e) + log(1 + 2*e)/(2*e) - (1 + 3*e)/((1 + 2*e)**2))
    dsigma = -diff(sigma, b)
    dsigma = dsigma.subs({b:cos(a)})
    #ddsigma = diff(dsigma, a)
    
    y = []
    for i in x:
        #y.append(ddsigma.subs({a:i})/(2*pi))
        y.append(dsigma.evalf(subs ={a:i})/(2*pi))
        
    X = y * np.cos(x)
    Y = y * np.sin(x)
    X_.append(X)
    Y_.append(Y)
    
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
plt.grid()

for i in range(len(EEE)):
    x = X_[i]
    y = Y_[i]
    x = x/max(x)
    y = y/max(y)
    ax.plot(x, y, label = str(EEE[i]) + ' Ev')

ax.legend()
plt.show()
