from sympy import integrate, symbols, diff, cos, log, pi

a, e = symbols('a e')

e = 0.6617/(0.511 + 0.6617*(1 - cos(a)))
r = 2.81e-15

sigma = 2*pi*(r**2) * ((1+e)/(e**2) * ((2 + 2 *e)/(1 + 2*e) - log(1 + 2*e)/e) + log(1 + 2*e)/(2*e) - (1 + 3*e)/((1 + 2*e)**2))

dsigma = diff(sigma, a)
integr = integrate(dsigma, (a, 0, 2*pi))
print(integr)




