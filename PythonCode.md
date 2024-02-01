# Tauc plots and calculation of band gap with uncertainty propagation

from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from scipy.stats import linregress
import math

def f(x, a, b):
    return a*x + b

#Open file with wavelength and absorbance
lamb, A = np.loadtxt(r'C:\Users\User\OneDrive\Desktop\Tauc.txt', delimiter='\t', unpack=True)

#Calculation of energy of incident photon and (F(R)E)**(1/2)
h = 6.62607004*(10**(-34))
c = 3.00*(10**(8))
j = 1.602177*(10**(-19))

x = ((h*c)/(lamb*(10**(-9))))/j
y = (x*A)**(1/2)

#Derivatives
dy = np.diff(y, 1)
dx = np.diff(x, 1)
dydx = dy / dx
maxindex = np.argmax(dydx)

#Band gap
x_linear = x[maxindex - 10: maxindex + 10]
y_linear = y[maxindex - 10: maxindex + 10]
a, b, r_value, p_value, stderr = linregress(x_linear, y_linear)
E_bandgap = round(-b/a, 3)

#Uncertainty 
n = 560 #number of points
delta_a = a*math.sqrt((((r_value)**(-2))-1)/(n-2))

def _sum(xi):
    sum=0
    for i in xi:
        sum = sum + i
    return(sum)

if __name__ == "__main__":
    xi = np.array(x)
    i = len(xi)
    ans = _sum(xi)

delta_b = delta_a*math.sqrt((1/n)*(ans))

delta_x = math.fabs(b/(a**2))*math.fabs(delta_a)+math.fabs(-1/a)*math.fabs(delta_b)

#Plots
visualization_x = np.linspace(E_bandgap, x[maxindex-60], 2)

plt.scatter(E_bandgap, 0, marker='x', color='k', label="Bandgap = " + str(E_bandgap) + "eV")
plt.plot(x, y)
plt.plot(dx, dy)
plt.plot(visualization_x, f(visualization_x, a, b), color='red')
plt.xlabel("Energia")
plt.ylabel("$(alpha h nu)^2$")
plt.title("Tauc plot")
plt.legend()
plt.grid()
plt.show()

#Values
print("a=",a,"b=", b,"r_value", r_value,"delta_a=", delta_a,"delta_b=", delta_b, "delta_x=", delta_x)
