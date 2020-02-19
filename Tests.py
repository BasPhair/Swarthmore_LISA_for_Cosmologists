"""
The following code is designed to test the results of the "Sensitivity to a
flat stochastic background" section of the code. The tests will compare the
results of the python calculations to the Mathematica functions () and
results derived from "the paper" (https://arxiv.org/pdf/1908.00546.pdf).
"""
from LISA_Object import LISA
import math
from scipy import integrate

test = LISA('LISA_Data.csv')



"""
The following code performs the calculations found in the "Sensitivity to a flat
stochastic background" section of the Mathematica code.
"""
print("The following code performs the calculations found in the 'Sensitivity to a flat\nstochastic background' section of the Mathematica code.")
fmin = 10**(-5)
fmax = 10**(0)

integrandApp = lambda x: (1/(test.SumOmegaApprox(x)**2))
integrandNum = lambda x: (1/(test.SumOmegaNum(x)**2))
def test_flat_background(fmin,fmax,integrand,SNR,T):
    return (SNR/math.sqrt(T*integrate.quad(integrand,fmin,fmax)[0]))


fmin = 10**(-5)
fmax = 10**(0)

print(test_flat_background(fmin,fmax,integrandApp,test.SNR5,test.T))
print('Actual: 4.79517x10^-13\n')

print(test_flat_background(fmin,fmax,integrandNum,test.SNR5,test.T))
print('Actual: 4.719x10^-13\n')

fmin = 10**(-4)
fmax = 10**(-1)

print(test_flat_background(fmin,fmax,integrandApp,test.SNR5,test.T))
print('Actual: 4.79517x10^-13\n')

print(test_flat_background(fmin,fmax,integrandNum,test.SNR5,test.T))
print('Actual: 4.719x10^-13\n')

print(test_flat_background(fmin,fmax,integrandNum,test.SNR5,test.T)*(.67**2))
print('Actual: 2.11836x10^-13\n')

print(test_flat_background(fmin,fmax,integrandNum,test.SNR5,(test.Year*3)))
print('Actual: 5.44903x10^-13\n')

print(test_flat_background(fmin,fmax,integrandNum,test.SNR5,(test.Year*3))*(.67**2))
print('Actual: 2.44607x10^-13\n')

print("Testing Amin(0)")
print(test.Amin(0,fmin,fmax))
print('Actual: 4.719x10^-13\n')

