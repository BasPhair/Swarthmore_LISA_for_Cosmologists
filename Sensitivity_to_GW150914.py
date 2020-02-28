"""
Sensitivity to GW150914
"""
from scipy import integrate
from LISA_Object import LISA
import math
import matplotlib.pyplot as plt
import numpy as np
lisa = LISA('LISA_Data.csv')


#Definitions, units
GN = 6.67*(10**-11) # m^3/(kg s^2)
c = 3*(10**8) #m/s
Msol = 2*(10**30) #kg
Mpc = 3.086*(10**22) #m

#GW150914 paramaters
M1 = 36*Msol
M2 = 29*Msol
Z = .09

#chirp mass in source frame
Mr = (((M1*M2)**3)/(M1+M2))**.2

#Distance
integrand = lambda x: 2997.9/.7/math.sqrt((.3*(1+x)**3)+.7)
Dcom = Mpc*integrate.quad(integrand,0,Z)[0]
DL = (1+Z)*Dcom

#Initial Radiating Frequency, in observer frame
fi = .01774/(1+Z) #Radiates for 4 years from this frequency, in observer frame

#Time to radiate
fr = fi*(1+Z)
Tc = lambda frad: ((5*(c**5))/(256*(fr**(8/3))*((GN*Mr)**(5/3))*(math.pi()**(8/3)))) - ((5*(c**5))/(256*(frad**(8/3))*((GN*Mr)**(5/3))*(math.pi**(8/3))))

#... in observer frame frame (eq 98)
#check that it radiates for 4 years
print(((5*(c**5))/(256*(fr**(8/3))*((GN*Mr)**(5/3))*(math.pi**(8/3))))/(lisa.Year)*(1+Z))

#characteristic strain (eq 99)
hc = ((math.sqrt(2/3)*(((GN*Mr*(1+Z))/(c**3))**(5/6)))*(fi**(-1/6)))/((math.pi**(2/3))*DL/c)
hcf = lambda f: hc*((f/fi)**(-1/6))
integrand1 = lambda f: ((hcf(f)/(2*f))**2)/lisa.SumHApprox(f)
integrand2 = lambda f: ((hcf(f)/(2*f))**2)/lisa.SumHNum(f)
SNR1 = math.sqrt(integrate.quad(integrand1, fi, 1)[0])
print(SNR1)
SNR2 = math.sqrt(integrate.quad(integrand2, fi, 1)[0])
print(SNR2)

x = np.logspace(-5,0,1000)
x2 = np.logspace(math.log(fi,10),math.log(.05,10),200)
plt.loglog(x,[math.sqrt(4*i*lisa.SumHApprox(i)) for i in x], x2,[hcf(i) for i in x2],'--r')
plt.title('Sensitivity of LISA to GW150914')
plt.xlabel('Frequency (Hz)')
plt.ylabel('hc[strain]')
plt.show()

