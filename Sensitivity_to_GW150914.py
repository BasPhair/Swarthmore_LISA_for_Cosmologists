"""
Sensitivity to GW150914
"""
from scipy import integrate
from LISA_Object import LISA
import math
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
integrand = lambda x: 2997.9/(.7/(math.sqrt(((.3*(1+x))**3)+.7)))
Dcom = Mpc*integrate.quad(integrand,0,Z)[0]
DL = (1+Z)*Dcom

#Initial Radiating Frequency, in observer frame
fi = .01774/(1+Z) #Radiates for 4 years from this frequency, in observer frame

#Time to radiate
fr = fi*(1+Z)
Tc = lambda frad: ((5*(c**5))/(256*(fr**(8/3))*((GN*Mr)**(5/3))*(math.pi()**(8/3)))) - ((5*(c**5))/(256*(frad**(8/3))*((GN*Mr)**(5/3))*(math.pi**(8/3))))

#... in observer frame frame (eq 98)
#check that it radiates for 4 years
print(((5*(c**5))/(256*(fr**(8/3))*((GN*Mr)**(5/3))*(math.pi**(8/3))))/((lisa.Year)*(1+Z)))

#characteristic strain (eq 99)
hc = (math.sqrt(2/3)*(((GN*Mr*(1+Z))/(c**3))**(5/6)))/((((math.pi**(2/3))*DL)/c)*(fi**(-1/6)))
hcf = lambda f: hc(f/fi)**(-1/6)
#integrand = lambda 
#SNR = math.sqrt(integrate.quad(
