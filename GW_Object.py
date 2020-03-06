"""
Need to ask about what potential inputs are.
"""
from scipy import integrate
import numpy as np
import math
import matplotlib.pyplot as plt
from LISA_Object import LISA
class GW():

    def __init__(self,M1,M2,Z,f=None):
        #Definitions, units
        self.GN = 6.67*(10**-11) # m^3/(kg s^2)
        self.c = 3*(10**8) #m/s
        self.Msol = 2*(10**30) #kg
        self.Mpc = 3.086*(10**22) #m
        self.M1 = M1*self.Msol
        self.M2 = M2*self.Msol
        self.Z = Z
        self.Mr = (((self.M1*self.M2)**3)/(self.M1+self.M2))**.2
        integrand = lambda x: 2997.9/.7/math.sqrt((.3*(1+x)**3)+.7)
        self.Dcom = self.Mpc*integrate.quad(integrand,0,Z)[0]
        self.DL = (1+Z)*self.Dcom
        self.fi = f/(1+Z) if f is not None else .01774/(1+Z)
        self.fr = self.fi*(1+Z)
        self.Tc = lambda frad: ((5*(self.c**5))/(256*(self.fr**(8/3))*((self.GN*self.Mr)**(5/3))*(math.pi()**(8/3)))) - ((5*(self.c**5))/(256*(frad**(8/3))*((self.GN*self.Mr)**(5/3))*(math.pi**(8/3))))
        self.hc = ((math.sqrt(2/3)*(((self.GN*self.Mr*(1+Z))/(self.c**3))**(5/6)))*(self.fi**(-1/6)))/((math.pi**(2/3))*self.DL/self.c)
        self.hcf = lambda f: self.hc*((f/self.fi)**(-1/6))

    def SNR_Num(self, lisa):
        integrand = lambda f: ((self.hcf(f)/(2*f))**2)/lisa.SumHNum(f)
        SNR = math.sqrt(integrate.quad(integrand, self.fi, 1)[0])
        return SNR

    def SNR_Approx(self, lisa):
        integrand = lambda f: ((self.hcf(f)/(2*f))**2)/lisa.SumHApprox(f)
        SNR = math.sqrt(integrate.quad(integrand, self.fi, 1)[0])
        return SNR

    def Plot_LISA_Sensitivity(self, lisa, title = None):
        x = np.logspace(-5,0,1000)
        x2 = np.logspace(math.log(self.fi,10),math.log(.05,10),200)
        plt.loglog(x,[math.sqrt(4*i*lisa.SumHApprox(i)) for i in x], x2,[self.hcf(i) for i in x2],'--r')
        if title == None:
            plt.title('Sensitivity of LISA to GW')
        else:
            plt.title(title)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('hc[strain]')
        plt.show()
