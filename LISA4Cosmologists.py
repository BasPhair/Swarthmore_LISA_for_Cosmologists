##This script provides the tools needed to calculate LISA sensitivity curves as
##described in "LISA for Cosomolgists"
import math
import csv
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
from numpy.polynomial import hermite
import numpy as np
pi = math.pi

##Retrieves the data to be analyzed. If the data file is not in the same directory as this code be sure to use the full path.
def get_data(file):
    with open(file) as data:
        readData = csv.reader(data, delimiter = "\t")
        data = list(readData)
        data = [(float(data[i][0]),float(data[i][1]),float(data[i][2])) for i in range(len(data))]
    return data


"""
Sensitivity Tables
"""

##Retrieve our data.
Rtab = get_data('LISA_Data.csv')

##Plot the Sensitivity Tables
def SensitivityPlot():
    plt.loglog([float(Rtab[i][0]) for i in range(len(Rtab))], [float(Rtab[i][1])-float(Rtab[i][2]) for i in range(len(Rtab))], 'o', markersize = 1)
    plt.title(r'$\mathsf{LogLog\;\;Plot\;\;of\;\;R}_{A,E}\mathsf{\;\;=\;\;R}_1\mathsf{\;\;-\;\;R}_2$')
    plt.show()  


"""
Definitions
"""

##Constants
#Number of Seconds in a year
Year = 365*24*60*60#Seconds
#Hubble Constant
H0 = (100*.67*(10**3))/(3.086*(10**22))#1/second
#We set the Hubble parameter to h=0.67 (despite the Hubble tension!)

##Mision Criteria
#This is the 6-link design with L = 2.5e6 km arms.
T = 4*Year
SNR5 = 5


"""
Noise Model
"""

#Per the LISA Science Requirements document
fI = .4*(10**(-3))
def SI(f):
    return 5.76*(10**(-48))*(1+((fI/f)**2))
def SII(f):
    return 3.6*(10**(-41))
def Sa(f):
    return ((1/4)*SI(f))/((2*pi*f)**4)
def Ss(f):
    return SII(f)
f2 = 25*(10**(-3))

##Low Frequency Aproximation
#for a stochastic background, Equation 63
def SumI(f):
    return float(math.sqrt(2)*(20/3)*((SI(f)/((2*pi*f)**4))+SII(f))*(1+((f/f2)**2)))
def SumOmegaApprox(f):
    return ((4*(pi**2)*(f**3))/(3*(H0**2)))*SumI(f)
#for a deterministic source, Equation 87
def SumHApprox(f):
    return float((1/2)*(20/3)*((SI(f)/((2*pi*f)**4)) + SII(f))*(1 + ((f/f2)**2)))

##No approximation to the sensitivity curves. (We call it "numerical")
L = 25/3
fstar = 1/(2*pi*L)
def C1(f):
    return (4*Ss(f))+(8*Sa(f)*(1+math.cos(f/fstar)**2))
def C2(f):
    return -((2*Ss(f)+8*Sa(f))*math.cos(f/fstar))

#For a stochastic background, Equation 61
SumOmegaTab = [[Rtab[i][0] for i in range(len(Rtab))],[(1/math.sqrt(((((3*(H0**2))/(4*(pi**2)*((Rtab[i][0])**3)))**2)*(2*(((Rtab[i][1]-Rtab[i][2])/(C1(Rtab[i][0])- C2(Rtab[i][0])))**2)))))for i in range(len(Rtab))]]
SumOmegaNum = interpolate.interp1d(SumOmegaTab[0], SumOmegaTab[1],"cubic")

#for a deterministic source, Equation 85
SumHTab = [[Rtab[i][0] for i in range(len(Rtab))], [.5*((2*((Rtab[i][1] - Rtab[i][2])/(C1(Rtab[i][0]) - C2(Rtab[i][0]))))**(-1)) for i in range(len(Rtab))]]
SumHNum = interpolate.interp1d(SumHTab[0],SumHTab[1])



##Graph Curves
def PlotOmega():
    x = np.logspace(-5,0,1000)
    plt.loglog(x,[SumOmegaNum(i) for i in x], '--r',x,[SumOmegaApprox(i) for i in x])
    plt.title(r'$\Sigma\Omega\mathsf{\;num\;\;vs\;\;}\Sigma\Omega\mathsf{\;approx}$')
    plt.show()

def PlotH():
    x = np.logspace(-5,0,1000)
    plt.loglog(x,[SumHNum(i) for i in x], '--r',x,[SumHApprox(i) for i in x])
    plt.title(r'$\Sigma\mathsf{h}\mathsf{\;num\;\;vs\;\;}\Sigma\mathsf{h}\mathsf{\;approx}$')
    plt.show()


SensitivityPlot()
PlotOmega()
PlotH()



"""
Sensitivity to a flat Stochastic background
"""

fmin = 10**(-5) #.01mHz
fmax = 10**(0) # 1Hz

def integrandNum(f):
    return (1/(SumOmegaNum(f)**2))

def integrandApp(f):
    return (1/(SumOmegaApprox(f)**2))

def Omega_GW_h2(fmin,fmax,t):
    return ((SNR5/math.sqrt(t*integrate.quad(integrandNum,fmin,fmax)[0]))*(.67**2))

##The following is to compare results with the Mathematica notebook with the the mathematica
##commands as comments to increase ease of reference

# SNR5/Sqrt[T NIntegrate[1/(\[CapitalSigma]\[CapitalOmega]approx[f])^2, {f, fmin, fmax}]]
print(SNR5/math.sqrt(T*integrate.quad(integrandApp,fmin,fmax)[0]))

# SNR5/Sqrt[T NIntegrate[1/(\[CapitalSigma]\[CapitalOmega]num[f])^2, {f, fmin, fmax}]]
print(SNR5/math.sqrt(T*integrate.quad(integrandNum,fmin,fmax)[0]))

"""
Num is stilled fucked
"""

#The above is repeated with new fmin and fmax

#(* \[CapitalOmega]_{GW} h^2 = *) % 0.67^2
fmin = 10**(-4)
fmax = 10**(-1)

#SNR5/Sqrt[T NIntegrate[1/(\[CapitalSigma]\[CapitalOmega]approx[f])^2, {f, fmin, fmax}]]
print(SNR5/math.sqrt(T*integrate.quad(integrandApp,fmin,fmax)[0]))

# SNR5/Sqrt[T NIntegrate[1/(\[CapitalSigma]\[CapitalOmega]num[f])^2, {f, fmin, fmax}]]
print(SNR5/math.sqrt(T*integrate.quad(integrandNum,fmin,fmax)[0]))


#SNR5/Sqrt[T3 NIntegrate[1/(\[CapitalSigma]\[CapitalOmega]num[f])^2, {f, fmin, fmax}]]
T3 = 3*Year
print(SNR5/math.sqrt(T3*integrate.quad(integrandNum,fmin,fmax)[0]))



"""
Sensitivity to an arbitrary stochastic background
"""



"""
Sensitivity to a power law stochastic background

Use the method presented in "Sensitivity curves for searches for gravitational-wave backgrounds,"
by Thrane and Romano, Phys.Rev. D88 (2013) no.12, 124032; https://arxiv.org/abs/1310.5300

If Subscript[\[CapitalOmega], GW]= A ((f/fstar))^nt, then what is the minimum detectable A for each value of n ?
"""
fmin = 10**-4 #Hz
fmax = 10**-1 #Hz



def Amin(nt,fmin,fmax):
    x = lambda x: (((x/fstar)**nt)/SumOmegaApprox(x))**2 ##Currently using Approx not Num, Num is used in .nb
    return SNR5/math.sqrt((T*integrate.quad(x, fmin, fmax)[0]))

ntmin = -7/2
ntmax = 9/2

Atab = [[nt, Amin(nt,fmin,fmax)] for nt in np.linspace(ntmin, ntmax, num = 41)]

def Ftab(f):
    return [Atab[i][1]*((f/fstar)**Atab[i][0]) for i in range(len(Atab))]


#Visualize the power-law functions for all values of nt
def Fplot(x_range):
    x = x_range
    ys = []
    for n in range(41):
        ys.append([Ftab(i)[n] for i in x])
    for y in ys:
        plt.loglog(x,y,'k-')
    plt.ylim(10**-14,10**-6)
    plt.show()

#Fplot([10**-4,.1])

"""
Sensitivity to GW150914
"""

"""
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
Mr = ((M1*M2)**3)/((M1+M2)**.2)

#Distance
integrand = lambda x: 2997.9/(.7/(math.sqrt(((.3*(1+x))**3)+.7)))
Dcom = Mpc*integrate.quad(integrand,0,Z)[0]
DL = (1+Z)*Dcom

#Initial Radiating Frequency, in observer frame
fi = .01774/(s/(1+Z)) #Radiates for 4 years from this frequency, in observer frame

#Time to radiate
fr = fi*(1+Z)
Tc = lambda frad: ((5*(c**5))/(256*(fr**(8/3))*((GN*Mr)**(5/3))*(math.pi()**(8/3)))) - ((5*(c**5))/(256*(frad**(8/3))*((GN*Mr)**(5/3))*(math.pi**(8/3))))

#... in observer frame frame (eq 98)
#check that it radiates for 4 years
print(((5*(c**5))/(256*(fr**(8/3))*((GN*Mr)**(5/3))*(math.pi**(8/3))))/((Year)*(1+Z)))
"""




