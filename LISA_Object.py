##This script provides the tools needed to calculate LISA sensitivity curves as
##described in "LISA for Cosomolgists"
import math
import csv
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
from numpy.polynomial import hermite
import numpy as np
pi = math.pi

class LISA():


    def __init__(self,file):
        self.Rtab = []
        self.Year = 365*24*60*60#Seconds
        self.H0 = (100*.67*(10**3))/(3.086*(10**22))#(1/second) We set the Hubble parameter to h=0.67 (despite the Hubble tension!)

        ##Mision Criteria
        #This is the 6-link design with L = 2.5e6 km arms.
        self.T = 4*self.Year
        self.SNR5 = 5
        self.fI = .4*(10**(-3))
        self.f2 = 25*(10**(-3))
        self.L = 25/3
        self.fstar = 1/(2*pi*self.L)

        ##Tabs
        self.get_data(file)
        self.SumOmegaTab = [[self.Rtab[i][0] for i in range(len(self.Rtab))],[(1/math.sqrt(((((3*(self.H0**2))/(4*(pi**2)*((self.Rtab[i][0])**3)))**2)*(2*(((self.Rtab[i][1]-self.Rtab[i][2])/(self.C1(self.Rtab[i][0])- self.C2(self.Rtab[i][0])))**2)))))for i in range(len(self.Rtab))]]
        self.SumHTab = [[self.Rtab[i][0] for i in range(len(self.Rtab))], [.5*((2*((self.Rtab[i][1] - self.Rtab[i][2])/(self.C1(self.Rtab[i][0]) - self.C2(self.Rtab[i][0]))))**(-1)) for i in range(len(self.Rtab))]]

        self.Atab = [[nt, self.Amin(nt,10**-4,10**-1)] for nt in np.linspace(-7/2, 9/2, num = 41)]

    ##Retrieves the data to be analyzed. If the data file is not in the same directory as this code be sure to use the full path.
    def get_data(self, file):
        with open(file) as LISA_data:
            readData = csv.reader(LISA_data, delimiter = "\t")
            Read_data = list(readData)
            self.Rtab = [(float(Read_data[i][0]),float(Read_data[i][1]),float(Read_data[i][2])) for i in range(len(Read_data))]


    def SensitivityPlot(self):
        plt.loglog([float(self.Rtab[i][0]) for i in range(len(self.Rtab))], [float(self.Rtab[i][1])-float(self.Rtab[i][2]) for i in range(len(self.Rtab))], 'o', markersize = 1)
        plt.title(r'$\mathsf{LogLog\;\;Plot\;\;of\;\;R}_{A,E}\mathsf{\;\;=\;\;R}_1\mathsf{\;\;-\;\;R}_2$')
        plt.show()
        
    def PlotOmega(self):
        x = np.logspace(-5,0,1000)
        plt.loglog(x,[self.SumOmegaNum(i) for i in x],label=(r'$\Sigma\Omega\mathsf{\;num}$'))
        plt.loglog(x,[self.SumOmegaApprox(i) for i in x],'--r',label=(r'$\Sigma\Omega\mathsf{\;approx}$'))
        plt.title(r'$\Sigma\Omega\mathsf{\;num\;\;vs\;\;}\Sigma\Omega\mathsf{\;approx}$')
        plt.xlabel('f(Hz)')
        plt.ylabel('Strain')
        plt.legend()
        plt.show()

    def PlotH(self):
        x = np.logspace(-5,0,1000)
        plt.loglog(x,[self.SumHNum(i) for i in x],label=r'$\Sigma\mathsf{h}\mathsf{\;num}$')
        plt.loglog(x,[self.SumHApprox(i) for i in x],'--r',label=r'$\Sigma\mathsf{h}\mathsf{\;approx}$')
        plt.title(r'$\Sigma\mathsf{h}\mathsf{\;num\;\;vs\;\;}\Sigma\mathsf{h}\mathsf{\;approx}$')
        plt.xlabel('f(Hz)')
        plt.ylabel('Strain')
        plt.legend()
        plt.show()

    """
    Noise Model
    """
    def SI(self,f):#equation 53
        return 5.76*(10**(-48))*(1+((self.fI/f)**2))

    def SII(self,f):#equation 54
        return 3.6*(10**(-41))
    
    def Sa(self, f):#equation 52
        return self.SI(f)/(4*((2*f*math.pi)**4))
    
    def Ss(self, f):#equation 52
        return self.SII(f)

    ##Low Frequency Approximation
    def SumI(self, f):
        return float(math.sqrt(2)*(20/3)*((self.SI(f)/((2*math.pi*f)**4))+self.SII(f))*(1+((f/self.f2)**2)))

    def SumOmegaApprox(self, f):
        return self.SumI(f)*((4*(math.pi**2)*(f**3))/(3*(self.H0**2)))

    def SumHApprox (self,f):
        return (1/2)*(20/3)*((self.SI(f)/((2*math.pi*f)**4))+self.SII(f))*(1+((f/self.f2)**2))

    def C1(self,f):
        return (4*self.Ss(f))+(8*self.Sa(f)*(1+math.cos(f/self.fstar)**2))

    def C2(self,f):
        return -((2*self.Ss(f)+8*self.Sa(f))*math.cos(f/self.fstar))

    
    def SumOmegaNum(self,f):
        return interpolate.interp1d(self.SumOmegaTab[0], self.SumOmegaTab[1],"cubic")(f)

    def SumHNum(self,f):
        return interpolate.interp1d(self.SumHTab[0],self.SumHTab[1])(f)

    def Amin(self,nt,fmin,fmax):
        integrand = lambda x: (((x/self.fstar)**nt)/self.SumOmegaNum(x))**2
        return self.SNR5/math.sqrt((self.T*integrate.quad(integrand, fmin, fmax)[0]))

    def Ftab(self,f):
        return [self.Atab[i][1]*((f/self.fstar)**self.Atab[i][0]) for i in range(len(self.Atab))]

    def PlotF(self,x_range):
        x = x_range
        ys = []
        for n in range(len(self.Atab)):
            ys.append([self.Ftab(i)[n] for i in x])
        for y in ys:
            plt.loglog(x,y,'k-')
        plt.ylim(10**-14,10**-6)
        plt.title('Plot of Power Law Functions')
        plt.xlabel('f(Hz)')
        plt.ylabel('Strain')
        plt.show()

    def FLogOmega(self,x_range,inc):
        ys = []
        xs = np.logspace(-4,-1,inc)
        for x in xs:
            ys.append(math.log(max(self.Ftab(x)),10))
        return interpolate.interp1d(xs,ys)

    def OmegaFPlot(self,x_range,inc):
        FLogO = self.FLogOmega(x_range,inc)
        ys = []
        xs = np.logspace(-4,-1,inc)
        ys = [FLogO(f) for f in xs]
        ys = [10**y for y in ys]
        plt.loglog(xs,ys)
        yval = self.Amin(0,10**-4,.1)
        plt.loglog(xs,[yval for x in xs],'r--')
        plt.title('Plot "graph (like Figure 3)"')
        plt.xlabel('f(Hz)')
        plt.ylabel('Strain')
        plt.show()





        
