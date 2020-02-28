class GW():

    def __init__(self,M1,M2,Z,f):
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
    self.Dcom = Mpc*integrate.quad(integrand,0,Z)[0]
    self.DL = (1+Z)*Dcom
    
