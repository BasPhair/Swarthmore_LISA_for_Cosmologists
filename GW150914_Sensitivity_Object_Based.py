"""
Sensitivity of LISA to GW150914 using LISA class and Binary class.
"""
from LISA_Object import LISA
from GW_Object import GW
M1 = 36
M2 = 29
Z = .09

GW150914 = GW(M1,M2,Z)
lisa = LISA('LISA_Data.csv')

print(GW150914.SNR_Num(lisa))

print(GW150914.SNR_Approx(lisa))

GW150914.Plot_LISA_Sensitivity(lisa, title = 'Sensitivity of LISA to GW150914')
