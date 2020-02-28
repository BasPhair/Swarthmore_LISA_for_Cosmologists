from LISA_Object import LISA

lisa = LISA('LISA_Data.csv')


lisa.SensitivityPlot()

lisa.PlotOmega()

lisa.PlotH()

lisa.PlotF([10**-4,.1])

lisa.OmegaFPlot([10**-4,.1],100)
