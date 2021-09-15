# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 20:26:51 2021

@author: Nikit
"""
#amplifier/laser parameters

Length = 1 #m
NumberOfStepZ = 1000
ActiveAreaRadia = 6e-6 #m
ActiveAreaRadiaPump = 62e-6 #m
PowerSignalIn = 0.01 #W
PowerPumpIn = 10 #W

import numpy as np
import scipy.constants as scc
import matplotlib.pyplot as plt

#Signal derived parameters

emissionSignalCrossection = 3.2e-25 #1/m^2
absorptionSignalCrossection = 1.7e-26 #1/m^2
signalWavelength = 1030e-9 #m

signalArea = ActiveAreaRadia**2 * scc.pi #m^2
signalCoef = 1 / signalArea / scc.h / scc.c**2 * signalWavelength
signalConcentration= (PowerSignalIn * signalCoef) #1/m^3

#Pump derived parameters

emissionPumpCrossection = 1.3e-24 #1/m^2
absorptionPumpCrossection = 1.45e-24 #1/m^2
pumpWavelength = 975e-9 #m

pumpArea = ActiveAreaRadiaPump**2 * scc.pi #m^2
pumpCoef = 1 / pumpArea / scc.h / scc.c**2 * pumpWavelength
pumpConcentration= PowerPumpIn * pumpCoef #1/m^3

# magic

concentrationppm = 1000 #ppm
concentration = concentrationppm * 6.62e22 #1/m^3

levelLifespan = 1.45e-3 #s (from arxiv:1502.02885)

# length discretisation

Z = np.linspace(0, Length, NumberOfStepZ)
StepZ = Z[1] #m

signalZ = signalConcentration * np.ones(NumberOfStepZ)
pumpZ = pumpConcentration * np.ones(NumberOfStepZ)
inversionZ = np.zeros(NumberOfStepZ)

#in the folder should be pic related to this part

for i in range(1, NumberOfStepZ):
    inversionZ[i-1] = (concentration * 
                       (signalZ[i-1]*absorptionSignalCrossection + 
                        pumpZ[i-1]*absorptionPumpCrossection) /
                       (1/scc.c/levelLifespan + 
                        signalZ[i-1]*(absorptionSignalCrossection + emissionSignalCrossection) + 
                        pumpZ[i-1]*(absorptionPumpCrossection + emissionPumpCrossection)))
    signalZ[i] = signalZ[i-1] * np.exp((
        (absorptionSignalCrossection + emissionSignalCrossection)*inversionZ[i-1] - 
        concentration*absorptionSignalCrossection)*StepZ)
    pumpZ[i] = pumpZ[i-1] * (1 + signalArea/pumpArea * np.expm1((
        (absorptionPumpCrossection + emissionPumpCrossection)*inversionZ[i-1] - 
        concentration*absorptionPumpCrossection)*StepZ))
    
# Guess what

plt.plot(Z, signalZ/signalCoef)
plt.plot(Z, pumpZ/pumpCoef)

    
