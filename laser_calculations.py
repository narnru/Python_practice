# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 20:26:51 2021

@author: Nikit
"""
#amplifier/laser parameters
#%%
Length = 1 #m
NumberOfStepZ = 1000
ActiveAreaRadia = 6e-6 #m
ActiveAreaRadiaPump = 62e-6 #m
PowerSignalIn = 0.01 #W
PowerPumpIn = 10 #W
levelLifespan = 1.45e-3 #s (from arxiv:1502.02885)
#%%

import numpy as np
import scipy.constants as scc
import matplotlib.pyplot as plt
import time
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

# length discretisation
#%%
# Z = np.linspace(0, Length, NumberOfStepZ)
# StepZ = Z[1] #m

# signalZ = signalConcentration * np.ones(NumberOfStepZ)
# pumpZ = pumpConcentration * np.ones(NumberOfStepZ)
# inversionZ = np.zeros(NumberOfStepZ)

# #in the folder should be pics related to this part
# #one way calculation

# start = time.time()

# for i in range(1, NumberOfStepZ):
#     inversionZ[i-1] = (concentration * 
#                         (signalZ[i-1]*absorptionSignalCrossection + 
#                         pumpZ[i-1]*absorptionPumpCrossection) /
#                         (1/scc.c/levelLifespan + 
#                         signalZ[i-1]*(absorptionSignalCrossection + emissionSignalCrossection) + 
#                         pumpZ[i-1]*(absorptionPumpCrossection + emissionPumpCrossection)))
#     signalZ[i] = signalZ[i-1] * np.exp((
#         (absorptionSignalCrossection + emissionSignalCrossection)*inversionZ[i-1] - 
#         concentration*absorptionSignalCrossection)*StepZ)
#     pumpZ[i] = pumpZ[i-1] * (1 + signalArea/pumpArea * np.expm1((
#         (absorptionPumpCrossection + emissionPumpCrossection)*inversionZ[i-1] - 
#         concentration*absorptionPumpCrossection)*StepZ))
    
# # Guess what
# print('time', time.time()- start)
# plt.plot(Z, signalZ/signalCoef)
# # plt.plot(Z, pumpZ/pumpCoef)
# testDataSave = signalZ/signalCoef
#%%
Z = np.linspace(0, Length, NumberOfStepZ)
StepZ = Z[1] #m

signalZ = signalConcentration * np.ones(NumberOfStepZ)
pumpZ = pumpConcentration * np.ones(NumberOfStepZ)
signalZback = signalConcentration * np.ones(NumberOfStepZ)
pumpZback = pumpConcentration * np.ones(NumberOfStepZ)
inversionZ = np.zeros(NumberOfStepZ)
inversionStart = np.zeros(NumberOfStepZ)
start = time.time()
lastAttempt = [0]

class crossSection:
    absorptionCrossection = 0
    emissionCrossection = 0
    def __init__(self, absor, emis):
        self.absorptionCrossection = absor
        self.emissionCrossection = emis
    
signalCS = crossSection(absorptionSignalCrossection, emissionSignalCrossection)
pumpCS = crossSection(absorptionPumpCrossection, emissionPumpCrossection)

def calculateInv(signal, pump, sCS:crossSection, pCS:crossSection):
    signal1 = np.zeros(NumberOfStepZ)
    for sig in signal:
        signal1 += sig
    pump1 = np.zeros(NumberOfStepZ)
    for pum in pump:
        pump1 += pum
    return ((signal1*sCS.absorptionCrossection + 
                            pump1*pCS.absorptionCrossection) /
                           (1/scc.c/levelLifespan + 
                            signal1*(sCS.absorptionCrossection + sCS.emissionCrossection) + 
                            pump1*(pCS.absorptionCrossection + pCS.emissionCrossection)))

#attempt to do two way amplifier

        
for j in range(100):
    inversionStart = inversionZ
    for i in range(1, NumberOfStepZ):        
        signalZ[i] = signalZ[i-1] * np.exp((
            (absorptionSignalCrossection + emissionSignalCrossection)*(inversionZ[i] + inversionZ[i-1])/2 - 
            concentration*absorptionSignalCrossection)*StepZ)
        pumpZ[i] = pumpZ[i-1] * (1 + signalArea/pumpArea * np.expm1((
            (absorptionPumpCrossection + emissionPumpCrossection)*(inversionZ[i] + inversionZ[i-1])/2 - 
            concentration*absorptionPumpCrossection)*StepZ))

    inversionZ = (concentration * calculateInv([signalZ, signalZback], [pumpZ, pumpZback], signalCS, pumpCS))
        
    for i in range(1, NumberOfStepZ):
        signalZback[-i-1] =  signalZback[-i] * np.exp((
            (absorptionSignalCrossection + emissionSignalCrossection)*(inversionZ[-i-1] + inversionZ[-i])/2 - 
            concentration*absorptionSignalCrossection)*StepZ)
        pumpZback[-i-1] =  pumpZback[-i] * (1 + signalArea/pumpArea * np.expm1((
            (absorptionPumpCrossection + emissionPumpCrossection)*(inversionZ[-i] + inversionZ[-i-1])/2 - 
            concentration*absorptionPumpCrossection)*StepZ))
    inversionZ = (concentration * calculateInv([signalZ, signalZback], [pumpZ, pumpZback], signalCS, pumpCS))
    # if (j != 0 and j%2 == 1):
    #     inversionZ = (2*inversionZ/3 + inversionStart/3)
    lastAttempt += [max(abs(inversionStart - inversionZ))]
    print(lastAttempt[-1]/max(inversionStart))
    plt.plot(Z, (inversionZ)/max(inversionZ))
    if lastAttempt[-1]/max(inversionStart) < 1e-3:
        break
fig1 = plt.figure()
plt.plot(Z, signalZ/signalCoef)
plt.plot(Z, signalZback/signalCoef)
fig2 = plt.figure()
plt.plot(Z, pumpZ/pumpCoef)
plt.plot(Z, pumpZback/pumpCoef)
fig3 = plt.figure()
plt.plot(Z, inversionZ)
print('time ', time.time()- start)
#%%
test = [signalZ]
for each in test:
    test1 = each