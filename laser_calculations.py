# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 20:26:51 2021
"""

import numpy as np
import scipy.constants as scc
import matplotlib.pyplot as plt
import time
import copy
import math

class Emission:
    """Class for storing everything connected to emission pump or signal
    
    todo - Probably it is better to go for intensities instead of photon concentrations
    at least easier to understand.
    """
    
    Area : float() #m**2
    WtoConcCoef : float() #1/m**3/W
    Wavelength : float() #m
    StartPower : float() #W
    Distrib : np.ndarray(100) #1/m**3
    AbsCS : float() #1/m**2
    EmiCS : float() #1/m**2
    
    def __init__(self, AreaRadia, Wavelength, Power, NumberOfStepZ, absorptionCS,
                 emissionCS):
        """Setuping a starting parameters of emission
        
        Conventions:
        AreaRadia - m, Wavelength - m, Power - W, crossections - 1/m**2
        """
        self.Area = AreaRadia**2 * scc.pi
        self.Wavelength = Wavelength
        self.WtoConcCoef = 1 / self.Area / scc.h / scc.c**2 * self.Wavelength
        self.StartPower = Power
        self.Distrib = Power * self.WtoConcCoef * np.ones(NumberOfStepZ, dtype = np.float32())
        self.AbsCS = absorptionCS
        self.EmiCS = emissionCS
    
    
    
    def set_power(self, Power):
        """For power changes, since it affects photons concentrations"""
        self.Distrib = Power * self.WtoConcCoef * np.ones(self.Distrib.size)
        self.StartPower = Power
       
    def calculateDistrib(self, inversion, concentration, StepZ, AreaCoef = 1):
        for i in range(1, self.Distrib.size):    
            self.Distrib[i] = self.Distrib[i-1] * (1 + AreaCoef * ((
                (self.AbsCS + self.EmiCS)*(inversion[i]+inversion[i-1])/2 - 
                concentration*self.AbsCS)*StepZ))
            
    
    def calculateDistribRev(self, inversion, concentration, StepZ, AreaCoef = 1):
        for i in range(1, self.Distrib.size):    
            self.Distrib[-i-1] = self.Distrib[-i] * (1 + AreaCoef * ((
                (self.AbsCS + self.EmiCS)*(inversion[-i]+inversion[-i-1])/2 - 
                concentration*self.AbsCS)*StepZ))

        
class Media:
    """Class for storing everything connected to the active media
    
    todo - hardcoding that it is quartz in concentration isn't good idea
    """
    Length : float() #m
    NumZ : int()
    StepZ : float() #m
    Points : np.ndarray(100) #m
    Inversion : np.ndarray(100) #1/m**3
    Concentration : float() #1/m**3
    Tau : float() #s
    def __init__(self, Length, NumberOfStepZ, ConcentrationPPM, levelLifespan):
        """Setuping parameters of media
        
        Conventions:
            Length - m, Concentration - ppm, Lifespan - s
        """
        self.Length = Length
        self.NumZ = NumberOfStepZ
        self.Points = np.linspace(0, Length, self.NumZ)
        self.StepZ = self.Points[1] 
        self.Inversion = np.zeros(NumberOfStepZ)
        self.Concentration = ConcentrationPPM * 6.62e22 #1/m^3
        self.Tau = levelLifespan 
    
    def calculateInv(self, signal : list):
        """Calculation of inversion based on signal information. 
        Fun fact - you just make a list from signal and pump and it works"""
        signalAbs = np.zeros(self.NumZ)
        signalEmi = np.zeros(self.NumZ)
        for sig in signal:
            signalAbs += sig.Distrib * sig.AbsCS
            signalEmi += sig.Distrib * sig.EmiCS
        self.Inversion = (self.Concentration * signalAbs / (1/scc.c/self.Tau + 
                                signalAbs + signalEmi))

        
time_in_one_way = 0
def one_way_amplif(Signal : Emission, Pump : Emission, media : Media):
    """ Calculation of one way amplifier. In the folder with this file should be
    picture about physics. Basically it is fully defined from one side differential
    equation. Use only for tests.
    """
    start = time.time()
    for i in range(1, media.NumZ):
        media.Inversion[i-1] = (media.Concentration * 
                            (Signal.Distrib[i-1]*Signal.AbsCS + 
                            Pump.Distrib[i-1]*Pump.AbsCS) /
                            (1/scc.c/media.Tau + 
                            Signal.Distrib[i-1]*(Signal.AbsCS + Signal.EmiCS) + 
                            Pump.Distrib[i-1]*(Pump.AbsCS + Pump.EmiCS)))
        Signal.Distrib[i] = Signal.Distrib[i-1] * math.exp(media.StepZ*(
            (Signal.AbsCS + Signal.EmiCS)*media.Inversion[i-1] 
            - media.Concentration*Signal.AbsCS))
        Pump.Distrib[i] = Pump.Distrib[i-1] * math.exp(Signal.Area/Pump.Area*((
            (Pump.AbsCS + Pump.EmiCS)*media.Inversion[i-1] 
            - media.Concentration*Pump.AbsCS)*media.StepZ))
    global time_in_one_way 
    time_in_one_way += -start + time.time()
    return


def two_way_amplifier(Signal : Emission, SignalBack : Emission,  
                      Pump : Emission, PumpBack : Emission, media : Media):
    """Calculation of two-way amplifier
    
    todo - documentation, adequate reason for finishing calculations
    """    
    lastAttempt = []
    
    for j in range(100):
        inversionStart = media.Inversion
    
        Signal.calculateDistrib(media.Inversion, media.Concentration, media.StepZ)
        Pump.calculateDistrib(media.Inversion, media.Concentration, media.StepZ, Signal.Area/Pump.Area)
        
        media.calculateInv([Signal, SignalBack, Pump, PumpBack])
            
        SignalBack.calculateDistribRev(media.Inversion, media.Concentration, media.StepZ)
        PumpBack.calculateDistribRev(media.Inversion, media.Concentration, media.StepZ, Signal.Area/Pump.Area)
        
        media.calculateInv([Signal, SignalBack, Pump, PumpBack])
    
        lastAttempt += [max(abs(inversionStart - media.Inversion))]
        print(lastAttempt[-1]/max(inversionStart))
        if lastAttempt[-1]/max(inversionStart) < 1e-5:
            break
  
def main():
  
    #amplifier/laser parameters
    
    Length = 5 #m
    NumberOfStepZ = 10000
    ActiveAreaRadia = 6e-6 #m
    ActiveAreaRadiaPump = 62.5e-6 #m
    PowerSignalIn = 0.01 #W
    PowerPumpIn = 10 #W
    levelLifespan = 1.45e-3 #s (from arxiv:1502.02885)
    
    emissionSignalCrossection = 3.2e-25 #1/m^2
    absorptionSignalCrossection = 1.7e-26 #1/m^2
    signalWavelength = 1030e-9 #m
    
    emissionPumpCrossection = 1.3e-24 #1/m^2
    absorptionPumpCrossection = 1.45e-24 #1/m^2
    pumpWavelength = 975e-9 #m
    
    concentrationppm = 1000 #ppm
      
    #derived parameters 
    Signal = Emission(ActiveAreaRadia, signalWavelength, PowerSignalIn, NumberOfStepZ,
                      absorptionSignalCrossection, emissionSignalCrossection)
    Pump = Emission(ActiveAreaRadiaPump, pumpWavelength, PowerPumpIn, NumberOfStepZ,
                    absorptionPumpCrossection, emissionPumpCrossection)
    media = Media(Length, NumberOfStepZ, concentrationppm, levelLifespan)
    
    one_way_amplif(Signal, Pump, media)
    old_attempt = Signal.Distrib/Signal.WtoConcCoef
    
    Signal = Emission(ActiveAreaRadia, signalWavelength, PowerSignalIn, NumberOfStepZ,
                      absorptionSignalCrossection, emissionSignalCrossection)
    Pump = Emission(ActiveAreaRadiaPump, pumpWavelength, PowerPumpIn, NumberOfStepZ,
                    absorptionPumpCrossection, emissionPumpCrossection)
    media = Media(Length, NumberOfStepZ, concentrationppm, levelLifespan)
    
    
    SignalBack = copy.deepcopy(Signal)
    SignalBack.set_power(Signal.StartPower*0)
    PumpBack = copy.deepcopy(Pump) 
    PumpBack.set_power(Pump.StartPower*0)
    
    #attempt to do two way amplifier
    
    # start = time.time()
    # two_way_amplifier(Signal, SignalBack, Pump, PumpBack, media)
    # plt.figure()
    # plt.plot(media.Points, Signal.Distrib/Signal.WtoConcCoef)
    # plt.plot(media.Points, old_attempt)
    # plt.figure()
    # plt.plot(media.Points, Pump.Distrib/Pump.WtoConcCoef)
    # plt.plot(media.Points, PumpBack.Distrib/PumpBack.WtoConcCoef)
    # plt.figure()
    # plt.plot(media.Points, media.Inversion)
    # print('time ', time.time()- start)

def test_stepz_one_way():
    Length = 5 #m
    ActiveAreaRadia = 6e-6 #m
    ActiveAreaRadiaPump = 62.5e-6 #m
    PowerSignalIn = 0.01 #W
    PowerPumpIn = 10 #W
    levelLifespan = 1.45e-3 #s (from arxiv:1502.02885)
    
    emissionSignalCrossection = 3.2e-25 #1/m^2
    absorptionSignalCrossection = 1.7e-26 #1/m^2
    signalWavelength = 1030e-9 #m
    
    emissionPumpCrossection = 1.3e-24 #1/m^2
    absorptionPumpCrossection = 1.45e-24 #1/m^2
    pumpWavelength = 975e-9 #m
    
    concentrationppm = 1000 #ppm
    resultsSP = []
    resultsInt = [[],[],[]]
    start = time.time()
    for N in np.linspace(1, 6, 20):
        Signal = Emission(ActiveAreaRadia, signalWavelength, PowerSignalIn, round(10**N),
                          absorptionSignalCrossection, emissionSignalCrossection)
        Pump = Emission(ActiveAreaRadiaPump, pumpWavelength, PowerPumpIn, round(10**N),
                        absorptionPumpCrossection, emissionPumpCrossection)
        media = Media(Length, round(10**N), concentrationppm, levelLifespan)
        one_way_amplif(Signal, Pump, media)
        resultsSP += [[media.Points, Signal.Distrib/Signal.WtoConcCoef]]
        resultsInt[0] += [N]
        # resultsInt[1] += [np.trapz(Signal.Distrib/Signal.WtoConcCoef, media.Points)]
        resultsInt[1] += [Signal.Distrib[-1]/Signal.WtoConcCoef]
        if len(resultsInt[1]) > 1:
            resultsInt[2] += [abs(resultsInt[1][-1] - resultsInt[1][-2])] 
    for res in resultsSP:
        plt.plot(res[0], res[1])
    plt.figure()
    plt.yscale("log")
    # plt.xscale("log")
    plt.plot(resultsInt[0][1:], resultsInt[2])
    print(resultsInt[2])
    print(resultsInt[1])
    print('time ', time.time()- start)
        
test_stepz_one_way() #looks like second order on z

def test_stepz_two_way():
    Length = 5 #m
    ActiveAreaRadia = 6e-6 #m
    ActiveAreaRadiaPump = 62.5e-6 #m
    PowerSignalIn = 0.01 #W
    PowerPumpIn = 10 #W
    levelLifespan = 1.45e-3 #s (from arxiv:1502.02885)
    
    emissionSignalCrossection = 3.2e-25 #1/m^2
    absorptionSignalCrossection = 1.7e-26 #1/m^2
    signalWavelength = 1030e-9 #m
    
    emissionPumpCrossection = 1.3e-24 #1/m^2
    absorptionPumpCrossection = 1.45e-24 #1/m^2
    pumpWavelength = 975e-9 #m
    
    concentrationppm = 1000 #ppm
    # resultsSP = []
    resultsInt = [[],[],[]]
    start = time.time()
    for N in range(1, 6):
        
        Signal = Emission(ActiveAreaRadia, signalWavelength, PowerSignalIn, 10**N,
                          absorptionSignalCrossection, emissionSignalCrossection)
        Pump = Emission(ActiveAreaRadiaPump, pumpWavelength, PowerPumpIn, 10**N,
                        absorptionPumpCrossection, emissionPumpCrossection)
        SignalBack = copy.deepcopy(Signal)
        SignalBack.set_power(Signal.StartPower*0)
        PumpBack = copy.deepcopy(Pump) 
        PumpBack.set_power(Pump.StartPower*0)
        media = Media(Length, 10**N, concentrationppm, levelLifespan)
        two_way_amplifier(Signal, SignalBack, Pump, PumpBack, media)
        # resultsSP += [[media.Points, Signal.Distrib/Signal.WtoConcCoef]]
        resultsInt[0] += [N]
        # resultsInt[1] += [np.trapz(Signal.Distrib/Signal.WtoConcCoef, media.Points)]
        resultsInt[1] += [Signal.Distrib[-1]/Signal.WtoConcCoef]
        if len(resultsInt[1]) > 1:
            resultsInt[2] += [abs(resultsInt[1][-1] - resultsInt[1][-2])] 
    # for res in resultsSP:
    #     plt.plot(res[0], res[1])
    plt.yscale("log")
    plt.plot(resultsInt[0][1:], resultsInt[2])
    plt.figure()
    plt.plot(resultsInt[0], resultsInt[1])
    print(resultsInt[2])
    print('time ', time.time()- start)
# test_stepz_two_way()
# RefractionLeft = 1
# RefractionRight = 0.6
# start = time.time()
# att = [0];

# for j in range(100):
#     inversionStart = inversionZ
#     for i in range(1, NumberOfStepZ):        
#         signalZ[i] = signalZ[i-1] * np.exp((
#             (absorptionSignalCrossection + emissionSignalCrossection)*(inversionZ[i] + inversionZ[i-1])/2 - 
#             concentration*absorptionSignalCrossection)*StepZ)
#         pumpZ[i] = pumpZ[i-1] * (1 + signalArea/pumpArea * np.expm1((
#             (absorptionPumpCrossection + emissionPumpCrossection)*(inversionZ[i] + inversionZ[i-1])/2 - 
#             concentration*absorptionPumpCrossection)*StepZ))

#     inversionZ = (concentration * calculateInv([signalZ, signalZback], [pumpZ, pumpZback], signalCS, pumpCS))
#     signalZback[-1] = signalZ[-1] * RefractionRight
#     for i in range(1, NumberOfStepZ):
#         signalZback[-i-1] =  signalZback[-i] * np.exp((
#             (absorptionSignalCrossection + emissionSignalCrossection)*(inversionZ[-i-1] + inversionZ[-i])/2 - 
#             concentration*absorptionSignalCrossection)*StepZ)
#         pumpZback[-i-1] =  pumpZback[-i] * (1 + signalArea/pumpArea * np.expm1((
#             (absorptionPumpCrossection + emissionPumpCrossection)*(inversionZ[-i] + inversionZ[-i-1])/2 - 
#             concentration*absorptionPumpCrossection)*StepZ))
#     inversionZ = (concentration * calculateInv([signalZ, signalZback], [pumpZ, pumpZback], signalCS, pumpCS))
#     signalZ[0] = signalZback[0]*RefractionLeft
#     lastAttempt += [max(abs(inversionStart - inversionZ))]
#     att += [(signalZ[-1]-signalZback[-1])/signalCoef]
#     print(abs(att[-1]-att[-2]))
#     if abs(att[-1]-att[-2]) < 1e-10:
#         break
# totalEnergy = ((pumpZ[-1]+pumpZback[0])/pumpCoef +
#            (-signalZ[0]+signalZback[0])/signalCoef + 
#            (signalZ[-1]-signalZback[-1])/signalCoef
#            + inversionZ.sum()/levelLifespan*StepZ*signalArea*scc.h*scc.c/signalWavelength)

# fig1 = plt.figure()
# plt.plot(Z, signalZ/signalCoef)
# plt.plot(Z, signalZback/signalCoef)
# fig2 = plt.figure()
# plt.plot(Z, pumpZ/pumpCoef)
# plt.plot(Z, pumpZback/pumpCoef)
# fig3 = plt.figure()
# plt.plot(Z, inversionZ)
# print('time ', time.time()- start)