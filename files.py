# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 19:23:52 2021
"""



import csv
import numpy as np
import matplotlib.pyplot as plt

def csv_crossection_reader(file_name:str, wavelengthPosition = 0, 
                           emissionCrossSectionPosition = 1, 
                           absorptionCrossSectionPosition = 2, 
                           delim = ';', start = 1):
    """Handling a reading crossections data from csv files. Expected to be used for amplifier calculations
    
    Returns a list of 3 ndarrays. Zero one is wavelength, first is absorption, second is emission.
    Default delimiter for excel created csv.
    Default start expects non float values in first line
    Everything else is self explanatory.
    """
    
    data = csv.reader(open(file_name), delimiter = delim)
    wavelength = list()
    absorptionCrossSection = list()
    emissionCrossSection = list()
    for i in range(start):
        data.__next__()
    for elem in data:
        print(elem)
        wavelength += [float(elem[wavelengthPosition])] 
        absorptionCrossSection += [float(elem[absorptionCrossSectionPosition])]
        emissionCrossSection += [float(elem[emissionCrossSectionPosition])]
    
    return [np.array(wavelength), np.array(absorptionCrossSection),
            np.array(emissionCrossSection)]


def reader_test():
    crossections = csv_crossection_reader('CrossSection.csv', 0, 9, 10, start = 44)
    plt.plot(crossections[0], crossections[1])
    plt.plot(crossections[0], crossections[2])
    
#reader_test()