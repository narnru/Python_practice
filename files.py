# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 19:23:52 2021

@author: Nikit
"""
import csv
import numpy as np

data = csv.reader(open('CrossSection.csv'), delimiter = ';', dialect = 'excel')
wavelength = list()
absorptionCrossSection = list()
emissionCrossSection = list()
data.__next__()
for elem in data:
    print(elem)
    wavelength += [float(elem[0])] 
    absorptionCrossSection += [float(elem[10])]
    emissionCrossSection += [float(elem[9])]


import matplotlib.pyplot as plt

plt.plot(wavelength, absorptionCrossSection)
plt.plot(wavelength, emissionCrossSection)