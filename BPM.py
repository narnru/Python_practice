# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 16:53:15 2021

@author: Nikit
"""

#parameters of the calculations

CalculationRadia = 1000 #um
NumberOfPointsRadia = 501
Wavelength = 1030*1e-3 #um
RefractionIndex = 1
CalculationLength = 1e5 #um
NumberOfPointsLength = 1001
Waist = 100 #um
waistPosition = 0


#meaningful code starts here

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spl
import mpl_interactions.ipyplot as iplt
import time
# derived parameters

k = 2*np.pi/Wavelength #1/mm wavenumber in vacuum
beta = k*RefractionIndex # wavenumber in media
waistRadia = Waist/2
raylaightLength = waistRadia**2 * beta/2 #for 1030nm and 100um waist should be 7.625
StartRadia = waistRadia * (1 + (waistPosition/raylaightLength)**2) ** 0.5

R = np.linspace(0, CalculationRadia, NumberOfPointsRadia)
StepR = R[1] # (i kind of too lazy to calculate it proper way)
Z = np.linspace(0, CalculationLength, NumberOfPointsLength)
StepZ = Z[1]

# curvature radia can be quite tricky if you need to divide by 0
# formula from rp-photonics

if waistPosition == 0:
    curvatureRadia = np.inf
else:
    curvatureRadia = waistPosition*(1+(raylaightLength/waistPosition)**2)

#analytical distribution of gaussian beam from wiki

startDist = (waistRadia/StartRadia
             *np.exp(-R**2/StartRadia**2
             + 1j*(-R**2*beta/curvatureRadia/2 - beta*waistPosition)
             -1j*np.arctan(waistPosition/raylaightLength)))

#for calculation data storage

matrixRZ = np.zeros((NumberOfPointsLength, NumberOfPointsRadia), dtype=complex)
matrixRZ[0,:] = startDist

#okay that's a little bit too hard to explain

A = StepZ/(2*1j*beta)
rp1z = A/R[1:NumberOfPointsRadia-1]/StepR/2 + A/StepR**2
rm1z = -A/R[1:NumberOfPointsRadia]/StepR/2 + A/StepR**2
rz = -A*2*np.ones(NumberOfPointsRadia)/StepR**2
mat = np.diagflat(np.block([0, rp1z]), 1) + np.diagflat(rm1z, -1) + np.diagflat(rz, 0)
mat[1,1] = (mat[1,1] + 4/3*mat[1,0])
mat[1,2] = (mat[1,2] - 1/3*mat[1,0])
matOnes = np.diagflat(np.ones(NumberOfPointsRadia))

matbackward = matOnes[1:,1:] - mat[1:, 1:]/2
matforward = matOnes[1:,1:] + mat[1:, 1:]/2

#matforward = sps.csr_matrix(matforward)
#matbackward = sps.csr_matrix(matbackward)

#x = spl.spsolve(mat, matrixRZ[0,:])

#Нужно аккуратно разобраться что на что умножается. А то трэш какой то

start = time.time()


for i in range(1, NumberOfPointsLength):
    #matrixRZ[i, 1:] = np.linalg.solve((matOnes[1:,1:] - mat[1:, 1:]), matrixRZ[i-1, 1:])
    #matrixRZ[i, 1:] = np.dot((matOnes[1:,1:] + mat[1:, 1:]), matrixRZ[i-1, 1:])
    #matrixRZ[i, 1:] = np.linalg.solve(matbackward, 
    #                                 np.dot(matforward, matrixRZ[i-1, 1:]))
    
    
    
    matrixRZ[i, 0] = 18/11*matrixRZ[i,1] - 9/11*matrixRZ[i,2] + 2/11*matrixRZ[i,3]
print(time.time()- start)
def plotter(l):
    return abs(matrixRZ[round(l)])

fig1, ax1 = plt.subplots()
controls1 = iplt.plot(plotter, l = np.linspace(0, NumberOfPointsLength), )