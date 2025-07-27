import numpy as np
from Particle import Particle
import math
import copy
from matplotlib import pyplot as plt

from matplotlib import rcParams

rcParams["font.family"] = "serif"  # change to default to a serif font

rcParams["font.size"] = 18  # change the default font size
rcParams["figure.figsize"] = [9.7, 6]  # change the default figure size
rcParams["figure.autolayout"] = True  # automatically apply tight_layout
rcParams["figure.figsize"] = [9.7, 6]  # change the default figure size

deltaT= 3600*24 #Time step
tend = 31536000*10 # time for 1 year x number of years
niter = int(tend/deltaT) #number of itterations

#Plotting Angular Momentum with euler method

Data1 = np.load("simulation_results_euler.npy", allow_pickle=True) #open data from simulation that used euler method
fig1, ax1 = plt.subplots()
time_L1 = []
FracChangeL_1 = []
for i in range (0,niter): # make a list of all the time values and fractional change at each iteration 
    t = Data1[i][0]
    time_L1.append(t)
    L = Data1[i][6]
    FracChangeL_1.append(L)

ax1.plot(time_L1, FracChangeL_1, linestyle ="-" , color ="k", label ="Euler")

#Plotting Angular Momentum using Euler-Cromer Method
Data2 = np.load("simulation_results_euler_cromer.npy", allow_pickle=True) # open data from simulation that used euler cromer method
time_L2 = []
FracChangeL_2 = []
for i in range (0,niter): 
    t = Data2[i][0]
    time_L2.append(t)
    L = Data2[i][6]
    FracChangeL_2.append(L)

ax1.plot(time_L2, FracChangeL_2, linestyle ="-" , color ="r", label ="Euler Cromer") #plot to the same axes for comparison
ax1.set_xlabel(r"$t[s]$")
ax1.set_ylabel(r"$ΔL/L$")
plt.legend()
plt.show()