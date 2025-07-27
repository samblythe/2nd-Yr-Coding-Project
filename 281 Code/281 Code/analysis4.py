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

deltaT= 3600 #Time step
tend = 31536000*10 # time for 1 year x number of years
niter = int(tend/deltaT)

#Plotting Total Energy with euler method
Data1 = np.load("simulation_results_euler.npy", allow_pickle=True) #open data from simulation that used euler method
fig1, ax1 = plt.subplots()
time_1 = []
FracChangeE_1 = []
for i in range (0,niter): #make a list of time and fractional change to plot
    t = Data1[i][0]
    time_1.append(t)
    E = Data1[i][7]
    FracChangeE_1.append(E)

ax1.plot(time_1, FracChangeE_1, linestyle ="-" , color ="c", label ="Euler") #plot time and fractional change


Data2 = np.load("simulation_results_euler_cromer.npy", allow_pickle=True) #open data from simulation that used euler cromer method

time_2 = []
FracChangeE_2 = []
for i in range (0,niter): #make a list of time and fractional change to plot
    t = Data2[i][0]
    time_2.append(t)
    E = Data2[i][7]
    FracChangeE_2.append(E)

ax1.plot(time_2, FracChangeE_2, linestyle ="-" , color ="r", label ="Euler Cromer") #plot on same axis for comparison with other method

plt.legend()
plt.show()