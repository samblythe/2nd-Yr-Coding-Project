import numpy as np
#from Particle import Particle
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
niter = int(tend/deltaT) #Nnumber of itterations


#Plotting Orbits using Euler Method
Data1 = np.load("simulation_results_euler.npy", allow_pickle=True) #load the data file from simulation using the euler method

#Plotting Earth Orbit with Euler Method
fig1, ax1 = plt.subplots()
xe = []
ye = []
for i in range (0,niter): #this loop creates a list of all the x and y values at each itteration step so they can be plotted
    earth_position = Data1[i][2].position 
    xe.append(earth_position[0])
    ye.append(earth_position[1])
ax1.plot(xe, ye, linestyle ="-" , color ="b", label ="Earth") #plotting x and y values

#Plotting Sun Orbit with Euler Method
xs=[]
ys=[]
for i in range (0,niter):
    sun_position = Data1[i][1].position
    xs.append(sun_position[0])
    ys.append(sun_position[1])
ax1.plot(xs, ys,linestyle = "-", color="k", label="Sun")

#Plotting Jupiter Orbit with Euler Method
xj =[]
yj =[]
for i in range (0,niter):
    jupiter_position = Data1[i][3].position
    xj.append(jupiter_position[0])
    yj.append(jupiter_position[1])
ax1.plot(xj, yj,linestyle="-",color="g", label="Jupiter")

#Plotting Mars Orbit with Euler Method
xm = []
ym = []
for i in range (0,niter): 
    mars_position = Data1[i][4].position
    xm.append(mars_position[0])
    ym.append(mars_position[1])
ax1.plot(xm, ym, linestyle ="-" , color ="r", label ="Mars")

ax1.set_xlim([-1e12,1e12]) #setting limits of axis
ax1.set_ylim([-1e12,1e12])
ax1.set_xlabel(r"$x[m]$") # labelling axis
ax1.set_ylabel(r"$y[m]$")
plt.legend()
plt.show()
fig1.tight_layout()
fig1.savefig("Euler_Method_Trajectories.png", dpi=150) #save figure 


#Plotting oribits using Euler-Cromer Method

Data2 = np.load("simulation_results_euler_cromer.npy", allow_pickle=True) #open data from simulation using the euler cromer method

fig2, ax2 = plt.subplots()
#Plotting Earth Orbit with Euler-Cromer Method
xe = []
ye = []
for i in range (0,niter): 
    earth_position = Data2[i][2].position
    xe.append(earth_position[0])
    ye.append(earth_position[1])
ax2.plot(xe, ye, linestyle ="-" , color ="b", label ="Earth")

#Plotting Sun Orbit with Euler-Cromer Method
xs=[]
ys=[]
for i in range (0,niter):
    sun_position = Data2[i][1].position
    xs.append(sun_position[0])
    ys.append(sun_position[1])
ax2.plot(xs, ys,linestyle = "-", color="k", label="Sun")

#Plotting Jupiter Orbit with Euler-Cromer Method
xj =[]
yj =[]
for i in range (0,niter):
    jupiter_position = Data2[i][3].position
    xj.append(jupiter_position[0])
    yj.append(jupiter_position[1])
ax2.plot(xj, yj,linestyle="-",color="g", label="Jupiter")

#Plotting Mars Orbit with Euler-Cromer Method
xm = []
ym = []
for i in range (0,niter): 
    mars_position = Data2[i][4].position
    xm.append(mars_position[0])
    ym.append(mars_position[1])
ax2.plot(xm, ym, linestyle ="-" , color ="r", label ="Mars")

ax2.set_xlim([-1e12,1e12])
ax2.set_ylim([-1e12,1e12])
ax2.set_xlabel(r"$x[m]$")
ax2.set_ylabel(r"$y[m]$")
plt.legend()
plt.show()
fig2.tight_layout()
fig2.savefig("Euler_Cromer_Method_Trajectories.png", dpi=150)




