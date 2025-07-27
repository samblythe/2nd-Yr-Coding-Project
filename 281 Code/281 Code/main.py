import numpy as np
from Particle import Particle
import copy


import astropy
from poliastro import constants
from astropy.constants import G  # Newton's gravitational constant

# Sun mass (converting to kg)
msun = (constants.GM_sun / G).value

# Earth mass (converting to kg)
mearth = (constants.GM_earth / G).value 

# Jupiter mass (converting to kg)
mjupiter = (constants.GM_jupiter / G).value

# Mars mass (converting to kg)
mmars = (constants.GM_mars / G).value


from astropy.time import Time

# get the time at 12am on 14th Dec 2022 - this is the time for which all initial positions and velocities of boides are taken from the JPL ephemeris
t = Time("2022-12-14 00:00:00.0", scale="tdb")

from astropy.coordinates import get_body_barycentric_posvel

# get positions of velocities for the Sun
sun_pos, sun_vel = get_body_barycentric_posvel("sun", t, ephemeris="jpl")

from spiceypy import sxform, mxvg

# make a "state vector" of positions and velocities (in metres and metres/second, respectively)
sun_statevec = [
    sun_pos.xyz[0].to("m").value,
    sun_pos.xyz[1].to("m").value,
    sun_pos.xyz[2].to("m").value,
    sun_vel.xyz[0].to("m/s").value,
    sun_vel.xyz[1].to("m/s").value,
    sun_vel.xyz[2].to("m/s").value,
]

# get transformation matrix to the ecliptic (use time in Julian Days)
trans = sxform("J2000", "ECLIPJ2000", t.jd)

# transform state vector to ecliptic
sun_statevececl = mxvg(trans, sun_statevec)

# get positions and velocities
sun_position = [sun_statevececl[0], sun_statevececl[1], sun_statevececl[2]]
sun_velocity = [sun_statevececl[3], sun_statevececl[4], sun_statevececl[5]]


# get positions of velocities for the Earth
earth_pos, earth_vel = get_body_barycentric_posvel("earth", t, ephemeris="jpl")

from spiceypy import sxform, mxvg

# make a "state vector" of positions and velocities (in metres and metres/second, respectively)
earth_statevec = [
    earth_pos.xyz[0].to("m").value,
    earth_pos.xyz[1].to("m").value,
    earth_pos.xyz[2].to("m").value,
    earth_vel.xyz[0].to("m/s").value,
    earth_vel.xyz[1].to("m/s").value,
    earth_vel.xyz[2].to("m/s").value,
]

# get transformation matrix to the ecliptic (use time in Julian Days)
trans = sxform("J2000", "ECLIPJ2000", t.jd)

# transform state vector to ecliptic
earth_statevececl = mxvg(trans, earth_statevec)

# get positions and velocities
earth_position = [earth_statevececl[0], earth_statevececl[1], earth_statevececl[2]]
earth_velocity = [earth_statevececl[3], earth_statevececl[4], earth_statevececl[5]]

# get positions of velocities for Jupiter
jup_pos, jup_vel = get_body_barycentric_posvel("jupiter", t, ephemeris="jpl")

from spiceypy import sxform, mxvg

# make a "state vector" of positions and velocities (in metres and metres/second, respectively)
jup_statevec = [
    jup_pos.xyz[0].to("m").value,
    jup_pos.xyz[1].to("m").value,
    jup_pos.xyz[2].to("m").value,
    jup_vel.xyz[0].to("m/s").value,
    jup_vel.xyz[1].to("m/s").value,
    jup_vel.xyz[2].to("m/s").value,
]

# get transformation matrix to the ecliptic (use time in Julian Days)
trans = sxform("J2000", "ECLIPJ2000", t.jd)

# transform state vector to ecliptic
jup_statevececl = mxvg(trans, jup_statevec)

# get positions and velocities
jupiter_position = [jup_statevececl[0], jup_statevececl[1], jup_statevececl[2]]
jupiter_velocity = [jup_statevececl[3], jup_statevececl[4], jup_statevececl[5]]

#get positions of velocities for Jupiter
mar_pos, mar_vel = get_body_barycentric_posvel("mars", t, ephemeris="jpl")

from spiceypy import sxform, mxvg

# make a "state vector" of positions and velocities (in metres and metres/second, respectively)
mar_statevec = [
    mar_pos.xyz[0].to("m").value,
    mar_pos.xyz[1].to("m").value,
    mar_pos.xyz[2].to("m").value,
    mar_vel.xyz[0].to("m/s").value,
    mar_vel.xyz[1].to("m/s").value,
    mar_vel.xyz[2].to("m/s").value,
]

# get transformation matrix to the ecliptic (use time in Julian Days)
trans = sxform("J2000", "ECLIPJ2000", t.jd)

# transform state vector to ecliptic
mar_statevececl = mxvg(trans, mar_statevec)

# get positions and velocities
mars_position = [mar_statevececl[0], mar_statevececl[1], mar_statevececl[2]]
mars_velocity = [mar_statevececl[3], mar_statevececl[4], mar_statevececl[5]]

#get positions of velocities for Neptune
nep_pos, nep_vel = get_body_barycentric_posvel("neptune", t, ephemeris="jpl")

from spiceypy import sxform, mxvg

# make a "state vector" of positions and velocities (in metres and metres/second, respectively)
nep_statevec = [
    nep_pos.xyz[0].to("m").value,
    nep_pos.xyz[1].to("m").value,
    nep_pos.xyz[2].to("m").value,
    nep_vel.xyz[0].to("m/s").value,
    nep_vel.xyz[1].to("m/s").value,
    nep_vel.xyz[2].to("m/s").value,
]

# get transformation matrix to the ecliptic (use time in Julian Days)
trans = sxform("J2000", "ECLIPJ2000", t.jd)

# transform state vector to ecliptic
nep_statevececl = mxvg(trans, nep_statevec)

# get positions and velocities
neptune_position = [nep_statevececl[0], nep_statevececl[1], nep_statevececl[2]]
neptune_velocity = [nep_statevececl[3], nep_statevececl[4], nep_statevececl[5]]



#Initialising each body in the simulation using the Particle class and positions and velocities from the JPL ephemeris (see above):

Earth = Particle(
    position=np.array(earth_position),
    velocity=np.array(earth_velocity),
    acceleration=np.array([0, 0, 0]),
    name="Earth",
    mass=mearth
)

Sun = Particle(
    position=np.array(sun_position),
    velocity=np.array(sun_velocity),
    acceleration=np.array([0, 0, 0]),
    name="Sun",
    mass=msun
)

Jupiter = Particle(
    position=np.array(jupiter_position),
    velocity=np.array(jupiter_velocity),
    acceleration=np.array([0, 0, 0]),
    name="Jupiter",
    mass=mjupiter
)

Mars = Particle(
    position=np.array(mars_position),
    velocity=np.array(mars_velocity),
    acceleration=np.array([0, 0, 0]),
    name="Mars",
    mass=mmars
)

bodies =[ Earth, Sun, Jupiter, Mars]  # make a list of bodies in order to loop through them all when calculating accelerations and positions

#Euler Method:


DataMethod1 = [] # empty list to save data to by appending

deltaT= 3600 #time step 
tend = 31536000*10 #year x number of years

niter = int(tend/deltaT) #number of itterations

#Finding the total linear momentum of the system at the begining (time = 0)
p0 = np.array([0, 0, 0], dtype=float)
for body in bodies:
    p0 += body.linearMomentum()    
p0_mag = np.sqrt((p0[0])**2  + (p0[1])**2   +(p0[2])**2)
p = np.array([0, 0, 0], dtype=float)

#Finding the total angular momentum of the system at the begining (time = 0)
L0 = np.array([0,0,0], dtype = float)
for body in bodies:
    L0 += body.angularMomentum()
L0_mag = np.sqrt((L0[0])**2  + (L0[1])**2   +(L0[2])**2)
L = np.array([0, 0, 0], dtype=float)

#Finding total energy of the system at the begining (time = 0)
e0 = 0.0
gE0 = 0.0
kE0 = 0.0
for body in bodies:
    gE0 += body.gravPotentialEngery(bodies)
    kE0 += body.kineticEnergy()
e0 = gE0 + kE0
e = 0.0
DataMethod1 = [] # empty list to save data to by appending


for  i in range (0, niter): # this loop will append the DataMethod1 list with another list containing all the information of the simulation for that itteration, niter times
    method = 1 # euler method in use

    #initialse all values backt to zero for each itteration:
    p = np.array([0, 0, 0], dtype=float)
    L = np.array([0, 0, 0], dtype=float)
    e = 0.0
    gE = 0.0
    kE = 0.0

    time = deltaT*i + 6 #total time elapsed
    for body in bodies: # loop through all the bodies in the list and calculate their accelerations, velocities, positions angular Momentums, linear Momentums and energies
        
        body.updateGravitationalAcceleration(bodies)
        body.update(deltaT, method)
        p += body.linearMomentum() #adds on the momenutum of the body each time round the loop to get a total for all bodies at the end
        L += body.angularMomentum() # ^^^^ same process for angular momentum
        gE += body.gravPotentialEngery(bodies)
        kE += body.kineticEnergy()
    e = gE + kE0
    p_mag = np.sqrt((p[0])**2  + (p[1])**2   +(p[2])**2) # finding magnitude of P needed for calculations
    fracChangeP = np.sqrt(((p_mag - p0_mag)/p0_mag )**2) # finding fractional change usuing inital value at t=0 and current value. getting modulus using sqrt of squar

    L_mag = np.sqrt((L[0])**2  + (L[1])**2   +(L[2])**2) # ^^^ same process as above
    fracChangeL = np.sqrt(((L_mag - L0_mag)/L0_mag )**2) # 

    fracChangeE = np.sqrt(((e-e0)/e0)**2)

    #append the list DataMethod1 with a list of data from this itteration
    DataMethod1.append([time, copy.deepcopy(Sun), copy.deepcopy(Earth), copy.deepcopy(Jupiter), copy.deepcopy(Mars), fracChangeP, fracChangeL, fracChangeE]) 
np.save("simulation_results_euler", DataMethod1, allow_pickle=True) # save resutls to an npy file for use later in analysis



#Euler-Cromer Method:


#need to reinitialise all particles and values again beofore doing next method:


Earth = Particle(
    position=np.array(earth_position),
    velocity=np.array(earth_velocity),
    acceleration=np.array([0, 0, 0]),
    name="Earth",
    mass=mearth
)

Sun = Particle(
    position=np.array(sun_position),
    velocity=np.array(sun_velocity),
    acceleration=np.array([0, 0, 0]),
    name="Sun",
    mass=msun
)

Jupiter = Particle(
    position=np.array(jupiter_position),
    velocity=np.array(jupiter_velocity),
    acceleration=np.array([0, 0, 0]),
    name="Jupiter",
    mass=mjupiter
)

Mars = Particle(
    position=np.array(mars_position),
    velocity=np.array(mars_velocity),
    acceleration=np.array([0, 0, 0]),
    name="Mars",
    mass=mmars
)




bodies =[ Earth, Sun, Jupiter , Mars]



#Finding the total linear momentum of the system at the begining (time = 0)
p0 = np.array([0, 0, 0], dtype=float)
for body in bodies:
    p0 += body.linearMomentum()   
p0_mag = np.sqrt((p0[0])**2  + (p0[1])**2   +(p0[2])**2)
p = np.array([0, 0, 0], dtype=float)

#Finding the total angular momentum of the system at the begining (time = 0)
L0 = np.array([0,0,0], dtype = float)
for body in bodies:
    L0 += body.angularMomentum()
L0_mag = np.sqrt((L0[0])**2  + (L0[1])**2   +(L0[2])**2)
L = np.array([0, 0, 0], dtype=float)

#Finding total energy of the system at the begining (time = 0)
e0 = 0.0
gE0 = 0.0
kE0 = 0.0
for body in bodies:
    gE0 += body.gravPotentialEngery(bodies)
    kE0 += body.kineticEnergy()
e0 = gE0 + kE0

DataMethod2 = []

#same process as above but using method = 2 so the program runs the euler-cromer method:

for  i in range (0, niter):
    method = 2
    p = np.array([0, 0, 0], dtype=float)
    L = np.array([0, 0, 0], dtype=float)
    e = 0.0
    gE = 0.0
    kE = 0.0
    time = deltaT*i + 3600
    for body in bodies:
        
        body.updateGravitationalAcceleration(bodies)
        body.update(deltaT, method)
        p += body.linearMomentum()
        L += body.angularMomentum()
        gE += body.gravPotentialEngery(bodies)
        kE += kE + body.kineticEnergy()
    e = gE + kE0

    p_mag = np.sqrt((p[0])**2  + (p[1])**2   +(p[2])**2)
    fracChangeP = np.sqrt(((p_mag - p0_mag)/p0_mag )**2) # getting modulus using sqrt

    L_mag = np.sqrt((L[0])**2  + (L[1])**2   +(L[2])**2)
    fracChangeL = np.sqrt(((L_mag - L0_mag)/L0_mag )**2) # getting modulus using sqrt

    fracChangeE = np.sqrt(((e-e0)/e0)**2)

    DataMethod2.append([time, copy.deepcopy(Sun), copy.deepcopy(Earth), copy.deepcopy(Jupiter), copy.deepcopy(Mars), fracChangeP, fracChangeL, fracChangeE])

np.save("simulation_results_euler_cromer", DataMethod2, allow_pickle=True) 
#save results to a new npy file for analysis




print("done") # this program takes a while to run through all the itterations so its useful to have an output for when it has reached the end



