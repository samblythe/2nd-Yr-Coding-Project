main.py  main simulation file where instances of the class Particle are created using initial values from the JPL ephemeris. 
Bodies are advanced, conserved values are calculated and results are saved as npy files. 
Requires numpy, copy, astropy and "pip install jplephem spiceypy poliastro" in the terminal. 
Tested on Python 3.9.12.
No user inputs required.

particle.py  Particle class containing methods to update the gravitational acceleration and advance the system. 
Also contains methods for calculating the Linear Momentum, Angular Momentum and Energy.
Requires numpy.

analysis.py  Post-processing script to plot the orbits using data from simulations using Euler and Euler-Cromer methods.
Requires matplotlib and numpy

analysis2.py  Post-processing script to plot the fractional change in angular momentum to compare Euler and Euler cromer methods. 
Requires matplotlib and numpy

analysis3.py  Post-processing script to plot the fractional change in linear momentum to compare Euler and Euler cromer methods. 
Requires matplotlib and numpy

analysis4.py  Post-processing script to plot the fractional change in total energy to compare Euler and Euler cromer methods. 
Requires matplotlib and numpy