
import math
import numpy as np 

class Particle:
    import math
    import numpy as np
    #initialising all class attributes
    def __init__(
    self,
    position=np.array([0, 0, 0], dtype=float),
    velocity=np.array([0, 0, 0], dtype=float),
    acceleration=np.array([0, -10, 0], dtype=float),
    name='Ball',
    mass=1.0,
    G=6.67408E-11
    ):
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.acceleration = np.array(acceleration, dtype=float)
        self.name = name
        self.mass = mass
        self.G = G
    #defining how to print
    def __str__(self):
        return "Particle: {0}, Mass: {1:.3e}, Position: {2}, Velocity: {3}, Acceleration: {4}".format(
        self.name, self.mass,self.position, self.velocity, self.acceleration
        )
    #Updating position and velocity using two different methods:
    def update(self, deltaT, method):
        if method == 1: # Euler Method
            self.position = self.position + self.velocity*deltaT
            self.velocity = self.velocity + self.acceleration*deltaT

        if method == 2: #Euler-Cromer Method
            self.velocity = self.velocity + self.acceleration*deltaT
            self.position = self.position + self.velocity*deltaT
            


    #this method updates the gravitational acceleration of a body by finding the acceleration caused by all the bodies interacting with it and adding them together
    def updateGravitationalAcceleration(self, bodies):
        self.acceleration = np.array([0,0,0],dtype=float) #initialise acceleration to 0
        for body in bodies: #loop through every body in the list bodies
            if body == self: # a body cant accelerate itself
                self.acceleration += np.array([0,0,0],dtype=float)
            else:
                dir = self.position - body.position #direction vector
                r2 = (dir[0])**2 + (dir[1])**2 + (dir[2])**2
                r = np.sqrt(r2) #magnitude of distance between bodies
                accel = -((self.G * body.mass) / r2 ) * (dir / r) #calculating acceleration
                self.acceleration += accel # add the calculated acceleration to the acceleration

        
    
        

    #self explainatory method for calculating the kinetic energy of a body:
    def kineticEnergy(self):
        velocityMagnitude = math.sqrt(((self.velocity[0])**2) + ((self.velocity[1])**2) + ((self.velocity[2])**2))
        kineticEnergy = 0.5*self.mass*(velocityMagnitude**2)
        return kineticEnergy
    
    #using a similar method to gravitational acceleration to find gravitational potential energy
    def gravPotentialEngery(self,bodies):
        gravPotentialEnergy = 0
        for body in bodies:
            if body == self:
                gravPotentialEnergy += 0
            else:
                dir = self.position - body.position
                r2 = (dir[0])**2 + (dir[1])**2 + (dir[2])**2
                r = np.sqrt(r2)
                gravPotentialEnergy += -((self.G*self.mass*body.mass)/r)
        return gravPotentialEnergy
    #linear momentum calculations:
    def linearMomentum(self):
        p_x = self.mass * self.velocity[0]
        p_y = self.mass * self.velocity[1]
        p_z = self.mass * self.velocity[2]
        linearMomentum_vector = np.array([p_x, p_y,p_z], dtype= float)
        return linearMomentum_vector
    #angular momentum calculations
    def angularMomentum(self):
        L_x = self.mass * self.velocity[0]
        L_y = self.mass * self.velocity[1]
        L_z = self.mass * self.velocity[2]
        linearMomentum_vector = np.array([L_x, L_y,L_z], dtype= float)

        angularMomentum_vector = np.array(np.cross(self.position, linearMomentum_vector), dtype= float)
        return angularMomentum_vector












