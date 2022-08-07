import pygame as pg
import numpy as np
from settings import *
import orbital_functions
import trymsLib as tl

def getMagnitude(vector):
    return np.sqrt(np.sum(vector ** 2))

def getNormalised(vector):
    return vector / getMagnitude(vector)

def calcOrbitVel(distance, mass, G):
    return np.sqrt(G * mass / distance)

def img(path, sizeX, sizeY):
    pic = pg.image.load(path).convert_alpha()
    pic = pg.transform.scale(pic, (sizeX, sizeY))
    return pic

def getZoomedPos(x, y, zoom):
    x = x - WIDTH // 2
    y = y - HEIGHT // 2
    x, y = x*zoom + WIDTH // 2, y*zoom + HEIGHT // 2
    return x, y

def getReducedMass(mass, mass2):
    return mass*mass2/(mass+mass2)

def getOrbitalEnergy(mass, velocity, distance):
    #use total mass here
    return 0.5 * velocity ** 2 - G * mass / distance

def getSemiMajorAxis(mass, energy):
    return -G * mass / (2 * energy)

def getSpesificAngularMomentum(relativeSpeed, relativePosition):
    relativeSpeed = np.array([relativeSpeed[0], relativeSpeed[1], 0])
    relativePosition = np.array([relativePosition[0], relativePosition[1], 0])

    return getMagnitude(np.cross(relativeSpeed, relativePosition))

def getEccentricity(totMass, relativeVelocity, relativePosition):
    relativeSpeed = getMagnitude(relativeVelocity)
    distance = getMagnitude(relativePosition)
    energy = getOrbitalEnergy(totMass, relativeSpeed, distance)
    angularMomentum = getSpesificAngularMomentum(relativeVelocity, relativePosition)
    gravParameter = G * totMass
    return np.sqrt(1 + 2*energy*(angularMomentum**2)/(gravParameter**2))

def calcOrbitInTime(time, massOrbitingBody, massOrbitedBody, posOrbitingBody, posOrbitedBody, velOrbitingBody, velOrbitedBody, angle):
    global M0
    totalMass = massOrbitedBody#+massOrbitingBody

    eccentricity = getEccentricity(totalMass, velOrbitingBody, posOrbitingBody-posOrbitedBody)

    relativeSpeed = getMagnitude(velOrbitingBody)
    distance = getMagnitude(posOrbitingBody-posOrbitedBody)

    energy = getOrbitalEnergy(totalMass, relativeSpeed, distance)

    a = getSemiMajorAxis(totalMass, energy) # semi-major axis
    if a < 0:
        return 0, 0
    M0 = 0 # constant, starting point of the orbit

    #rot
    rot = -angle
    pos = tl.coordinates(rot, a, eccentricity, G, totalMass, time, M0, 0) + posOrbitedBody
    return pos, 1

def calcExpectedOrbit(orbitingBody, orbitedBody):
    orbit = []
    relativeVelocity = orbitingBody.vel
    relativePosition = orbitingBody.pos - orbitedBody.pos
    orbitObject = orbital_functions.Orbit(relativePosition, relativeVelocity, G*orbitedBody.mass)
    for i in range(200): #Adjust this to do one full orbit
        #pos, inOrbit = calcOrbitInTime(i*25, orbitingBody.mass, orbitedBody.mass, orbitingBody.pos, orbitedBody.pos,
        #                               orbitingBody.vel, orbitedBody.vel, orbitingBody.startingAngleInOrbit)
        #if inOrbit:
        #    orbit.append(pos)
        orbitObject.time_update(i*30)
        orbitObject.update_state()
        orbit.append(orbitObject.position+orbitedBody.pos)

    return orbit

def getAngle(v1, v2):
    return np.arccos(np.dot(v1, v2) / (getMagnitude(v1) * getMagnitude(v2)))

def getSOIradius(orbitingBody, orbitedBody):
    relativeVelocity = orbitingBody.vel
    relativePosition = orbitingBody.pos - orbitedBody.pos

    orbitObject = orbital_functions.Orbit(relativePosition, relativeVelocity, G*orbitedBody.mass)
    a = orbitObject.semimajor_axis()
    return a*(orbitingBody.mass/orbitedBody.mass)**(2/5)

class CelestialBody():
    def __init__(self, pos, vel, mass, radius, color=(255, 255, 255), acc=(0.0, 0.0), locked=False, orbit=[], gravity=[], bodies=[]):
        self.pos = np.array([float(pos[0]), float(pos[1])])
        self.mass = mass
        self.vel = np.array([float(vel[0]), float(vel[1])])
        self.acc = np.array([float(acc[0]), float(acc[1])])
        self.color = color
        self.radius = radius
        self.locked = locked

        self.pos.astype(np.float64)
        self.vel.astype(np.float64)
        self.acc.astype(np.float64)

        self.gravity = gravity

        self.orbits = orbit

        if orbit != [] and not self.locked:
            #self.vel[0], self.vel[1] = 0, 0
            for o in orbit:
                o = bodies[o]
                mass = o.mass+self.mass
                pos = o.pos
                dist = self.pos - pos
                dist_mag = getMagnitude(dist)
                vel = calcOrbitVel(dist_mag, mass, G)

                #self.vel += vel * getNormalised(dist)
                self.vel[1] += -vel

        self.updateStartingAngle(bodies)

    def updateStartingAngle(self, bodies):
        ##set this when switching orbits
        if self.orbits != []:
            relPos = bodies[self.orbits[0]].pos - self.pos
            relVel = self.vel
            default = np.array([1, 0])

            self.startingAngleInOrbit = getAngle(relPos, default)
            if getMagnitude(relVel) > calcOrbitVel(getMagnitude(relPos), bodies[self.orbits[0]].mass, G):
                self.startingAngleInOrbit += np.pi

        else:
            self.startingAngleInOrbit = 0

    def update(self, bodies):
        self.vel += self.acc

        for o in self.orbits:
            self.pos += bodies[o].vel #The linked object's velocity is added to the current object's position

        self.pos += self.vel

    def handleGravity(self, bodies):
        if not self.locked:
            for body in self.gravity:
                body = bodies[body]
                dist = self.pos - body.pos
                dist_mag = getMagnitude(dist)
                if dist_mag > self.radius + body.radius:
                    dist_norm = getNormalised(dist)
                    acc = dist_norm * (G * body.mass / (dist_mag ** 2))
                    self.vel -= acc

    def draw(self, screen, offset, zoom):

        x, y = getZoomedPos((self.pos[0]+offset[0]), (self.pos[1]+offset[1]), zoom)
        pg.draw.circle(screen, (255, 255, 255), (x-1, y-1), self.radius*zoom)

        pg.draw.circle(screen, self.color, (x, y), self.radius*zoom)

    def __repr__(self):
        return 'CelestialBody(pos={}, vel={}, mass={}, radius={}, color={}, acc={}, locked={})'.format(self.pos, self.vel, self.mass, self.radius, self.color, self.acc, self.locked)

class ManouverBody(CelestialBody):
    def __init__(self, pos, vel, mass, radius, color=(255, 255, 255), acc=(0.0, 0.0), locked=False, orbit=[], gravity=[], bodies=[], spritePath=None):
        super().__init__(pos, vel, mass, radius, color, acc, locked, orbit, gravity, bodies)
        self.spritePath = spritePath
        self.angle = 0
        self.vel[1] -= 0.01
        self.updateStartingAngle(bodies)

    def baseUpdate(self, bodies):
        self.vel += self.acc

        for o in self.orbits:
            self.pos += bodies[o].vel  # The linked object's velocity is added to the current object's position

        self.pos += self.vel

    def update(self, bodies):
        distances = []
        if self.gravity != []:
            for i, body in enumerate(bodies):
                if body != self:
                    dist = self.pos - body.pos
                    dist = getMagnitude(dist)
                    dist_val = dist*(self.mass/body.mass)**(2/5)
                    #getSOIradius(self, body)
                    distances.append([dist_val, i])


            best = min(distances)
            bestIndex = distances[distances.index(best)][1]
            if bestIndex != self.gravity[0]:
                self.gravity = [bestIndex]
                self.orbits = [bestIndex]
                self.updateStartingAngle(bodies)

        self.baseUpdate(bodies)


    def draw(self, screen, offset, zoom):
        if self.spritePath != None:
            sprite = img(self.spritePath, int(self.radius*2*zoom), int(self.radius*2*zoom))
            sprite = pg.transform.rotate(sprite, self.angle*180/np.pi)


            x, y = getZoomedPos((self.pos[0]//1+offset[0]-self.radius), (self.pos[1]//1+offset[1]-self.radius), zoom)
            screen.blit(sprite, (x, y))

