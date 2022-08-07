# -*- coding: utf-8 -*-

#This file was made by Trym Vegard Sæle Marton

import numpy as np
from scipy.spatial.transform import Rotation

"""
INSTRUCTIONS FOR USE:

    
In order to initialize an orbit,
assign a variable to the "Orbit" class, for example:

orbit1 = Orbit()

This will initialize an orbit, computing the necessary
orbital elements, etc.

The arguments are:
relative position, relative velocity
(both in 2D, relative to the parent object),
mu (G times mass of parent object),
and the current global timestamp.


In order to determine the position and velocity
in the orbit at a certain moment, first call the
.time_update() and then the .update_state() functions.

.time_update() takes two arguments:
the current global timestamp,
and the maximum error of the newton-raphson solver
(0.0001 is a value that works, i hardly believe that
any larger accuracy is necessary)

.update_state() takes no arguments.


When both of those functions are called, in that order,
relative position and velocity in the orbit can be called by the
.position and .velocity properties of the orbit object
"""

def state_vectors_to_orbital_elements(r_vec, v_vec, K_vec, mu):
    v_rad = np.dot(v_vec, (r_vec / np.linalg.norm(r_vec)))
    v_per = np.sqrt((np.linalg.norm(v_vec))**2 - v_rad**2)
    
    h_vec = np.cross(r_vec, v_vec)
    
    i = np.arccos(h_vec[2] / np.linalg.norm(h_vec))
    
    N_vec = np.cross(K_vec, h_vec)
    
    if N_vec[1] >= 0:
        asn = np.arccos(N_vec[0] / np.linalg.norm(N_vec))
    
    else:
        asn = (2 * np.pi) - np.arccos(N_vec[0] / np.linalg.norm(N_vec))
    
    e_vec = (np.cross(v_vec, h_vec) / mu) - (r_vec / np.linalg.norm(r_vec))
    
    if e_vec[2] >= 0:
        arp = np.arccos(np.dot(e_vec, N_vec) / (np.linalg.norm(e_vec) * np.linalg.norm(N_vec)))
    
    else:
        arp = (2 * np.pi) - np.arccos(np.dot(e_vec, N_vec) / (np.linalg.norm(e_vec) * np.linalg.norm(N_vec)))
    
    if v_rad >= 0:
        tra = np.arccos(np.dot(e_vec, r_vec) / (np.linalg.norm(e_vec) * np.linalg.norm(r_vec)))
    
    else:
        tra = (2 * np.pi) - np.arccos(np.dot(e_vec, r_vec) / (np.linalg.norm(e_vec) * np.linalg.norm(r_vec)))
    
    return np.linalg.norm(h_vec), np.linalg.norm(e_vec), i, asn, arp, tra

def orbital_elements_to_state_vectors(mu, tra, h, e, i, asn, arp):
    r_vec = h ** 2 / mu / (1 + e * np.cos(tra)) * np.array((np.cos(tra), np.sin(tra), 0))
    v_vec = mu / h * np.array((-np.sin(tra), e + np.cos(tra), 0))
    
    R = Rotation.from_euler("ZXZ", [-arp, -i, -asn])
    r_rot = r_vec @ R.as_matrix()
    v_rot = v_vec @ R.as_matrix()
    
    return r_rot, v_rot

def kepler_elliptical(m, e, error):
    E = m
    
    while True:
        if abs(E - (e * np.sin(E)) - m) < error:
            return E
        
        E = E - ((E - (e * np.sin(E)) - m) / (1 - (e * np.cos(E))))

def kepler_hyperbolic(m, e, error):
    #m = m * 1j
    H = m
    
    while True:
        if abs((e * np.sinh(H)) - H - m) < error:
            return H
        
        H = H - (((e * np.sinh(H)) - H - m) / ((e * np.cosh(H)) - 1))

class Orbit:
    def __init__(self, relative_position, relative_velocity, mu, t=0):
        self.mu = mu
        self.t0 = t
        
        K = np.array([0, 0, 1])
    
        dpos = np.array([relative_position[0], relative_position[1], 0])
        dvel = np.array([relative_velocity[0], relative_velocity[1], 0.0000000001])
    
        self.h_vec, self.e_vec, self.i, self.asn, self.arp, self.tra = state_vectors_to_orbital_elements(dpos, dvel, K, self.mu)
        
        e = np.linalg.norm(self.e_vec)
        
        if e < 1:
            E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(self.tra / 2))
            self.m0 = E - (e * np.sin(E))
            #self.m0 = np.arctan2(- np.sqrt(1 - (np.linalg.norm(self.e_vec))**2) * np.sin(self.tra), - np.linalg.norm(self.e_vec) - np.cos(self.tra)) + np.pi + (np.linalg.norm(self.e_vec) * np.sqrt(1 - (np.linalg.norm(self.e_vec))**2) * np.sin(self.tra) / (1 + (np.linalg.norm(self.e_vec) * np.cos(self.tra))))
        
        if e > 1:
            H = 2 * np.arctanh(np.sqrt((e - 1) / (e + 1)) * np.tan(self.tra / 2))
            self.m0 = (e * np.sinh(H)) - H
        
        self.m = self.m0
        
    def update_state(self):
        position3D, velocity3D = orbital_elements_to_state_vectors(self.mu, self.tra, np.linalg.norm(self.h_vec), np.linalg.norm(self.e_vec), self.i, self.asn, self.arp)
        self.position = np.array([position3D[0], position3D[1]])
        self.velocity = np.array([velocity3D[0], velocity3D[1]])

    def semimajor_axis(self):
        e = np.linalg.norm(self.e_vec)
        a = abs((np.linalg.norm(self.h_vec)**2) / (self.mu * (1 - e**2)))
        return a
    
    def time_update(self, time, error=0.0001):
        e = np.linalg.norm(self.e_vec)
        a = abs((np.linalg.norm(self.h_vec)**2) / (self.mu * (1 - e**2)))
        self.m = self.m0 + (np.sqrt(self.mu / (a**3)) * (time - self.t0))
        
        if e < 1:
            E = kepler_elliptical(self.m, e, error)
            beta = e / (1 + np.sqrt(1 - (e**2)))
            self.tra = E + (2 * np.arctan((beta * np.sin(E)) / (1 - (beta * np.cos(E)))))
        
        if e >= 1:
            H = kepler_hyperbolic(self.m, e, error)
            self.tra = 2 * np.arctan(np.sqrt((e + 1) / (e - 1)) * np.tanh(H / 2))
            #http://control.asu.edu/Classes/MAE462/462Lecture05.pdf
            #https://space.stackexchange.com/questions/27602/what-is-hyperbolic-eccentric-anomaly-f