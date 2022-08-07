import math
import numpy as np

def radius(a, e, v):
    return a * (1 - (e**2)) / (1 + (e * math.cos(v)))

def true_anomaly(M, e):
    return M + (((2 * e) - ((e**3) / 4)) * math.sin(M)) + (5 * (e**2) * math.sin(2 * M) / 4) + (13 * (e**3) * math.sin(3 * M) / 12)

def mean_anomaly(M0, t, t0, G, m, a):
    return M0 + ((t - t0) * math.sqrt(G * m / (a**3)))

def coordinates(rot, a, e, G, m, t, M0, t0):
    v = true_anomaly(mean_anomaly(M0, t, t0, G, m, a), e)
    r = radius(a, e, v)
    return np.array([math.cos(v + rot), math.sin(v + rot)]) * r