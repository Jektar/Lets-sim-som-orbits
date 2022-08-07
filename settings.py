FPS = 30

caption = "Let's Sim Some Orbits"

WIDTH, HEIGHT = 640, 640

drawTrace = False

BG_COLOR = (0, 0, 0)

G = 1

shipSpeed = 0.005

turnSpeed = 0.1

panspeed = 8

nStars = (WIDTH * HEIGHT) // (3*10 ** 3)

import random, getStarColours
nStars = random.randint(nStars // 2, nStars)
colours = getStarColours.getStarColours()
stars = [[random.randint(0, WIDTH), random.randint(0, HEIGHT),
          random.randint(1, 1), random.choice(colours)] for i in range(nStars)]

