import pygame as pg
from pygame.locals import *
from sys import exit
from settings import *
import numpy as np
import Body

pg.init()
screen = pg.display.set_mode((WIDTH, HEIGHT))
pg.display.set_caption(caption)
clock = pg.time.Clock()

offset = [0, 0]
zoom = 1
lock = None

arrowKeys = {'up': False, 'down': False, 'left': False, 'right': False, 'zoomIn': False, 'zoomOut': False}
keyBindings = {K_UP:'up', K_DOWN:'down', K_LEFT:'left', K_RIGHT:'right', K_z:'zoomIn', K_x:'zoomOut'}
numberKeys = {K_0:0, K_1:1, K_2:2, K_3:3, K_4:4, K_5:5, K_6:6, K_7:7, K_8:8, K_9:9}

otherKeys = {K_w:'forward', K_s:'backward', K_a:'turnLeft', K_d:'turnRight'}
keysHolding = {'turnLeft': False, 'turnRight': False, 'forward': False, 'backward': False}

trail = []

def main():
    global offset, zoom, panspeed, lock
    bodies = []
    #The sun

    bodies = [Body.CelestialBody([WIDTH // 2, HEIGHT // 2], [0, 0], 50, 20, (255, 255, 0), acc=[0, 0], locked=True)]

    #The earth
    bodies.append(Body.CelestialBody([(WIDTH // 2+180), (HEIGHT // 2)], [0, 0], 10, 10, (0, 0, 255),
                                     acc=[0, 0], orbit=[0], gravity=[0], bodies=bodies))

    #The moon
    #bodies.append(Body.CelestialBody([(WIDTH // 2), (HEIGHT // 2)+180], [0, 0], 0.1, 5, (128, 128, 128),
    #                                 acc=[0, 0], orbit=[0], gravity=[0], bodies=bodies))

    #The spaceship
    bodies.append(Body.ManouverBody([(WIDTH // 2+240), (HEIGHT // 2)], [0, 0], 0.0001, 5, (255, 255, 255),
                                     acc=[0, 0], orbit=[1], gravity=[1], bodies=bodies, spritePath='spaceShip.png'))

    for b in bodies:
        print(b)
    maneouverBody = 2

    while True:
        for event in pg.event.get():
            if event.type == QUIT:
                pg.quit()
                exit()

            if event.type == KEYDOWN:
                if event.key in keyBindings.keys():
                    key = keyBindings[event.key]
                    arrowKeys[key] = True

                if event.key in numberKeys.keys():
                    number = numberKeys[event.key] -1
                    lock = number
                    if number == -1 or number >= len(bodies):
                        lock = None

                if event.key in otherKeys.keys():
                    key = otherKeys[event.key]
                    keysHolding[key] = True

            if event.type == KEYUP:
                if event.key in keyBindings.keys():
                    key = keyBindings[event.key]
                    arrowKeys[key] = False

                if event.key in otherKeys.keys():
                    key = otherKeys[event.key]
                    keysHolding[key] = False

        for key in arrowKeys.keys():
            if arrowKeys[key]:
                change = panspeed/zoom
                if key == 'up':
                    offset[1] += change
                elif key == 'down':
                    offset[1] -= change
                elif key == 'left':
                    offset[0] += change
                elif key == 'right':
                    offset[0] -= change

                elif key == 'zoomIn':
                    zoom += 0.1

                elif key == 'zoomOut':
                    if zoom > 0.1:
                        zoom -= 0.1

        for key in keysHolding.keys():
            if keysHolding[key]:
                if key == 'turnLeft':
                    bodies[maneouverBody].angle += turnSpeed
                elif key == 'turnRight':
                    bodies[maneouverBody].angle -= turnSpeed

                elif key == 'forward':
                    bodies[maneouverBody].vel[0] -= np.cos(bodies[maneouverBody].angle) * shipSpeed
                    bodies[maneouverBody].vel[1] += np.sin(bodies[maneouverBody].angle) * shipSpeed
                    #bodies[maneouverBody].updateStartingAngle(bodies)

        if lock != None:
            offset = -(bodies[lock].pos - np.array([WIDTH // 2, HEIGHT // 2]))

        for body in bodies:
            body.handleGravity(bodies)

        screen.fill(BG_COLOR)
        for star in stars:
            pg.draw.circle(screen, star[3], star[:2], star[2])

        for body in bodies:
            body.update(bodies)
            body.draw(screen, offset, zoom)

        trail.append([bodies[-1].pos[0], bodies[-1].pos[1]])

        for orbitTarget in range(len(bodies)):
            if orbitTarget != None and len(bodies[orbitTarget].orbits) == 1:

                for coord in Body.calcExpectedOrbit(bodies[orbitTarget], bodies[bodies[orbitTarget].orbits[0]]):
                    x, y = coord[0]+offset[0], coord[1]+offset[1]
                    x, y = Body.getZoomedPos(x, y, zoom)
                    pg.draw.circle(screen, bodies[orbitTarget].color, (x, y), 1)

        if drawTrace:
            for t in trail:
                x, y = t[0]+offset[0], t[1]+offset[1]
                x, y = Body.getZoomedPos(x, y, zoom)
                pg.draw.circle(screen, (255, 255, 255), (x, y), 1)

        pg.display.update()
        clock.tick(FPS)

if __name__ == '__main__':
    main()