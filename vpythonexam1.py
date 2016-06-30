from __future__ import division
from visual import *
sphere()

ball = sphere (pos=vector(-5,0,0), radius=.5, color=color.cyan,make_trail=True,trail_type="points",
              interval=10, retain=50)

wallR = box(pos=(6,0,0), size=(0.2,12,12), color=color.green)
wallL = box(pos=(-6,0,0), size=(0.2,12,12), color=color.green)
wallU = box(pos=(0,6,0), size=(12,0.2,12), color=color.green)
wallD = box(pos=(0,-6,0), size=(12,0.2,12), color=color.green)
wallB = box(pos=(0,0,-6), size=(12,12,0.2), color=color.red)
#wallF = box(pos=(0,0,6), size=(12,12,0.2), color=color.green)
ball.velocity = vector(25,5,-5)

vscale = 0.1

#varr = arrow(pos=ball.pos, axis=vscale*ball.velocity, color=color.yellow, make_trail=True,trail_type="points",
#              interval=10, retain=50)
#ball.trail.append(pos=ball.pos)

deltat = 0.005
t = 0
scene.autoscale = False
while t <3e10:
    rate(100)
    if ball.pos.x > wallR.pos.x:
        ball.velocity.x = -ball.velocity.x
    ball.pos = ball.pos + ball.velocity*deltat
    #varr.pos = ball.pos + ball.velocity*deltat
    t =t + deltat
    if ball.pos.x < wallL.pos.x:
            ball.velocity.x = - ball.velocity.x
            #arrow.trail.append(ball.pos)

    if ball.pos.y > wallU.pos.y:
        ball.velocity.y = - ball.velocity.y
    if ball.pos.y < wallD.pos.y:
        ball.velocity.y = - ball.velocity.y
    if ball.pos.z < wallB.pos.z:
        ball.velocity.z = - ball.velocity.z

    #if ball.pos.z > wallU.pos.z:
    #    ball.velocity.z = - ball.velocity.z

    #if ball.pos.x < wallR.pos.x:
    #        ball.velocity.x = -ball.velocity.x
    #ball.pos = ball.pos + ball.velocity*deltat
    #ball.trail.append(pos=ball.pos)
    #varr.pos = ball.pos + ball.velocity*deltat
    #t =t + deltat
#import ROOT
