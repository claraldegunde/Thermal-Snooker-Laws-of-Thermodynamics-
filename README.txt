#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 13:36:23 2021

@author: claraaldegundemanteca
"""
'''
My code is formed by 4 modules and 2 scripts:
    
The modules:
    
    - BallClass.py is the class used to create the balls (also the container)
    for the simulation. It contains their characterisics, as well as the 
    information needed to perform the collisions (time to collision, move and 
    collide methods). It will be the base for all later simulations.
    
    - SingleBallSimulation.py creates the simulation for one single ball and a 
    container. This was used as a first step to develop my code, 
    to make testing easier and check that the ball was bouncing as expected.
    
    - MultimpleBallsSimulation.py extends the Single Ball Simulation code to a
    list of balls. This was done to check the behaviour of multiple balls
    before creating a simulation for n balls. Notice that the container is 
    one element of the ball list.
    
    - NBallsSimulation.py is the simulation for n balls, the one that was
    used for the investigations. It is in this module where we can find
    methods for temperature, pressure, total kinetic energy...
    
The scripts:

    - Tests: used to check the written code (velocities being changed
    correctly after each collision, time between collisions calculated as
    expected...).
    
    - Investigations: using the results of our simiulations, questions on 
    tasks 10-14 were answered here.

'''