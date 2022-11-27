#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 18:56:09 2021

@author: claraaldegundemanteca
"""

import pylab as pl  
from BallClass import *
               
class SimulationSingleball:
    '''
    Creates a simulation for a single ball and a container, both defined 
    using BallClass.
    '''
    def __init__(self, ball, container): 
        self._ball=ball
        self._container=container

    def next_collision(self):
        ''' 
        Moves the system to the next collision and performs it.
        Pressure will be calculated in the NBallsSimulation module.
        '''
        t_collide=self._ball.time_to_collision(self._container)
        self._ball.move(t_collide)
        self._ball.collide(self._container)        
        
    def run(self, num_frames, animate=False):
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
            ax.add_artist(self._container.get_patch())
            ax.add_patch(self._ball.get_patch())
            ax.set_aspect('equal')
        for frame in range(num_frames):
            self.next_collision()
            if animate:
                pl.pause(0.001)
        if animate:
            pl.show()
