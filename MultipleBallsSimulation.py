#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 21:07:20 2021

@author: claraaldegundemanteca
"""

import pylab as pl  
import numpy as np
from BallClass import *

class SimulationMultipleBalls:
    '''
    Creates a simulation for a manually created list of balls. The container
    is one element of the list, with the characteristics given by the 
    BallClass (negative radius,large mass, 
    centered at the origin and not moving).
    '''
    
    def __init__(self, ball_list): 
        #Container is an element of the ball list
        self._time=0
        self._delta_p_total=0
        self._ball_list=ball_list
        self._collision_n=0
        for i in range (len(self._ball_list)):
            if self._ball_list[i].iscontainer()==True:
                self._container=self._ball_list[i]

    def next_collision(self):
        ''' 
        Moves the system to the next collision and performs it.
        '''
        self._collision_n+=1 
        #Useful to check if collisions are happening
        t_min=1e6 
        #Give a large time, re write it if t to collision is less
        for i in np.arange(0,len(self._ball_list)):
            for k in np.arange(i+1,len(self._ball_list)):
                t_collide=self._ball_list[i].\
                time_to_collision(self._ball_list[k]) 
                if t_collide!= None and t_collide > 0:
                    if t_collide < t_min:
                        t_min=t_collide
                        index1=i 
                        index2=k
                        #Indeces of balls that give minimum time to collision
        ball_colliding1=self._ball_list[index1]
        ball_colliding2=self._ball_list[index2]
        
        for ball in self._ball_list:
            ball.move(t_min)
            self._time+=t_min          
        ball_colliding1.collide(ball_colliding2)

    def run(self, num_frames, animate=False):
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
            ax.set_aspect('equal')
            ax.add_artist(self._container.get_patch())
            for i in range(0,(len(self._ball_list))):
                if self._ball_list[i]!=self._container:
                    ax.add_patch(self._ball_list[i].get_patch()) 
        for frame in range(num_frames): #Each frame is a collision
            self.next_collision()
            if animate:
                pl.pause(.1)
        if animate:
            pl.show()
            
