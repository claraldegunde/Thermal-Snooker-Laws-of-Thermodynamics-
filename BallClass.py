#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 20:10:29 2021

@author: claraaldegundemanteca
"""
import numpy as np
import pylab as pl   

class Ball:  
    '''
    Creates a ball given the mass, radius, position and velocity.
    Need to specify if it will be a container, so its mass, 
    velocity and position are defined (1e10, [0,0] and [0,0] respectively).
    The collor, filling of the patch can also be changed. 
    '''
    def __init__(self, m=0., rad=0., pos=[0,0], v=[0,0],iscontainer=False,\
                 color='b', fill=False,):
        self._m=m
        self._rad=rad
        self._pos=np.asarray(pos)
        self._v=np.asarray(v)
        self._color=color
        self._fill=fill
        self._iscontainer=iscontainer
        self._n=0 
        #Collision counter
        if iscontainer == True:
            self._m=1e10 
            self._fill=False
            self._color='black'
            self._rad=-rad 
        self._patch=pl.Circle([self._pos[0], self._pos[1]], self._rad,\
                              color=self._color,fill=self._fill) 
        if len(pos)>2:
            raise Exception("Position should be 2 dimensional")
            
        if len(v)>2:
            raise Exception("Velocity should be 2 dimensional")
            
    def pos(self):
        return  self._pos

    def v(self):
        return  self._v

    def set_v(self, v):
        self._v = v

    def m(self):
        return self._m
    
    def rad(self):
        return self._rad
    
    def iscontainer(self):
        return self._iscontainer

    def get_patch (self):
        return self._patch
    
    def number_of_collisions(self):
        return self._n
    
    def time_to_collision(self,other): #Calculated using equation on script
        R=self._rad+other.rad()
        rel_pos=self._pos-other._pos
        rel_v=self._v-other._v        
        v2=np.dot(rel_v,rel_v)
        r2=np.dot(rel_pos,rel_pos)
        rv=np.dot(rel_v,rel_pos)
        discr=4*rv*rv-4*v2*(r2-(R*R))
        if discr < 0:
            return None
        else:
            delta_t_minus=(-2*rv-discr**0.5)/(2*v2)
            delta_t_plus=(-2*rv+discr**0.5)/(2*v2)
            if v2 == 0:
                return None
            if (self.iscontainer()!=other.iscontainer()):
                return delta_t_plus
            if rv > 0:
                return None
            else:
                return delta_t_minus

    def move(self, dt):
        '''
        Moves the system to the next collision, after a time dt.
        '''
        if self._iscontainer==False: #Only moves if it is not a container
           self._pos=self._pos+dt* .9999999*self._v
        self._patch.center=self._pos #Moves the patch
        
    def collide (self, other): #0nly parallel components change
        '''
        Changes the velocities of the balls after the collision.
        '''
        rel_pos=self._pos-other.pos()
        parall_dir=rel_pos/(np.dot(rel_pos,rel_pos))**0.5 #Normalised
        v1_i_parall=np.dot(self._v,parall_dir)*parall_dir 
        v2_i_parall=np.dot(other._v,parall_dir)*parall_dir
        v1_perp=self._v-v1_i_parall 
        v2_perp=other._v-v2_i_parall #Constant, doesnt change
        v1_f_parall=((self._m -other.m())*v1_i_parall/(self._m +other.m())+\
                     (2*other.m()*v2_i_parall)/(self.m() +other.m()))  
        v2_f_parall=(((other.m() -self._m)*v2_i_parall)/(self._m +other.m())\
                     +2*self._m*v1_i_parall/(self._m +other.m()))
        self._v=v1_f_parall+v1_perp 
        other.set_v(v2_f_parall+v2_perp)   #Changes velocities
        
        self._n+=1
        other._n+=1 #Add 1 to the collision counter for both balls
        
        p_i=self._m*v1_i_parall+other.m()*v2_i_parall
        p_f=self._m*v1_f_parall+other.m()*v2_f_parall
        ke_i=0.5*(self._m*np.dot(v1_i_parall,v1_i_parall)+\
                  other.m()*np.dot(v2_i_parall,v2_i_parall))
        ke_f=0.5*(self._m*np.dot(v1_f_parall,v1_f_parall)+\
              other.m()*np.dot(v2_f_parall,v2_f_parall))
        if np.round(p_i[0],4) != np.round(p_f[0],4):
            raise Exception  ('Momentum not conserved') 
            #Up to 4 d.p. to avoid rounding
        if np.round(p_i[1],4) != np.round(p_f[1],4):
            raise Exception  ('Momentum not conserved')
        if np.round(ke_i,4) != np.round(ke_f,4):
            raise Exception  ('KE not conserved')
        return self._v,other.v()











                
               
              
                
                
                