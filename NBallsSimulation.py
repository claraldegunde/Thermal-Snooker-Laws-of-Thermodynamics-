#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Nov 29 14:46:16 2021
 
@author: claraaldegundemanteca
"""
 
import pylab as pl 
import matplotlib.pyplot as plt
import numpy as np
from BallClass import *
 
time=[]
class SimulationNBalls:
    def __init__(self, n_of_balls, m_balls, rad_balls, rad_container, k=1):
        '''
        Initialise the class for a number of balls, their mass, radius and 
        radius of the container. Therefore, only works for  balls with same 
        mass and radius. 
        Note that the container is an element of the ball list.
        
        The factor k multiplies the velocity, which will change the KE and 
        therefore the temeprature.
        
        Tha balls are initialised in a square inside the container, with
        random velocities that average to 0.
        
        This is the same code as SimulationNBalls but slightly modified.
        '''
        self._time=0
        self._n_of_balls=n_of_balls
        self._delta_p_total=0
        self._ball_list=[]
        self._collision_n=0
        self._m_balls=m_balls
        self._rad_balls=rad_balls
        self._rad_container=rad_container

        #Create array of possible postitions:
        max_num_balls=int(((rad_container*(2)**0.5+0.5*rad_balls)/\
                           rad_balls)/2.5)
        pos_list=[]
        for n in range (0,max_num_balls):
            y=((-rad_container*(2)**0.5/2)+rad_balls)+n*(2.5*(rad_balls))
            for n in range (0,max_num_balls):
                x=((-rad_container*(2)**0.5/2)+rad_balls)+n*(2.5*(rad_balls))
                pos_list.append([x,y])
        possible_positions=np.array(pos_list) 
        
        #Create the ball list:
        container=Ball(iscontainer=True, rad=rad_container)
        self._ball_list.append(container)
        for i in range (self._n_of_balls):
            pos_pos=possible_positions[i,:]
            rand_v=[k*(np.random.normal(0, 20)), k*(np.random.normal(0, 20))]
            #So the average velocity is 0
            ball=Ball(m=m_balls, rad=rad_balls,pos=pos_pos, v=rand_v)           
            self._ball_list.append(ball)                                   
        for i in range (len(self._ball_list)):
            if self._ball_list[i].iscontainer()==True:
                self._container=self._ball_list[i]
   
    def number_of_collisions(self):
        '''
        Returns the total number of collisions ocurred in the simulation
        '''
        n_of_collisions=0
        for i in np.arange(0,len(self._ball_list)):
           n_of_collisions+=self._ball_list[i].number_of_collisions()
        return n_of_collisions
       
    def next_collision(self):
        '''
        Moves everything to the collision point and performs it
        '''
        #Calculate the time to the next collision
        t_min=1e6
        for i in np.arange(0,len(self._ball_list)):
            for k in np.arange(i+1,len(self._ball_list)):
                t_collide=self._ball_list[i].time_to_collision(self._ball_list[k])            
                if t_collide!= None and t_collide >0:
                    if t_collide < t_min:
                        t_min=t_collide
                        index_1=i
                        index_2=k    
        self._time+=t_min
        time.append(self._time) 
        ball_colliding1=self._ball_list[index_1]
        ball_colliding2=self._ball_list[index_2]
        
        for ball in self._ball_list:
            ball.move(t_min)
        ball_colliding1.collide(ball_colliding2)
        
        #Add changes of momentum to calculate pressure
        if ball_colliding1 == self._container:
            parall_dir=ball_colliding2.pos()/\
                (np.dot(ball_colliding2.pos(),ball_colliding2.pos()))**0.5
            v_parall=np.dot(ball_colliding2.v(),parall_dir)*parall_dir
            self._delta_p_total += 2*ball_colliding2.m()*\
                (np.dot(v_parall,v_parall))**0.5
                
        if ball_colliding2 == self._container:
            parall_dir=ball_colliding1.pos()/\
                (np.dot(ball_colliding1.pos(),ball_colliding1.pos()))**0.5
            v_parall=np.dot(ball_colliding1.v(),parall_dir)*parall_dir
            self._delta_p_total += 2*ball_colliding1.m()*\
                (np.dot(v_parall,v_parall))**0.5
   
    def run(self, num_frames, animate=False):
        '''
        Runs the simulation for a given number of frames (collisions). Default 
        animate is False.
        '''
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
            ax.set_aspect('equal')
            ax.add_artist(self._container.get_patch())
            for i in range(0,(len(self._ball_list))):
                if self._ball_list[i]!=self._container:
                    ax.add_patch(self._ball_list[i].get_patch())       
        for frame in range(num_frames):
            self.next_collision()
            if animate:
                pl.pause(0.1)
        if animate:
            pl.show()   
        
    def distance_to_centre(self, n_of_collisions, n_of_bins): 
        '''
        Plots an histogram of the distance from the centre for every ball,
        after a number of collisions.
        '''
        distances_centre=[]
        for i in range (n_of_collisions):
            print(i)
            self.next_collision()
            for i in range (0,len(self._ball_list)):
                if self._ball_list[i].iscontainer() == False:
                    distance_centre=(np.dot(self._ball_list[i].pos(),\
                                           self._ball_list[i].pos()))**0.5
                    distances_centre.append(distance_centre) 
        distances_centre=np.array(distances_centre)
        
        for i in range (len(distances_centre)):
            if distances_centre[i]>(self._rad_container-self._rad_balls):
                raise Exception ('A ball is outside the container')
        
        plt.figure()
        plt.grid()
        plt.hist(distances_centre,bins=n_of_bins, color='grey')
        plt.xlabel('Distance to centre (m)')
        plt.title('Histogram of distances to the centre')
        plt.ylabel('Number of particles')
        plt.xlim(-1,10)
        plt.xticks(np.arange(-1,10,step=1))
               
    def distance_between_balls(self, n_of_collisions, n_of_bins): 
        '''
        Plots an histogram of the distance from the centre for every ball,
        after a number of collisions.
        '''   
        distances_ball=[]
        for i in range (n_of_collisions):
            print(i)
            self.next_collision()
            for i in np.arange(0,len(self._ball_list)):
                if self._ball_list[i].iscontainer()==False:
                    for j in np.arange(i+1,len(self._ball_list)):
                        if self._ball_list[j].iscontainer()== False:
                            rel_pos=self._ball_list[i].pos()-\
                                self._ball_list[j].pos()
                            distance_balls=(np.dot(rel_pos,rel_pos))**0.5
                            distances_ball.append(distance_balls)
        distances_balls=np.array(distances_ball)
        
        plt.figure()
        plt.grid()
        plt.hist(distances_balls,bins=n_of_bins, color='grey')
        plt.xlabel('Distance between particles (m)')
        plt.title('Histogram of distances between particles')
        plt.ylabel('Number of particles')
        plt.xlim(-1,20)
        plt.xticks(np.arange(-1,20,step=1))
        
        for i in range (len(distances_ball)):
            if distances_ball[i]>2*(self._rad_container-self._rad_balls):
                raise Exception ('A ball is outside the container')       
   
    def total_momentum(self, n_of_collisions): 
        '''
        Momentum of all the particles after a given number of collisions
        Should stay constant regardless of n_of_collisions.
        '''
        total_momentum=0
        self.run(n_of_collisions)
        for i in np.arange(0,len(self._ball_list)):
            v=self._ball_list[i].v()
            p=self._ball_list[i].m()*v
            total_momentum += p
        return total_momentum    
    
    def total_mass(self):
        '''
        Total mass of all the balls inside the container
        '''
        return self._n_of_balls*self._m_balls
            
    def total_KE(self, n_of_collisions):
        '''
        Total kinetic energy of the system. Should stay constant.
        '''
        total_KE=0
        self.run(n_of_collisions)
        for i in np.arange(0,len(self._ball_list)):
            if self._ball_list[i].iscontainer() == False:
                v=(np.dot(self._ball_list[i].v(),self._ball_list[i].v()))**0.5
                total_KE += 0.5*(self._ball_list[i].m()*v*v)
        return total_KE      
    
    def temperature(self, n_of_collisions): 
        '''
        Temperature of the system, related to the average kinetic energy of 
        the particles.
        Should stay constant after a large number of collisions.
        '''
        kb=1.38e-23
        T=(self.total_KE(n_of_collisions)/self._n_of_balls)/kb
        #Dropped the 3/2 factor because we work in 2D
        return T
       
    def pressure(self, n_of_collisions):
        '''
        Returns the pressure of the system, calcuated using the total change 
        of momentum of the balls colliding with the container divided by both
        the time elapsed and the perimeter of the container. 
        Notice the force  resemblance.
        '''
        self.run(n_of_collisions) 
        total_time=self._time
        circumference=2*np.pi*self._rad_container
        p=self._delta_p_total/(total_time*circumference)
        return p
    
    def time_list(self):
        '''
        Returns a list of all the times to the next collision.
        '''
        return time
    
    def KE_conservation(self, n_of_collisions):
        '''
        Checks the conservation of kinetic energy, plotting its value against
        time up to a certain number of collisions. 
        Should get a straight horizontal line.
        '''
        self.run(n_of_collisions)
        KE=[]
        times=time[0:(n_of_collisions)]
        for i in range(n_of_collisions):
            KE.append(np.round(self.total_KE(i),4)) #To avoid rounding errors
        plt.figure()
        plt.grid()
        plt.plot(times,KE, color='#CC7722')
        plt.xlabel('Time (s)')
        plt.ylabel('Kinetic energy (J)')
        plt.title('Conservation of KE')
        
    def momentum_conservation(self, n_of_collisions): 
        '''
        Checks the conservation of momentum, plotting its value against
        time up to a certain number of collisions. 
        Should get a straight horizontal line.
        '''
        self.run(n_of_collisions)
        momentum=[]
        times=time[0:(n_of_collisions)]
        for i in range(n_of_collisions):
            p_vector=self.total_momentum(i)
            p=(np.dot(p_vector,p_vector))**0.5
            momentum.append(p) #To avoid rounding errors
        plt.figure()
        plt.grid()
        plt.plot(times,momentum,color='#CC7722')
        plt.xlabel('Time (s)')
        plt.ylabel('Momentum (kg m s$^{-1}$)')
        plt.title('Conservation of momentum')
 
    def velocities (self): 
        '''
        Returns a list of the velocities of all particles, used to get the 
        velocity distribution (Maxwell Boltzmann),
        '''
        velocities=[]
        for i in np.arange(0,len(self._ball_list)):
            if self._ball_list[i].iscontainer()==False:
                v=np.linalg.norm(self._ball_list[i].v())
                velocities.append(v)
        return velocities

