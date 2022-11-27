#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 18:49:27 2021

@author: claraaldegundemanteca
"""
from BallClass import *
from SingleBallSimulation import *
from MultipleBallsSimulation import *
from NBallsSimulation import *
from scipy.optimize import curve_fit
# import pylab as pl 
import matplotlib.pyplot as plt
import numpy as np
                       
#%% TASK 10 (a)(b)

sim10=SimulationNBalls(10,1e-27,0.1,10,k=5)   
#Mass order of magnitude of a hydrogen atom

collisions=500 #Large number of collisions to reach stable values of T,P
temperature=[]
KE=[]
time=sim10.time_list()

for i in range(1,collisions):
    print(i)
    temperature.append(np.round(sim10.temperature(i),4))
    KE.append(sim10.total_KE(i)) 

plt.figure()
plt.grid()
plt.plot(time[:499], temperature) #Only want the first 500 times
plt.plot(time[:499], KE)
plt.xlabel('Time')
plt.title('KE and temperature relation')
plt.legend(['Temperature','Kinetic energy'])

'''
Expected result, T is a KE multiplied by a constant, and both should be 
constant in time. 

Momentum should also be conserved.

Time other than when 2 balls collide being considered in the graph.
'''
#%%TASK 10 (c)

collisions=100  
sim_vx1=SimulationNBalls(20,1e-27,0.1,10,k=1)   
sim_vx2=SimulationNBalls(20,1e-27,0.1,10,k=2)   #Doubled velocities here
temperature1=[]
temperature2=[]
pressure1=[]
pressure2=[]
time1=sim_vx1.time_list()
time2=sim_vx2.time_list()

for i in range(1,collisions):
    print(i)
    temperature1.append(np.round(sim_vx1.temperature(i),4))
    temperature2.append(np.round(sim_vx2.temperature(i),4))
    pressure1.append(sim_vx1.pressure(i))
    pressure2.append(sim_vx2.pressure(i))

                    #Temperature changing with doubling velocities
plt.figure()
plt.grid()
plt.plot(time1[10:(collisions-1)], temperature1[10:],'x') 
#Slicing to get stable values
plt.plot(time2[10:(collisions-1)], temperature2[10:],'x')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('Temperature dependence on initial velocities of the balls')
plt.legend(['Velocity x1','Velocity x2'])

                    #Pressure changing with doubling velocities
plt.figure()
plt.plot(time1[10:(collisions-1)], pressure1[10:],'x') 
#As before, we only want times for the 500 collisions
plt.plot(time2[10:(collisions-1)], pressure2[10:],'x')
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Pa)')
plt.title('Pressure dependence on initial velocities of the balls')
plt.legend(['Velocity x1','Velocity x2'])

'''

As expected the temperature gets multiplied by around a factor of 4, and so 
does pressure.

Check animation, balls move faster.
'''

#%% Task 11: (a) 
     
                 #Testing conservation of KE
                 
sim11=SimulationNBalls(20,1,0.1,10)        
print(sim11.total_KE(5))
print(sim11.total_KE(10))
print(sim11.total_KE(45)) #Get the same value no matter the number of frames

sim11.KE_conservation(20) #Doesnt change with n of collisions, seeen in graph

                # Conservation of momentum 
                
print(sim11.total_momentum(14))
print(sim11.total_momentum(10))
print(sim11.total_momentum(20))
print(sim11.total_momentum(30)) # Checking it manually


sim11.momentum_conservation(20) # Straight line so correct
                
#%% Task 11:  (b) 
 
'''
We can change T changing the velocity of the particles, via the k parameter.
Pressure should change by roglhy the same factor (see ideal gas equation).

Checking it manually gives us that P does depend on T. This
will be observed in more detail in the P against T graph.
'''
sim11b=SimulationNBalls(20,1e-27,0.1,10,k=2)  

print('Changing the velocity by a factor of 4 give us: P initial=',\
      sim11.pressure(100),'Pa and P final=',sim11b.pressure(100),' Pa')
    
print('So the pressure increased by roughly the same factor')
        
#%% Task 11:  (c) 

                        # Changing radius of container
                       
collisions=500
number_of_radius= 20

P=[]
volume=[] #Area of container
T=[]
for i in range (number_of_radius):
    print(i)
    sim11c=SimulationNBalls(50,1e-27,0.1,10+10*i) 
    pressure=sim11c.pressure(collisions)
    temperature=np.round(sim11c.temperature(collisions),4)
    P.append(pressure)
    T.append(temperature)
    volume.append(np.pi*(10+10*i)*(10+10*i))

def P_V (x, A):
    return 1.3*50*1.38e-23*T[5]*(1/x)

x=np.linspace(volume[0],volume[number_of_radius-1],100)
fit, cov=curve_fit(P_V, volume, P) 
data_fit =  P_V (x, *fit)
plt.figure()
plt.plot(volume,P,'x', color='black')
plt.plot(x, data_fit,color='#CC7722', linestyle='dashed')
plt.grid()
plt.ylabel('Pressure (Pa)')
plt.xlabel('Volume (m$^{3}$)')
plt.title('Pressure dependence on volume of container')
plt.legend(['Results from the simulation', 'Fit'])

'''
P chould decrease as we increase the volume of the container, as the 
frequency at which particles hit the walls decreases. 
'''

# plt.figure()
# plt.plot(volume,T,'x')
# plt.grid()
# plt.ylabel('Temperature (K)')
# plt.xlabel('Volume (m$^{3}$)')
# plt.title('Temperature dependence on volume of container')

'''
T shouldn't depend on the radius of the container, a it is a function 
of the average KE of the system, which only depends on the random velocities 
we inititialize the system with.

The fluctuations are due to the random nature of the velocities, but observe 
it takes roughly the same value regardless of the volume.
'''
                    #Changing number of balls
                    
collisions=500
number_of_balls= 10

P=[]
N=[] #List of number of balls

for i in range (number_of_balls):
    print(i)
    sim11c=SimulationNBalls(10+10*i,1e-27,0.1,10) 
    pressure=sim11c.pressure(collisions)
    P.append(pressure)
    N.append(10+10*i)
  
    
# plt.figure()
# fit,cov=np.polyfit(T,P, 1, cov=True)
# pfit=np.poly1d(fit)
# plt.grid()
# plt.plot(T,P,'x')
# plt.plot(T,pfit(T1))

fit, cov=np.polyfit(N,P, 1, cov=True)
pfit=np.poly1d(fit)
plt.figure()
plt.plot(N,P,'x', color='black')
plt.plot(N, pfit(N), color='#CC7722', linestyle='dashed')
plt.grid()
plt.ylabel('Pressure (Pa)')
plt.xlabel('Number of particles')
plt.title('Pressure dependence on number of particles')
plt.legend(['Results from the simulation', 'Fit'])
print('The slope of the graph is:', fit[0],'+/-', cov[0,0]**0.5)
                        

    # Changing number of balls(2), observing the stable P  and T values
                    
# As above, for a given number of collisions after which T and P stabilize 
# This is because for small t not all the balls have reached the container

collisions=200

sim_ballsx1=SimulationNBalls(10,1e-27,0.1,10)   
sim_ballsx5=SimulationNBalls(50,1e-27,0.1,10)
# sim_ballsx10=SimulationNBalls(100,0.1,0.1,10) 

temperature1=[]
temperature2=[]
# temperature3=[]
pressure1=[]
pressure2=[]
# pressure3=[]
time1=sim_ballsx1.time_list()
time2=sim_ballsx5.time_list()
# time3=sim_ballsx10.time_list()

for i in range(1,collisions):
    temperature1.append(np.round(sim_ballsx1.temperature(i),4))
    temperature2.append(np.round(sim_ballsx5.temperature(i),4))
    # temperature3.append(np.round(sim_ballsx10.temperature(i),4))
    pressure1.append(sim_ballsx1.pressure(i),4)
    pressure2.append(sim_ballsx5.pressure(i),4)
    # pressure3.append(np.round(sim_ballsx10.pressure(i),4))

# plt.figure()
# plt.grid()
# plt.plot(time1[100:(collisions-1)], temperature1[100:],'x')
# plt.plot(time2[100:(collisions-1)], temperature2[100:],'x')
# plt.xlabel('Time (s)')
# plt.ylabel('Temperature (K)')
# plt.title('Temperature dependence on number of balls of the balls')
# plt.legend(['10 balls','50 balls'])

plt.figure()
plt.grid()
plt.plot(time1[100:(collisions-1)], pressure1[100:],)
plt.plot(time2[100:(collisions-1)], pressure2[100:])
# plt.plot(time3[100:499], pressure3[100:]) 
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Pa)')
plt.title('Pressure dependence on number of balls')
plt.legend(['10 balls','50 balls'])
                
'''
See that when P and T stabilize, they reach a greater constant value for an 
increasing number of particles.

We already know P T values stay constant for a large number of collisions. 
Further check with 3 different numbers of balls. 
Printing pressures for 10, 50 and 1000 balls:

sim_ballsx1.pressure(500),sim_ballsx2.pressure(500),sim_ballsx4.pressure(500)

gives us as result 

(7.353716438029267, 7.78116046892716, 9.130814732034226) 

which confirms our hypothesis of increasing number of particles means
an increase in pressure
'''

#%% Task 12 (a)

# Change parameter k to change velocities and therefore temperature

dif_temperatures=20
dif_radius=3

T=[]
P=[]

for i in range (1, dif_temperatures):
    print(i)
    sim12=SimulationNBalls(20,1e-27,0.1,10,k=i) 
    #Considerably small radius so it can be considered ideal gas
    temperature=sim12.temperature(500)
    pressure=sim12.pressure(500)
    T.append(temperature)
    P.append(pressure)

'''
Plotting and checking results fitting a straight line: expect a straight 
line passing throguh the origin with a slope of NKb/V
'''
plt.figure()
fit,cov=np.polyfit(T,P, 1, cov=True)
pfit=np.poly1d(fit)
plt.grid()
plt.plot(T,pfit(T), color='#CC7722', linestyle='--')
plt.plot(T,P,'x', color='k')
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')
# plt.title('Checking ideal gas behaviour')
plt.legend(['Linear fit','Results from the simulation'])
print('The slope of the graph is:',fit[0],'+/-',cov[0,0],'Pa/K')

                    #Changing radius of balls
T1=[]
P1=[]

for i in range (1, dif_temperatures):
    sim12=SimulationNBalls(20,1e-27,1,10,k=i) 
    temperature=sim12.temperature(500)
    print(i)
    pressure=sim12.pressure(500)
    T1.append(temperature)
    P1.append(pressure)
    
T2=[]
P2=[]
for i in range (1, dif_temperatures):
    print(i)
    sim12=SimulationNBalls(20,1e-27,0.1,10,k=i) 
    temperature=sim12.temperature(500)
    pressure=sim12.pressure(500)
    T2.append(temperature)
    P2.append(pressure)

T3=[]
P3=[]
for i in range (1, dif_temperatures):
    print(i)
    sim12=SimulationNBalls(20,1e-27,0.01,10,k=i) 
    temperature=sim12.temperature(500)
    pressure=sim12.pressure(500)
    T3.append(temperature)
    P3.append(pressure)

#Checking how the shape changes for different radius of ballls
plt.figure()
plt.grid()

fit1,cov1=np.polyfit(T1,P1, 1, cov=True)
pfit1=np.poly1d(fit1)
plt.plot(T1,pfit1(T1))
fit2,cov2=np.polyfit(T2,P2, 1, cov=True)
pfit2=np.poly1d(fit)
plt.plot(T2,pfit2(T2))
fit3,cov3=np.polyfit(T3,P3, 1, cov=True)
pfit3=np.poly1d(fit)
plt.plot(T3,pfit3(T3))

plt.plot(T1,P1,'x', color='blue')
plt.plot(T2,P2, 'x',color='red')
plt.plot(T3,P3, 'x',color='green')

plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')
# plt.title('Pessure vs temperature graph for different radius of balls')
plt.legend(['Radius 1','Radius 0.1', 'Radius 0.01'])

print('The slopes of the 3 different lines are:')
print('Radius 1:',fit1[0])
print('Radius 0.1:',fit2[0])
print('Radius 0.01:',fit3[0])

'''
To simulate an ideal gas, the radius of the balls has to be negligible 
compared to that of the container (point particles). 

Printing the values of the slopes we see that from a radius of 0.1 the slopes
start taking a similar value. We could affirm that at this point the radius
of the balls are small enough to observe ideal gas behaviour. 
This slope shouldten be NKb/V (V being the 'volume' area of the container)

'''

#%% Task 12 (further investigation)

                        #PT ratio depending on the number of particles
                        
collisions=500 #Large number so P, T are stable

n_of_balls=[]
PT_ratio=[]

for i in range (10):
    print(i)
    sim12=SimulationNBalls(10+10*i,1e-27,0.1,10) 
    pressure=sim12.pressure(collisions)
    temperature=np.round(sim12.temperature(collisions),3)
    n_of_balls.append(10+10*i)
    PT_ratio.append(pressure/temperature)
    
plt.figure()
plt.grid()
plt.plot(n_of_balls,PT_ratio, 'x')
fit,cov=np.polyfit(n_of_balls,PT_ratio, 1, cov=True)
pfit=np.poly1d(fit)
plt.plot(n_of_balls,pfit(n_of_balls))
plt.ylabel('P/T ratio')
plt.xlabel('Number of particles')
plt.title('Pressure vs Temperature graph for changing number of balls')
plt.legend(['Obtained results','Linear fit'])

'''
Remember from last task that for sufficiently small balls, they can be treated
as an ideal gas, so the PT ratio against number of particles should give us
a straight line passing through the origin with slope Kb/V.

'''

#%% Task 12 (b)
                        #PT ratio depending on the radius of particles
                        
collisions=500 #Large number so P, T are stable

r_of_balls=[1,0.9,0.8,0.7,0.6,0.5,0.3,0.2,0.1,0.05]
PT_ratio=[]

for i in range (len(r_of_balls)):
    print(i)
    sim7=SimulationNBalls(20,1e-27,r_of_balls[i],10) 
    pressure=sim7.pressure(collisions)
    temperature=np.round(sim7.temperature(collisions),4)
    PT_ratio.append(pressure/temperature)
    
def quad (x, A, B):
    return A*x**2+B

# fit, cov=curve_fit(Max_boltz, centres, heights, p0=guess) 
# data_fit =  Max_boltz (v, *fit)
# plt.plot(v,data_fit, color='#CC7722')
x=np.linspace(0,1.02,100)

fit, cov=curve_fit(quad, r_of_balls, PT_ratio, p0=[2e-24,0.5e-24]) 
data_fit =quad(x, *fit)
plt.figure()
plt.grid()
plt.plot(r_of_balls,PT_ratio, 'x', color='k')
plt.plot(x,data_fit, color='#CC7722', linestyle='--')
plt.ylabel('P/T ratio (Pa K$^{-1}$)')
plt.xlabel('Radius of particles (m)')
# plt.title('Pressure vs Temperature graph for changing radius of balls')
plt.legend(['Results from the simulation', 'Fit'])

'''
The temperature should not depend on the radius of particles (it is only a 
function of the average kinetic energy of the particles). Nonetheless, the 
pressure should increase as the time taken to reach the walls should decrease.
This is observed in the graph.
'''

#%% Task 13

sim13a=SimulationNBalls(50,1e-27,0.1,10)
sim13a.run(1000) 
#Want to reach thermal equilibrium with a large number of collisions
velocities13a=sim13a.velocities()

for i in range (1000):
    print(i)
    sim13a.run(50) 
#No need for a big number of collisions, already thermalised. Need more data
    velocities13a.extend(sim13a.velocities())

number_of_bins=30
plt.figure()
plt.grid()
plt.hist(velocities13a,number_of_bins, color='grey')
plt.xlabel('Velocities')
plt.ylabel('Nuber of particles')
plt.title('Distribution of velocities')

#Get coordinates of centres of bins for our v distribution to fit the Maxwell
ax = plt.gca()
p = ax.patches
low_left_corner=np.array([bins.get_xy() for bins in p])
width=np.array([bins.get_width() for bins in p])
height=np.array([bins.get_height() for bins in p])
centres=[]
heights=[]
hist_err=[] #Statistical error due to random v
for i in range(number_of_bins):
    print(i)
    centres.append(low_left_corner[:,0][i]+width[i]/2)
    heights.append(height[i])
    hist_err.append(height[i]**0.5)
#These are the points we want to fit Maxwell Boltzman to

def Max_boltz (v, A, B): #From script
    return   A*v*np.exp(B*v**2)

guess_A=1e-4
guess_B=1e-3 

''''
Calculated using the equation for a maxwell Boltzmann distribution and 
knowing the order of magnitude of the temperature, mass of the gas and 
Boltzmann constant
'''
guess=[guess_A,guess_B]
v=np.linspace(0,100,100)
fit, cov=curve_fit(Max_boltz, centres, heights, p0=guess) 
data_fit =  Max_boltz (v, *fit)
plt.plot(v,data_fit, color='#CC7722')
plt.errorbar(centres, heights, yerr=hist_err,fmt='x', color='red')
plt.xlabel('Velocities (m s$^{-1}$)')
plt.ylabel('Nuber of particles')
plt.title('Distribution of velocities')
plt.legend(['Maxwell-Boltzmann distribution', \
            'Velocities for 50 balls inside the container',\
                'Errorbars due to random initial velocities'])

variance=1.38e-23*sim13a.temperature(500)/sim13a.total_mass()
print('Error on the fit', cov)
'''
Calculated looking at the form of the MB distribution
'''
#%% Task 13 (further investigation)

                            #Now try for different temperatures
                            
sim13b=SimulationNBalls(10,1e-27,0.1,10, k=2)
sim13b.run(1000)

v_hist_t2=sim13b.velocities()

for i in range (100):
    print(i)
    sim13b.run(50) 
    v_hist_t2.extend(sim13b.velocities())

plt.figure()
plt.grid()
plt.hist(velocities13a)
plt.hist(v_hist_t2)
plt.xlabel('Velocities')
plt.ylabel('Nuber of particles')
plt.title('Distribution of velocities for changing temperatures')
plt.legend(['T1', 'T2=4T1'])

'''
Behaves just as expected for a Maxwell Boltzman distribution for changing
temperature: increasing T results in a flatter curve displaced to the right

'''

#%% Task 14

'''
- a accounts for the interaction between particles, which hasn't been included
in our coude. Therefore a=0.

- b accounts for the particles not being point particles (as asumed by the 
ideal gas equation). In reality, the particles occupy a volume and their 
centres never reach the container. 
It can be seen as a 'volume'correction factor.

Observing the equation and making a=0 we can see that if we plot KbT/P against
V, b will be the value of our intercept.

For a given number of collisions (large enough so the system is thermalised),
calculate P and T of the system for a changing radius of container 
(and therefore volume)

'''

n_of_radius=20
Kb=1.38e-23
KbT_P=[]
volume=[]

for i in range(n_of_radius):
    print(i)
    sim14=SimulationNBalls(20,1e-27,0.1,10+2*i)
    pressure=sim14.pressure(500)
    temperature=np.round(sim14.temperature(500),4)
    KbT_P.append(Kb*temperature/pressure)
    volume.append(np.pi*(10+2*i)**2)

plt.figure()
fit,cov=np.polyfit(volume,KbT_P, 1, cov=True)
pfit=np.poly1d(fit)
plt.plot(volume,KbT_P,'x', color='black')
plt.plot(volume,pfit(volume), linestyle='dashed', color='#CC7722')
plt.title('Van der Waals law')
plt.ylabel('K$_{b}$T/P ((m$^3$)')
plt.xlabel('V (m$^3$)')
plt.legend(['Simulated data','Linear fit'])
plt.grid()

print('b=',fit[1],'+/-',cov[1,1])

#Check if its corrrect cause we know the slope should be 1/number of balls
print('Slope of the fit=',fit[0],'+/-',np.sqrt(cov[0,0]))
print('Sould equal', 1/20)

#%% Task 14 (further investigation)

'''
b should increae if we increased the radius of the balls, as it would have to
account for move volume
'''
n_of_radius=20
Kb=1.38e-23
KbT_P=[]
volume=[]

for i in range(n_of_radius):
    print(i)
    sim14=SimulationNBalls(20,1e-27,1,10+2*i)
    pressure=sim14.pressure(500)
    temperature=np.round(sim14.temperature(500),4)
    KbT_P.append(Kb*temperature/pressure)
    volume.append(np.pi*(10+2*i)**2) #Area (volume) of container 

plt.figure()
fit,cov=np.polyfit(volume,KbT_P, 1, cov=True)
pfit=np.poly1d(fit)
plt.plot(volume,KbT_P,'x')
plt.plot(volume,pfit(volume))
plt.title('Van der Waals law')
plt.ylabel('K$_{b}$T/P')
plt.xlabel('V (m$^3$)')
plt.legend(['Simulated data','Linear fit'])
plt.grid()

print('b=',fit[1],'+/-',np.sqrt(cov[1,1]))


#%% Task 14 (further investigation 2)

'''
If the gas is made to be non ideal (radius of the particles=1 so they are 
not point partiles), the slope should stop being 1/N
'''
n_of_radius=20
Kb=1.38e-23
KbT_P=[]
volume=[]

for i in range(n_of_radius):
    print(i)
    sim14=SimulationNBalls(8,1e-27,2,10+2*i)
    pressure=sim14.pressure(500)
    temperature=np.round(sim14.temperature(500),4)
    KbT_P.append(Kb*temperature/pressure)
    volume.append(np.pi*(10+2*i)**2) #Area (volume) of container 

plt.figure()
fit,cov=np.polyfit(volume,KbT_P, 1, cov=True)
pfit=np.poly1d(fit)
plt.plot(volume,KbT_P,'x')
plt.plot(volume,pfit(volume))
plt.title('Van der Waals law')
plt.ylabel('K$_{b}$T/P')
plt.xlabel('V (m$^3$)')
plt.legend(['Simulated data','Linear fit'])
plt.grid()

print('b=',fit[1],'+/-',np.sqrt(cov[1,1]))
print('Slope=',fit[0],'+/-',np.sqrt(cov[0,0]))
print('Shouldnt be',1/20)

#%% Task 14 (trying to plot against 1/N) (used in report)

n_of_radius=20
Kb=1.38e-23
KbT_P=[]
N_inv=[]

for i in range(n_of_radius):
    print(i)
    sim14=SimulationNBalls(10+5*i,1e-27,0.1,10)
    pressure=sim14.pressure(500)
    temperature=np.round(sim14.temperature(500),4)
    KbT_P.append(Kb*temperature/pressure)
    N_inv.append(1/(10+5*i))

plt.figure()
fit,cov=np.polyfit(N_inv,KbT_P, 1, cov=True)
pfit=np.poly1d(fit)
plt.plot(N_inv,KbT_P,'x', color='black')
plt.plot(N_inv,pfit(N_inv), linestyle='dashed', color='#CC7722')
# plt.title('Van der Waals law')
plt.ylabel('K$_{b}$T/P (m$^3$)')
plt.xlabel('1/N')
plt.legend(['Simulated data','Linear fit'])
plt.grid()

print('b=',fit[1],'+/-',cov[1,1])

#Check if its corrrect cause we know the slope should be 1/number of balls
print('Slope of the fit=',fit[0],'+/-',np.sqrt(cov[0,0]))
print('Sould equal', 1/20)
