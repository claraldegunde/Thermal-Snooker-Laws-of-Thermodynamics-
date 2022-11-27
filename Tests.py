#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 14:39:18 2021

@author: claraaldegundemanteca
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 14:37:52 2021

@author: claraaldegundemanteca
"""

from BallClass import *
from SingleBallSimulation import *
from MultipleBallsSimulation import *
from NBallsSimulation import *
from scipy.optimize import curve_fit
import pylab as pl 
import matplotlib.pyplot as plt
import numpy as np

#%%Testing the ball class

                    #Colliding two balls
ball1=ball=Ball(1.,1.,[0,0],[0,0], iscontainer=False)
ball2=Ball(1.,1.,[-5,0],[1,0], iscontainer=False)

print('Time to collide, expect 3:',ball1.time_to_collision(ball2))
ball1.collide(ball2)
print('Velocity of ball1 after collision, expect [1,0]', ball1.v())

                    #Now x and y velocities
ball1=ball=Ball(1.,1.,[0,0],[0,0], iscontainer=False)
ball3=Ball(1.,1.,[-5,-5],[1,1], iscontainer=False)
ball1.collide(ball3)
print('Velocity of ball1 after collision, expect [1,1]', ball1.v())

                    #Colliding a ball with a container 
ball=Ball(1.,1.,[0,0],[3,5], iscontainer=False)
container=Ball(rad=10,iscontainer=True)
ball.move(ball.time_to_collision(container))
ball.collide(container)
print('Velocity of ball after collision, expect change signs',ball.v())
        
                    #Check that the p and KE conservation exceptions work                     
ball1=ball=Ball(1.,1.,[0,2],[0,0], iscontainer=False)
ball2=Ball(1.,1.,[-5,0],[1,0], iscontainer=False)

ball1.collide(ball2)

#Collide them changing if statement to ensure it brings up an error

#%%Testing the single ball simulation

ball1=Ball(1.,1.,[0,0],[1,0], iscontainer=False)
container=Ball(rad=10,iscontainer=True)
sim1=SimulationSingleBall(ball1,container)
print('Velocity after', ball1.number_of_collisions(),'th collision',ball1.v())
sim1.next_collision()
print('Velocity after', ball1.number_of_collisions(),'th collision',ball1.v())
sim1.next_collision()
print('Velocity after', ball1.number_of_collisions(),'th collision',ball1.v())
sim1.next_collision()
print('Velocity after', ball1.number_of_collisions(),'th collision',ball1.v())
print('Bounces as expected')

                    #Now with x and y components
ball1=Ball(1.,1.,[0,0],[2,3], iscontainer=False)
container=Ball(rad=10,iscontainer=True)
sim1=SimulationSingleBall(ball1,container)
print('Velocity after', ball1.number_of_collisions(),'th collision',ball1.v())
sim1.next_collision()
print('Velocity after', ball1.number_of_collisions(),'th collision',ball1.v())
sim1.next_collision()
print('Velocity after', ball1.number_of_collisions(),'th collision',ball1.v())
sim1.next_collision()
print('Velocity after', ball1.number_of_collisions(),'th collision',ball1.v())
print('Bounces as expected')
sim1.run(10, animate=True)

                            # Now not starting from origin
ball1=ball=Ball(1.,1.,[-3,5],[2,3], iscontainer=False)
container=Ball(rad=10,iscontainer=True)
sim1=SimulationSingleBall(ball1,container)
sim1.run(10, animate=True) #Animation looks as it should


#%%Testing the multiple ball simulation, creating a list of balls

container=Ball(rad=10.,iscontainer=True) 
ball1=ball=Ball(1.,1.,[0,0],[2,3], iscontainer=False)
ball2=Ball(1.,1.,[0,0],[0,1], iscontainer=False, color='brown')
ball3=Ball(1.,1.,[0,3],[0,1], iscontainer=False, color='purple')

sim2=SimulationMultipleBalls([container,ball1,ball2, ball3])
sim2.run(30, True)

#Try with different initial positions and velocities

#%%Testing n balls

sim1=SimulationNBalls(10,1,1,10) #For 200 balls, had to reduce radius, 0.3
#Change number of balls, works fine
sim1.run(10, True) #For ten collisions, change also to check if correct

                                #Histograms
sim1.distance_to_centre(10,10) 
sim1.distance_between_balls(10, 10)  
#No one going outside of the container, correct
                    
                        #Qucickly chechin KE conservation
print(sim1.total_KE(8))
print(sim1.total_KE(3))
print(sim1.total_KE(5))
#Same value, independent of the number of collisions occured, correct


        #Chech delta p increases each time it hits container, print statements
        
        
                                #Check pressure
print(sim1.pressure())


"""
compare histograms with different size of balls, 
see if they go more towards the container etc
"""

"""
KE related to temperature through Avogadro's number and constant of gases,
see NBallsSimulation method Temperature to see the exact relation.

Momentum is also conserved, no external forces

No attraction/interaction between particles (apart from collisions) so 
can model an ideal gas

If we double all velocities, the change in momentum will be greater, 
generating a greater pressure. This makes sense as the number of collisions 
with the container per unit time would increase. The temperature will also 
increase as the KE if the particles has risen as a consequence of the 
increment in velocity.
In the simulation we expect to see balls moving quicker. 

Compare two values for T and P for the same number of balls, number of
collisions, ball characteristics... But multiplying by 2 the random velocity 
generator.
"""

#%%Checking distance to centre and separation between balls histograms

sim2=SimulationNBalls(50,1e-27,0.1,10)

sim2.distance_to_centre(10000,20)
sim2.distance_between_balls(10000, 20)

# plt.hist(distances)
# plt.xlabel('Distance to centre')
# plt.title('Histogram of distances to the centre for 500 collisions and \
    # 20 balls')
# plt.ylabel('Number of particles')
# plt.xlim(-1,10)
# plt.xticks(np.arange(-1,10,step=1))

# # sim2.distance_between_balls(500, 20)





