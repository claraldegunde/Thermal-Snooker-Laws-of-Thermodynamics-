#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 15:52:26 2020

@author: Clara Aldegunde Manteca (ca520)
"""

#Importing data-----------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit #Importing needed modules


D=np.loadtxt('Distance_Mpc.csv', skiprows=1, delimiter=',', unpack=False)  #Distance data
spec_data=np.loadtxt('Halpha_spectral_data.csv',skiprows=2, delimiter=',', unpack=False)  #Spectral data (intensity against wavelength)


#Selecting data with good instrument response----------------------------------


                            #Distance data
                            
D_clean=[] 
obs=[]

for i in range(0,30):
    if (D[i,2])==1: #If valid instrument response is 1 (third column)
        D_clean.append(D[i,1]) #then add the value of the distance to the D_clean list
        obs.append(D[i,0]) #and record the observation number (first column)


D_clean=np.array(D_clean) #Array of correct distances
obs=np.array(obs) #Correct observation numbers
obs_sort=np.sort(obs) #Correct observations sorted


dist_indices=[]

for i in range (0,25):
    dist_indices.append(np.where(obs==obs_sort[i])[0])

dist_indices=np.array(dist_indices) #Positions of the distance data for each observation nummber

dist=D_clean[dist_indices] #Distances in ascendant order of observation number
dist=np.resize(dist, (25,))


#Distances are in megaparsec
#dist is the variable to represent against velocity


                              # Spectral data
                              

with open ('Halpha_spectral_data.csv', 'r') as file: #To select the needed rows
    line1=file.readline()
    line2=file.readline() #Observations line
line2=line2.replace('Observation: ','') #Removing the word observation
line2=line2.split(',') #, as separation
line2_array=np.array(line2[1::]) #Create an array (strings)

line2_int=[] #Strings to integers
for i in range(0, len(line2_array)): 
    line2_int.append(int(line2_array[i])) 
observation_number=np.array(line2_int) #Observation numbers from spectral data

spec_index=[] #Emtpy list to fill with the correct indices
for i in range (0,25):
    spec_index.append(np.where(observation_number==obs_sort[i])[0]) #Return the column indices of spec for which the observations are correct (sorted)
    
spec_index=np.array(spec_index) #Correct column indices to analyse
spec_index=np.resize(spec_index,(25,))

spec_sliced=spec_data[:,1:] #Skip first column:wavelength

spec=[]

for i in range (0,25):
   spec.append(spec_sliced[:,spec_index[i]])
   
spec=np.array(spec).T #Intensity values(each column is one observation)
                
#See that we sorted the observation numbers in both distance and spectral data, making sure we link the same observation numbers in both data sets


#Intensity against wavelenght plots-------------------------------------------------


in_var=spec_data[:,0]  #Independent variable (x). Wavelength in m
dep_var=spec #Dependent variable, intensity values


for i in range (0,25):
   plt.figure()    
   plt.plot(in_var*1e9,dep_var[:,i], color='brown')
   plt.xlabel('Wavelength (nm)') #In nm as we multiplied the independent variable x1e9
   plt.ylabel('Intensity (arbitrary units)')
   plt.title('Intensity against wavelength for observation %d' % obs[i])
   plt.grid()
    
    
#Guesses for the fit ----------------------------------------------------------------
#To represent the gausian+line fit we will need several rough guesses to introduce in the fit function


                          # Guessing slope and intercept of the linear background


m=[]
b=[]

for i in range (0,25):
    m.append((dep_var[-1,i]-dep_var[0,i])/(in_var[-1]-in_var[0])) #Guess for slope using the first and last points of spectral data
    b.append(dep_var[-1,i]-m[i]*in_var[-1]) #Guess for b knowing m and following the equation of a straight line (y=mx+b)
    
m=np.array(m) 
b=np.array(b)  #Arrays for slope and y intercept respectively


                          # Guessing amplitude of the Gaussian function

    
y_values=[] #List of the y values for each straight line

for i in range (0,25):
    def straight_line (x):
        y=m[i]*x+b[i]
        return y
    y_values.append(straight_line(in_var))

y_values=np.array(y_values) #Each column of this array records the y values of the guess line for the different observations
y_values=y_values.T

diff_max=[]
diff=[]

for i in range (0,25):
    diff_max.append(np.max(dep_var[:,i]-y_values[:,i])) #Each column of y_values is one observation 
    diff.append(dep_var[:,i]-y_values[:,i])
    
diff_max=np.array(diff_max).T  #Guessed amplitude for each observation
diff=np.array(diff).T 


                            # Guessing mu (location of the peak)
                            
                            
index_max=[]

for i in range (0,25):
    index_max.append(np.where(diff[:,i]==diff_max[i]))

index_max=np.array(index_max) #This array returns where the maximum difference was reached (peak)
index_max.resize((25,))

mu=[]

for i in range (0,25):
    mu.append(in_var[index_max[i]]) #x values with the calculated indices

mu=np.array(mu) #Guessed location for the peak on each observation


                                 #Standard deviation
                                 
                                 
std=[]
for i in range (0,25):
    std.append(1/(diff_max[i]*np.sqrt(2*np.pi))) #Taking into account the relationship between amplitude and standard deviation (Gaussian function)
    
std=np.array(std)
                          
#Gaussian and line fitting-------------------------------------------------             


def gauss_line (in_var, a, sigma, mu, m, b): #Where in_var (x), sigma (standard deviation), a (amplitude) and mu (mean) are the parameters for the gaussian fit 
                                             #and m, b are the slope and intercept of the line   
    gaussian = a*np.exp(-(in_var-mu)**2/2*sigma**2)
    line = in_var*m+b
    return gaussian + line #Gaussian with linear background

lambda_0=np.zeros((25)) 
lambda_err=np.zeros((25))


for i in range (0,25):
# Our initial guess for the parameters
    guess_amplitude=diff_max[i]
    guess_sigma=std[i]*1e10
    guess_mu=mu[i]
    guess_m=m[i]
    guess_b=b[i]
    guess=[guess_amplitude, guess_sigma, guess_mu, guess_m, guess_b]# Array of initial parameter values
    fit, cov=curve_fit(gauss_line, in_var, dep_var[:,i], p0=guess) #To fit the Gaussian and straight line
    data_fit = gauss_line(in_var, *fit)
    plt.figure() #Plotting
    plt.xlabel('Wavelegth (nm)')
    plt.ylabel('Intensity (arbitrary units)')
    plt.plot(in_var*1e9,dep_var[:,i], color='brown')
    plt.plot(in_var*1e9,data_fit, color='orange')
    plt.title('Intensity against wavelength for observation %d' % obs[i])
    plt.grid()
    plt.show()
    lambda_0[i]=fit[2] #Observed wavelength, third fit parameter
    lambda_err[i]=np.sqrt(cov[2,2]) #Its error
                                        

#Calculating velocity-------------------------------------------------------


c=3e8 #Speed of light on ms-1
lambda_e=656.28e-9 #H-alpha wavelength in m
lambda_0=lambda_0 #Observed wavelength



v=(c*((lambda_0**2-lambda_e**2)/(lambda_e**2+lambda_0**2))) #Velocity in ms-1 for each observation
v=v*1e-3 #Velocity in km/s

#Error propagation: v is a quotient times c. The error on v will be its value times the relative error of the quotient (c without error)
#We will therefore use the formula of the error on a quotient
#The error on numerator and denominator were calculated using the error of summations and powers

num=np.zeros((25))
denom=np.zeros((25))
quot=np.zeros((25))
num_err=np.zeros((25))
denom_err=np.zeros((25))
quot_err=np.zeros((25))
v_err=np.zeros((25))



for i in range (0,25):
    num[i]=(lambda_0[i]**2-lambda_e**2) #Numerator
    denom[i]=(lambda_0[i]**2+lambda_e**2) #Denominator
    quot[i]=num[i]/denom[i] #Quotient
    num_err[i]=2*lambda_0[i]*lambda_err[i] #Error of a power (as lambda_e is considered without error)
    denom_err[i]=num_err[i] #Using error propagation for a power (lambda0^2) and for a sum (lambda0 and lambdae)
    quot_err[i]=quot[i]*np.sqrt((num_err[i]/num[i])**2+(denom_err[i]/denom[i])**2) #Error of a quotient
    v_err[i]=c*(quot_err[i]) #v times relative error of the quotient (v/quot=c)

v_err=v_err*1e-3

#Fitting a straight line

fit_hubble,cov_hubble=np.polyfit(dist,v, 1, w=1/v_err,cov=True) 
pfithubble=np.poly1d(fit_hubble)
hubble=fit_hubble[0]
hubble_err=np.sqrt(cov_hubble[0,0])

#Velocity against distance-------------------------------------------------------


plt.figure()
plt.errorbar(dist,v, yerr=v_err, fmt='+', color='black') 
plt.plot(dist,pfithubble(dist), color='brown')
plt.xlabel('Distance(Mpc)')
plt.ylabel('Velocity(km/s)')
plt.title('Velocity against distance for different galaxies')
plt.grid()
plt.savefig("Clara Aldegunde Plot.png",dpi=300)

print('The obtained value for Hubbles constant is: %0.2f +/- %0.2f km/s/Mpc.' % (hubble, hubble_err))


#Discussion---------------------------------------------------------------------

#Hubble's constst real value, obtained from computing lecture 5: 68km/s/Mpc


print('It differs from the theoretical value in a %0.f' % (((68-hubble)/68)*100) 
      +str('%'))

print('The porcentual error on the result is %0.f' % ((hubble_err/hubble)*100) 
      +str('%, which means that the real value is within our error'))