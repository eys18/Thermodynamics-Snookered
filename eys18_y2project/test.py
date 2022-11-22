#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:12:10 2020

@author: ewansaw

Test script which was used to obtain plots and relations used in report
"""
import simulation as sim #importing necessary modules and packages
import numpy as np
import pylab as pl
from numpy.polynomial.polynomial import polyfit

params = { #parameters for figures
   'axes.labelsize': 15,
   'font.size': 15,
   'legend.fontsize': 14,
   'xtick.labelsize': 15,
   'ytick.labelsize': 15,
   'figure.figsize': [9, 6]
   } 

pl.rcParams.update(params)

#%%
x = sim.Simulation(num_balls = 64) #animation of simulation (will only work with square numbers of balls)
x.run(num_frames = 1000, animate = True) 

#%%
x = sim.Simulation(num_balls = 100,) #plotting the temperature of the system over time
x.temperature(num_frames = 100)

#%%
x = sim.Simulation(num_balls = 100) #plotting the kinetic energy of the system over time
x.kin_energy(num_frames = 100)

#%%
#x = sim.Simulation(num_balls = 100) #plotting the pressure of the system over time
x.pressure(step = 100) # number of impulse collisions to calculate over 

#%%
#100 balls, data obtained from above plots for varying temperatures
pressures = np.array([5.94166648, 22.47570136, 97.45781008, 293.49823669, 406.38919627, 575.60437631])
temperatures = np.array([2.36546145e+24, 9.89079943e+24, 3.85410721e+25, 1.29827583e+26, 1.63819444e+26, 2.32670371e+26])

m, b = np.polyfit(pressures, temperatures, 1)
b, m = polyfit(pressures, temperatures, 1)
pl.plot(pressures, b + m * pressures)
pl.plot(pressures, temperatures, 'o')

pl.grid()
pl.title("Pressure over Temperature")
pl.xlabel("Temperature (K)")
pl.ylabel("Pressure (N/m)")
pl.legend()



#%%
x = sim.Simulation(num_balls = 100)
x.run(num_frames = 7000, animate = False)
x.max_boltz()
    