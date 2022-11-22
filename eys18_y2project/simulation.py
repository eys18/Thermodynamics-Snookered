#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:05:54 2020

@author: ewansaw
"""

import ball as bl
import pylab as pl
import numpy as np
import copy
import scipy.stats as stats

class Simulation:

    def __init__(self, num_balls): #will only work with square numbers!!!
        bl.Ball.ballcount = 0
        self._N = num_balls
        self._container = bl.Container()
        self._simulation_time = 0 #absolute time of system
        self._simulation_timelist = [0]  #times of collisions
        self._distance = []  # list of distances from centre
        self._balldist = []  # list of distances between balls
        self._ball = []
        self._ke = [] #KEs of system
        self._impulses = [] #impulses of balls on container
        self._impulse_time = []
        self._velocities = []
        
        a = np.linspace(-10, 10, num = round(np.sqrt(self._N))) #initialisation of balls
        positions = []
        for j in a:
            for i in a:
                positioni = [j, i]
                positions.append(positioni)            
        for i in range(self._N):
            self._ball.append(bl.Ball(r = positions[i], v = np.array([np.random.uniform(-15, 15), np.random.uniform(-15, 15)]))) 
    
        KEs = []
        for i in range(self._N):
            KEs.append(self._ball[i].ke())
        TotalKE = round(sum(KEs), 5)
        self._ke.append(TotalKE)  # list of the total ke of the system after each collision

    def next_collision(self):
        for n in range(self._N):
            self._ball[n].move(0.00000001) #to prevent sticking of balls
        r = bl.Ball.ballcount
        c = r + 1
        collns = np.zeros(shape=[r,c]) #matrix of balls
        for i in range(r):
            for j in range(i + 1, c - 1):
                collns[i,j] = self._ball[i].time_to_collision(self._ball[j]) #ball with ball
            collns[i, c - 1] = self._ball[i].time_to_collision(self._container) #ball with container

        t = np.min(collns[np.nonzero(collns)]) #finding minimum dt
        for i in range(r):
            self._ball[i] = self._ball[i].move(t)
        self._simulation_time += t
        self._simulation_timelist.append(self._simulation_time)
    
        colln = np.where(collns == np.min(collns[np.nonzero(collns)])) #finding coordinates of smallest dt
        if colln[1][0] == (c - 1):
            r = self._ball[colln[0][0]]._position - self._container._position
            v = self._ball[colln[0][0]]._velocity
            rmag = np.sqrt(np.dot(r, r))
            vpar = (v.dot(r) / rmag)
            self._ball[colln[0][0]].collide(self._container)
            self._impulses.append(2 * self._ball[colln[0][0]]._mass * vpar) #updating list of impulses
            self._impulse_time.append(self._simulation_time)
        else:
            self._ball[colln[0][0]].collide(self._ball[colln[1][0]])     
        
        KEs = [] #updating KE list after each collision
        for i in range(self._N):
            KEs.append(self._ball[i].ke())
        TotalKE = round(sum(KEs), 5)
        self._ke.append(TotalKE)  
    
        for i in range(self._N): #updating velocity list
            self._velocities.append(np.sqrt(self._ball[i]._velocity
                                            .dot(self._ball[i]._velocity)))

    def run(self, num_frames, animate=False):
        if animate:
            f = pl.figure(figsize = (5, 5))
            ax = pl.axes(xlim=(-15, 15), ylim=(-15, 15))
            ax.add_artist(self._container.get_patch())
            for i in range(bl.Ball.ballcount):
                ax.add_patch(self._ball[i].get_patch())
        for frame in range(num_frames):
            self.next_collision()
            if animate:
                pl.pause(0.03)
        if animate:
            pl.show()
        
    def kin_energy(self, num_frames = 1000):
        self.run(num_frames)
        pl.plot(self._simulation_timelist, self._ke, label = 'Kinetic Energy')
        pl.grid()
        pl.title("Kinetic Energy over Time")
        pl.xlabel("Time (t)")
        pl.ylabel("Kinetic Energy")
        pl.legend()
        
    def temperature(self, num_frames = 1000):
        self.run(num_frames)
        temp = np.asarray(self._ke) / (1.38e-23 * self._N)
        print(temp)
        pl.plot(self._simulation_timelist, temp, label = 'Temperature')
        pl.grid()
        pl.title("Temperature over Time")
        pl.xlabel("Time (t)")
        pl.ylabel("Temperature (T)")
        pl.legend()
  
    def pressure(self, num_frames = 100, step = 50): #calculating pressure over step number of impulse collisions
        self.run(num_frames)
        time_interval = [self._impulse_time[step - 1]] #time interval over which set of impulses occur
        time_axis = [self._impulse_time[step - 1]]      
        for i in range(1, int(np.floor(len(self._impulse_time)) / step)): #appending both time and time interval lists
            time_interval.append(self._impulse_time[step * (i + 1) - 1] - self._impulse_time[step * i]) 
        for i in range(1, int(np.floor(len(self._impulse_time)) / step)):
            time_axis.append(self._impulse_time[step * (i + 1) - 1]) 
        cutoff_point = int(np.floor(int(len(self._impulses)) / step) * step) #to correct for proper indexing of lists
        imp = copy.copy(self._impulses)
        del imp[cutoff_point:] #deleting extra data of impulses that are recorded past the required index
        impulses = np.add.reduceat(imp, np.arange(0, int(np.floor(int(len(self._impulses)) / step)) * step, step)) #summing impulses
        pressures = (impulses / time_interval) / (2 * np.pi * self._container._radius) #calculating pressure
        pl.axes()
        pressures = np.insert(pressures, 0, 0)
        time_axis = np.insert(time_axis, 0, 0)
        print(pressures) 
        pl.plot(time_axis, pressures)
        pl.gca().set_ylim(ymin=0)
        pl.grid()
        pl.title("Pressure over Time")
        pl.xlabel("Time (t)")
        pl.ylabel("Pressure (N/m)")
        pl.plot(time_axis, pressures)
        
    def max_boltz(self):
        pl.hist(self._velocities, 50, (0,30), edgecolor = 'black', density = True)
        params = stats.maxwell.fit(self._velocities) #fit parameters for maxwell boltzmann distribution
        x = np.linspace(0,30,31)
        pl.plot(x, stats.maxwell.pdf(x, *params), lw = 2)
        pl.xlabel("Velocity (ms^-2)")
        pl.ylabel("Probability Density")
        pl.show()
        
    def centraldist(self, num_frames): #distances of balls from the centre of the container
        self._distance = []
        for frame in range(num_frames):
            self.next_collision()
            for i in range(self._N):
                self._distance.append(np.sqrt(self._ball[i]._position.dot(self._ball[i]._position)))
        pl.hist(self._distance, 15, (0,15), edgecolor = 'black', density = True)
        pl.show()
        
    def balldist(self): #distances of balls from each other
        balls = bl.Ball.ballcount
        distances = np.zeros(shape=[balls, balls]) #matrix of balls
        for i in range(balls):
            for j in range(i + 1, balls):
                r = self._ball[i]._position - self._ball[j]._position
                distances[i,j] = np.sqrt(r.dot(r))
        x, y = np.nonzero(distances)
        distances = distances[x, y]
        pl.hist(distances, 15, (0,15), density = True) 
        pl.show()