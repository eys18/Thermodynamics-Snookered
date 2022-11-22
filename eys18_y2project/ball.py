#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:00:20 2020

@author: ewansaw

Ball Class and Container Class
"""
import numpy as np
import pylab as pl

class Ball:
    ballcount = 0
    def __init__(self, m = 1.0, R = 0.5, r = np.array([0.0,0.0]), v = np.array([0.0,0.0])):
         self._mass = m
         self._radius = R
         self._position = r
         self._velocity = v
         self._patch = pl.Circle(self._position, self._radius, fc='r')
         Ball.ballcount += 1
         
    def __repr__(self):
        return "%s(mass=%s, radius=%s, position=array%s velocity=array%s)" % ("Ball", 
                  self._mass, self._radius, self._position, self._velocity)

    def __str__(self):
        return "(%r, %r, %r, %r)" % (self._mass, self._radius, self.position, self.velocity)
    
    def pos(self): #to return position of ball
        print(self._position)
        
    def vel(self): #to return velocity of ball
        print(self._velocity)
        
    def get_patch(self):
        return self._patch
        
    def move(self, dt):
        self._position = self._position + self._velocity * dt
        self._patch.center = self._position
        return self
        
    def time_to_collision(self, other):
        r = self._position - other._position
        v = self._velocity - other._velocity
        if isinstance(other, Container) or isinstance(self, Container): #to check for collision with container 
            R = self._radius - other._radius
        else:
            R = self._radius + other._radius
        a = v.dot(v)
        b = 2 * r.dot(v)
        c = r.dot(r) - R ** 2
        dis = b ** 2 - 4 * (a * c) #discriminant
        if dis < 0:
            return np.inf #so it doesnt get picked
        elif a == 0:
            return np.inf
        else:      
            dt1 = (-b + np.sqrt(dis)) / (2 * a)
            dt2 = (-b - np.sqrt(dis)) / (2 * a)
            
            if dis == 0:
                return  -b / (2 * a)
            elif dt2 > 0 and dt1 > 0:
               if dt2 > dt1:
                   return dt1 - 0.0000000001 #to stop balls from stop sticking
               else:
                   return dt2 
            elif dt1 <= 0 and dt2 <= 0:
                return np.inf
            elif dt1 > 0 and dt2 <= 0:
                return dt1 - 0.0000000001
            elif dt2 > 0 and dt1 <= 0:
                return dt2 - 0.0000000001
    
    
    def collide(self, other):
        ke1 = 0.5 * self._mass * self._velocity.dot(self._velocity) + 0.5 * other._mass * other._velocity.dot(other._velocity)
        r = self._position - other._position
        v = self._velocity
        v1 = other._velocity
        rmag = np.sqrt(np.dot(r, r))
        vpar = (v.dot(r) / rmag) * (r / rmag) #velocity of self parallel to vector between centres of two balls
        vpar1 = (v1.dot(r) / rmag) * (r / rmag) #velocity of other parallel to vector between centres of two balls
        vper = v - vpar #perpendicular velocity of self
        vper1 = v1 - vpar1 #perpendicular velocity of other ball
        
        self._velocity = ((self._mass - other._mass) / (self._mass + other._mass) * (vpar)
                            + (2 * other._mass) / (self._mass + other._mass) * (vpar1) + vper) #new velocity
        other._velocity = ((2 * self._mass) / (self._mass + other._mass) * (vpar)
                            - (self._mass - other._mass) / (self._mass + other._mass) * (vpar1) + vper1) #new velocity of other ball
        ke2 = 0.5 * self._mass * self._velocity.dot(self._velocity) + 0.5 * other._mass * other._velocity.dot(other._velocity)     
        
        if round(ke1, 5) != round(ke2, 5):
            raise TypeError("KE not conserved") #to check if collision is allowed
        else:
            return self
    
    def ke(self):
        KE = (1 / 2) * self._mass * self._velocity.dot(self._velocity)
        return KE


class Container(Ball):
    """
Inherited from Ball class
    """
    def __init__(self):
        m = 999999999999 #assumed large mass for container
        R = 15 
        v = np.array([0.0, 0.0])
        r = np.array([0.0, 0.0])
        Ball.__init__(self, m, R, r, v)
        self._patch = pl.Circle([0., 0.], 15, ec='b', fill=False, ls='solid')
        Ball.ballcount -= 1
        
    def __repr__(self):
        return "%s(mass=%s, radius=%s, position=array%s velocity=array%s)" % ("Container", 
                  self._mass, self._radius, self._position, self._velocity)
