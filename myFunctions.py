#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 22:38:13 2018

@author: christophernovitsky
"""

import math as m
import numpy as np 
import math as m
import numpy as np 
import  scipy.optimize


def yAxis(utm):
    count = 0
    ang = [];
    yAxis = np.zeros(shape=(96,32))
    for i in range(1,32,2):
    
        P0 = [utm.loc[0,i], utm.loc[0,i+1]] 
        P1 = [utm.loc[0,i], utm.loc[0,i+1]+20] 
        P = [utm.loc[1,i], utm.loc[1,i+1]] 
   
        v1 = np.asarray(P0)-np.asarray(P1)
        v2 = np.asarray(P0)-np.asarray(P)
        a = np.dot(v1,v2)
        b = np.linalg.norm(v1)*np.linalg.norm(v2)
  
        ang = np.arccos(a/b)*180/m.pi
        
        x = np.array(np.linspace(0,360,96)+ ang)
        for k in range(0,len(x)): 
            if (x[k] > 360): x[k] -= 360
     
        yAxis[:,count] = x
        count +=1
    return yAxis



    
    