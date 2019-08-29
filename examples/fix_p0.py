#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:37:26 2018

@author: linchen
"""

import numpy as np
from matplotlib import pyplot as plt 
from flow import Flow

# Load data
flow = Flow("fix_p0.txt","FixMassFlux")
flow.load()

#%% Plot bifurcation curves
fig_gbifur = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig_gbifur.gca()
ax.plot(flow.Q[0:], flow.height[0:],'b-o',markersize=4,label='fixed $p_0$')
ax.grid(linestyle="--")
ax.legend()
plt.xlabel(r'$Q$')
plt.ylabel('Wave height')
plt.tight_layout()

#%% Plot profiles
n_point = flow.n_q-1
flow.extract(flow.n_h-1)
flow.compute_derivative()
fig_profile1 = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig_profile1.gca()
ax.plot(flow.q2,flow.y2[:,0:-1:40],'b-',linewidth=1.2)
ax.plot(flow.q2,flow.y2[:,0],'-m',linewidth=1.2)
ax.grid(linestyle="--")
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.tight_layout()

#%% Plot velocity
flow.compute_velocity(flow.n_h-1)
fig_u = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig_u.gca()
ax.plot(flow.p,flow.c_u[0,:].T,'b-')
ax.grid(linestyle="--")
plt.ylabel('$c-u$',fontsize=12)
plt.xlabel('$p$',fontsize=12)
ax.grid(linestyle="--")
plt.tight_layout() 
print(min(flow.c_u[0,:]))

#%% Plot variation of depth
fig_depth = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig_depth.gca()
ax.plot(flow.Q[0:], flow.d[0:],'b-o',markersize=4)
ax.grid(linestyle="--")
plt.xlabel('$Q$',fontsize=12)
plt.ylabel('$d$',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()

#%% Plot pressure
flow.compute_pressure(0,flow.n_h-1)
fig = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig.gca()
im=ax.contourf(flow.x,flow.y2,flow.pres,20,cmap="RdBu_r")
ax.contour(flow.x,flow.y2,flow.pres,20,linewidths=0.5)
fig.colorbar(im)
plt.xlabel('x')
plt.xlabel('y')