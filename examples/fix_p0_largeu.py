#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:37:26 2018

@author: linchen
"""

import numpy as np
import matplotlib.pyplot as plt
from flow import Flow

SMALL_SIZE = 16
MEDIUM_SIZE = 17
BIGGER_SIZE = 18
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
plt.rc('text', usetex=False)
plt.rc('savefig',dpi=300) 
plt.rc('font',family='serif')
plt.rc('axes', edgecolor='k')
plt.rc('lines', linewidth=2)
plt.rc('grid', color='k',alpha=0.6)
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
MARKER_SIZE=5

#%% For the case of u<c
flow = Flow("largeu/u_less_c.txt","FixMassFlux")
flow.load()

#%%
fig = plt.figure(figsize=(6,3.6))
ax = fig.gca()
ax.plot(flow.Q[0:], flow.height[0:],'b-o',markersize=4)
ax.grid(linestyle="--")
plt.xlabel('$Q$')
plt.ylabel('Wave height')
plt.tight_layout()
plt.ylim([-0.02,0.62])
plt.savefig('fig3a-diagram.pdf', format='pdf')

#%%
flow.extract(flow.n_h-1)
flow.compute_derivative()
fig = plt.figure(figsize=(6,3.6))
ax = fig.gca()
ax.plot(flow.q2,flow.y2[:,100:-50:50],'b-')
ax.plot(flow.q2,flow.y2[:,-1],'-b')
ax.plot(flow.q2,flow.y2[:,0],'-b')

ax.grid(linestyle="--")
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.xlim([-np.pi,np.pi]) 
plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi], 
           ('$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'))
plt.ylim([-0.8,0.5])
plt.tight_layout()
plt.savefig('fig4a-prof.pdf', format='pdf')

#fig_angle = plt.figure(figsize=(6*1.,3.6*1.))
#ax = fig_angle.gca()
#ax.plot(flow.q,np.arctan(flow.hq[:,-1])/np.pi*180)
#print(np.arctan(np.min(flow.hq[:,-1]))/np.pi*180)
#plt.tight_layout()

#%%
flow.compute_velocity(flow.n_h-1)
fig_u = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig_u.gca()
ax.plot(flow.p,flow.c_u[0,:],'b-')
plt.ylabel('$c-u$')
plt.xlabel('$p$')
ax.grid(linestyle="--")
plt.tight_layout() 
plt.savefig('fig5a-c-u.pdf', format='pdf')

#%%
fig = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig.gca()
ax.plot(flow.Q[0:], flow.d[0:],'b-o',markersize=4)
ax.grid(linestyle="--")
plt.xlabel('$Q$')
plt.ylabel('$d$')
plt.tight_layout()
plt.savefig('fig6a-depth.pdf', format='pdf')

#%% Compute slope versus mu
k = np.array([0,2,5,10])
nk = k.shape[0]
mu = np.zeros((flow.n_h, nk))
c = np.zeros((flow.n_h, nk))
theta = np.zeros((flow.n_h, nk))
qc = np.zeros((flow.n_h,nk))
for j in range(0, nk):
    for i in range(0, flow.n_h):
        flow.extract(i)
        flow.compute_derivative()
        flow.compute_velocity(i)
        c[i,j] = k[j]+np.trapz(flow.c_u[:,0],x=flow.q)/np.pi
        qc[i,j] = np.abs(-flow.c_u[0,-1]-k[j])
        mu[i,j] = np.power(c[i,j]/qc[i,j], 3.0)
        theta[i,j] = np.abs(np.min(np.arctan(flow.hq[:,-1])))/np.pi*180
        print([c[i,j],qc[i,j]],theta[i,j])

#%%
plt.plot(mu[:,0], theta[:,0], 'k:v', markersize=4,label='$k$='+'{:.2f}'.format(k[0]))
plt.plot(mu[:,1], theta[:,1], 'b-o', markersize=4,label='$k$='+'{:.2f}'.format(k[1]))
plt.plot(mu[:,2], theta[:,2], 'm--*', markersize=5,label='$k$='+'{:.2f}'.format(k[2]))
plt.plot(mu[:,3], theta[:,3], 'y-.s', markersize=4,label='$k$='+'{:.2f}'.format(k[3])) 
plt.legend(loc=4, framealpha=0.4)
plt.grid(linestyle="--")
plt.xlabel('$\mu=c^3/q_C^3$')
plt.ylabel('$\\theta_m~(^\circ)$')   
plt.xlim([0.8,9.2])
plt.xticks(np.arange(1,10,1))
plt.tight_layout()
#plt.savefig('fig7a-mu-theta.pdf', format='pdf')
plt.savefig('fig7a-mu-theta.png', format='png')
#%% Particle path
flow.compute_velocity(flow.n_h-1)
k = 10
c = k + np.trapz(flow.c_u[:,0],x=flow.q)/np.pi

flow.prepare_path_computation(c)
ix = -1
n_per = 4
t_max = 2*np.pi/c*n_per
t = np.linspace(0,t_max,200*n_per)

idx = np.where(flow.p==-0.5)[0]
traj1 = flow.compute_path([flow.q[ix],flow.y[ix,flow.n_p-1]], t)
traj2 = flow.compute_path([flow.q[ix],flow.y[ix,450]], t)
traj3 = flow.compute_path([flow.q[ix],flow.y[ix,300]], t)
traj4 = flow.compute_path([flow.q[ix],flow.y[ix,150]], t)
traj5 = flow.compute_path([flow.q[ix],flow.y[ix,1]], t)
print(c+k)

#%
fig = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig.gca()
plt.plot(traj1[0,0],traj1[0,1],'ro', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj1[-1,0],traj1[-1,1],'rs', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj1[:,0],traj1[:,1],'r-', linewidth=2)
plt.plot(traj2[0,0],traj2[0,1],'go', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj2[-1,0],traj2[-1,1],'gs', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj2[:,0],traj2[:,1],'g-', linewidth=1.2)
plt.plot(traj3[0,0],traj3[0,1],'mo', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj3[-1,0],traj3[-1,1],'ms', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj3[:,0],traj3[:,1],'m-.', linewidth=1.2)
plt.plot(traj4[0,0],traj4[0,1],'ko', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj4[-1,0],traj4[-1,1],'ks', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj4[:,0],traj4[:,1],'k--', linewidth=2)
plt.plot(traj5[0,0],traj5[0,1],'yo', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj5[-1,0],traj5[-1,1],'ys', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj5[:,0],traj5[:,1],'y--', linewidth=1.2)

plt.fill_between(flow.xm+(c)*t_max, flow.ym[:,-1], flow.ym[:,0], 
                 where=flow.ym[:,-1] >= flow.ym[:,0], 
                 facecolor='deepskyblue', interpolate=True)

plt.text(2,0.27,r'$c=$'+'%1.3f' % c)
ax.arrow(2, 0.22, 6, 0, head_width=0.03, head_length=0.5, fc='k', ec='k')

plt.xlim([-0.1*np.pi,10.1*np.pi])
plt.xticks(np.arange(-0.*np.pi,10.1*np.pi,np.pi), 
           ('$0$', '$\pi$', '$2\pi$', '$3\pi$', '$4\pi$', '$5\pi$',
            '$6\pi$', '$7\pi$', '$8\pi$','$9\pi$', '$10\pi$'))
plt.yticks([-0.8,-0.4,0,0.4])
plt.xlabel('$X$')
plt.ylabel('$Y$')
plt.tight_layout()
#plt.savefig('fig8a-path.pdf', format='pdf')
plt.savefig('fig8a-path.png', format='png')

#%% For the case of u>c
flow = Flow("largeu/u_larger_c.txt","FixMassFlux")
flow.load()

#%
#fig_name_ext='gam='+'{:.2f}'.format(flow.gamma)
fig = plt.figure(figsize=(6,3.6))
ax = fig.gca()
ax.plot(flow.Q[0:], flow.height[0:],'b-o',markersize=4)
ax.grid(linestyle="--")
plt.xlabel('$Q$')
plt.ylabel('Wave height')
plt.tight_layout()
plt.ylim([-0.02,0.62])
plt.savefig('fig3b-diagram.pdf', format='pdf')
#plt.savefig('fig3b-diagram.png', format='png')

#%%
flow.extract(flow.n_h-1)
flow.compute_derivative()
fig = plt.figure(figsize=(6,3.6))
ax = fig.gca()
ax.plot(flow.q2,flow.y2[:,0:-1:50],'b-')
ax.plot(flow.q2,flow.y2[:,-1],'-b')

ax.grid(linestyle="--")
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.xlim([-np.pi,np.pi]) 
plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi], 
           ('$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'))
plt.ylim([-0.8,0.5])
plt.tight_layout()
plt.savefig('fig4b-prof.pdf', format='pdf')

#fig_angle = plt.figure(figsize=(6*1.,3.6*1.))
#ax = fig_angle.gca()
#ax.plot(flow.q,np.arctan(flow.hq[:,0])/np.pi*180)
#print(np.arctan(np.min(flow.hq[:,0]))/np.pi*180)
#plt.tight_layout()

#%%
flow.compute_velocity(flow.n_h-1)
fig_u = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig_u.gca()
ax.plot(flow.p,flow.c_u[0,:],'b-')
plt.ylabel('$c-u$')
plt.xlabel('$p$')
ax.grid(linestyle="--")
plt.tight_layout() 
plt.savefig('fig5b-c-u.pdf', format='pdf')
#%%
fig = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig.gca()
ax.plot(flow.Q[0:], flow.d[0:],'b-o',markersize=4)
ax.grid(linestyle="--")
plt.xlabel('$Q$')
plt.ylabel('$d$')
plt.tight_layout()
plt.savefig('fig6b-depth.pdf', format='pdf')

#%% Compute slope versus mu
k = np.array([3,5,10,20])
nk = k.shape[0]
mu = np.zeros((flow.n_h, nk))
c = np.zeros((flow.n_h, nk))
theta = np.zeros((flow.n_h, nk))
qc = np.zeros((flow.n_h,nk))
for j in range(0, nk):
    for i in range(0, flow.n_h):
        flow.extract(i)
        flow.compute_derivative()
        flow.compute_velocity(i)
        c[i,j] = k[j]+np.trapz(flow.c_u[:,-1],x=flow.q)/np.pi
        qc[i,j] = np.abs(-flow.c_u[0,0]-k[j])
        mu[i,j] = np.power(c[i,j]/qc[i,j], 3.0)
        theta[i,j] = np.abs(np.min(np.arctan(flow.hq[:,0])))/np.pi*180
        print([c[i,j],qc[i,j]],theta[i,j])

#%%
plt.plot(mu[:,0], theta[:,0], 'k:v', markersize=4,label='$k$='+'{:.2f}'.format(k[0]))
plt.plot(mu[:,1], theta[:,1], 'b-o', markersize=4,label='$k$='+'{:.2f}'.format(k[1]))
plt.plot(mu[:,2], theta[:,2], 'm--*', markersize=5,label='$k$='+'{:.2f}'.format(k[2]))
plt.plot(mu[:,3], theta[:,3], 'y-.s', markersize=4,label='$k$='+'{:.2f}'.format(k[3]))
plt.legend(loc=5, framealpha=0.4)
plt.grid(linestyle="--")
plt.xlabel('$\mu=c^3/q_C^3$')
plt.ylabel('$\\theta_m~(^\circ)$')    
plt.xlim([-0.025,1.025])
plt.tight_layout()
plt.savefig('fig7b-mu-theta.pdf', format='pdf')
#plt.savefig('fig7b-mu-theta.png', format='png')
#%% Particle path
flow.compute_velocity(flow.n_h-1)
k = 10
c = k + np.trapz(flow.c_u[:,-1],x=flow.q)/np.pi

flow.prepare_path_computation(c)
ix = -1
n_per = 3
t_max = 2*np.pi/c*n_per
t = np.linspace(0,t_max,200*n_per)

traj1 = flow.compute_path([flow.q[ix],flow.y[ix,0]], t)
traj2 = flow.compute_path([flow.q[ix],flow.y[ix,flow.n_p-1-450]], t)
traj3 = flow.compute_path([flow.q[ix],flow.y[ix,flow.n_p-1-300]], t)
traj4 = flow.compute_path([flow.q[ix],flow.y[ix,flow.n_p-1-150]], t)
traj5 = flow.compute_path([flow.q[ix],flow.y[ix,flow.n_p-1]], t)
print(c+k)

#%%
fig = plt.figure(figsize=(6*1.,3.6*1.))
ax = fig.gca()
plt.plot(traj1[0,0],traj1[0,1],'ro', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj1[-1,0],traj1[-1,1],'rs', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj1[:,0],traj1[:,1],'r-', linewidth=2)
plt.plot(traj2[0,0],traj2[0,1],'go', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj2[-1,0],traj2[-1,1],'gs', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj2[:,0],traj2[:,1],'g-', linewidth=1.2)
plt.plot(traj3[0,0],traj3[0,1],'mo', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj3[-1,0],traj3[-1,1],'ms', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj3[:,0],traj3[:,1],'m-.', linewidth=1.2)
plt.plot(traj4[0,0],traj4[0,1],'ko', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj4[-1,0],traj4[-1,1],'ks', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj4[:,0],traj4[:,1],'k--', linewidth=2)
plt.plot(traj5[0,0],traj5[0,1],'yo', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj5[-1,0],traj5[-1,1],'ys', linewidth=1.2, markersize=MARKER_SIZE)
plt.plot(traj5[:,0],traj5[:,1],'y--', linewidth=1.2)

plt.fill_between(flow.xm+c*t_max, flow.ym[:,-1], flow.ym[:,0], 
                 where=flow.ym[:,-1] >= flow.ym[:,0], 
                 facecolor='deepskyblue', interpolate=True)

plt.text(2, 0.27,r'$c=$'+'%1.3f' % c)
ax.arrow(2, 0.22, 5.5, 0, head_width=0.03, head_length=0.5, fc='k', ec='k')

plt.xlim([-0.1*np.pi,10.1*np.pi])
plt.xticks(np.arange(-0.*np.pi,10.1*np.pi,np.pi), 
           ('$0$', '$\pi$', '$2\pi$', '$3\pi$', '$4\pi$', '$5\pi$',
            '$6\pi$', '$7\pi$', '$8\pi$','$9\pi$', '$10\pi$'))
plt.xlabel('$X$')
plt.ylabel('$Y$')
plt.yticks([-0.8,-0.4,0,0.4])
plt.tight_layout()
plt.savefig('fig8b-path.pdf', format='pdf')