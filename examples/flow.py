#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23, 2018
@author: linchen
"""

#import matplotlib as mpl
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#import matplotlib.pyplot as plt
from scipy.integrate import odeint

class Flow(object):
    def __init__(self,f_name, opt):
        self.file_name=f_name
        self.type=opt.lower()
        
    def gen_grid(self,p0,n_q,n_p):
        self.p = np.linspace(p0,0,num=n_p,endpoint=True)
        self.q = np.linspace(0,np.pi,num=n_q,endpoint=True)
        self.q2 = np.linspace(-np.pi,np.pi,num=n_q*2-1,endpoint=True)
        
    def load(self):
        i_line = 0
        for line in open(self.file_name):
            i_line = i_line+1
            line_split = line.split()
            
            if i_line==1:
                pass
            elif i_line == 2:
                if len(line_split)==2:
                    self.grav = float(line_split[0])
                    self.gamma = np.zeros((1,1));
                    self.gamma[0] = float(line_split[1])
                elif len(line_split)==4:
                    self.grav = float(line_split[0])
                    self.gamma = np.zeros((2,1));
                    self.gamma[0] = float(line_split[1])
                    self.gamma[1] = float(line_split[2])
                    self.p1 = float(line_split[3])
                else:
                    pass
                
            elif i_line==3:
                self.q = np.asarray(line_split,dtype=np.float64)
            elif i_line==4:
                self.p = np.asarray(line_split,dtype=np.float64)
            else:
                break
            
        if self.gamma.shape[0]==2:
            self.ip1 = np.argmin(np.abs(self.p-self.p1))
        
        self.n_q = self.q.shape[0]
        self.n_p = self.p.shape[0]
        self.q2 = np.zeros((self.n_q*2-1,))
        self.q2[0:self.n_q]=-self.q[-1::-1]
        self.q2[self.n_q:]=self.q[1:]
        
        self.n_v = self.n_q*self.n_p+3
        
        self.vs = np.reshape(np.loadtxt(self.file_name,skiprows=4),(-1,self.n_v))
        self.n_h = self.vs.shape[0]
        self.d = np.zeros((self.n_h,1))
        self.p0 = np.zeros((self.n_h,1))
        self.Q = np.zeros((self.n_h,1))
        self.height = np.zeros((self.n_h,1))
        
        for i in range(0,self.n_h):
            self.d[i] = self.vs[i,-3]
            self.p0[i] = self.vs[i,-2]
            self.Q[i] = self.vs[i,-1]
            
        if self.type=="fixmassflux":
            for i in range(0,self.n_h):
                if self.p0[i]<0:
                    self.height[i]=self.vs[i,self.n_p-1]-self.vs[i,self.n_p*self.n_q-1]
                elif self.p0[i]>0:
                    self.height[i]=self.vs[i,0]-self.vs[i,self.n_p*(self.n_q-1)]
            
        elif self.type=="fixdepth":            
            for i in range(0,self.n_h):
                self.height[i]=(self.vs[i,self.n_p-1]-self.vs[i,self.n_p*self.n_q-1])*self.d[i]
                
    def extract(self, idx):
        hv = self.vs[idx,0:-3]
        self.h = np.reshape(hv,(self.n_q,self.n_p))
        
        if self.type=="fixmassflux":
            self.y = self.h-self.d[idx]
        elif self.type=="fixdepth":
            self.y = self.d[idx]*(self.h+self.p)
            
        self.y2 = np.zeros((self.n_q*2-1,self.n_p))
        self.y2[0:self.n_q,:]=self.y[self.n_q::-1,:]
        self.y2[self.n_q:,:]=self.y[1:,:]

    def compute_derivative(self):
        dq = self.q[1:]-self.q[0:-1]
        dp = self.p[1:]-self.p[0:-1]
        dq2i = 1/(dq[0:-1]+dq[1:])
        dp2i = 1/(dp[0:-1]+dp[1:])
        
        self.hq = np.zeros((self.n_q,self.n_p))
        self.hp = np.zeros((self.n_q,self.n_p))
        
        self.hq[1:-1,:] = (self.h[2:,:]-self.h[0:-2,:])*dq2i.reshape(self.n_q-2,1)
        self.hp[:,1:-1] = (self.h[:,2:]-self.h[:,0:-2])*dp2i.reshape(1,self.n_p-2)
        # Boundary condition
        self.hq[0,:] = (self.h[1,:]-self.h[0,:])/dq[0]
        self.hq[-1,:] = (self.h[-1,:]-self.h[-2,:])/dq[-1]
        # Boundary condition
        self.hp[:,0] = (self.h[:,1]-self.h[:,0])/dp[0]
        self.hp[:,-1] = (self.h[:,-1]-self.h[:,-2])/dp[-1]
        
    def compute_velocity(self, idx):
        self.extract(idx)
        self.compute_derivative()
        
        if self.type=="fixmassflux":
            self.c_u = 1/self.hp
            self.v = -self.c_u*self.hq
        elif self.type=="fixdepth":
            self.c_u = -self.p0[idx]/self.d[idx]/(self.hp+1)
            self.v = -self.c_u*self.hq*self.d[idx]
        
    def compute_pressure(self, Patm, idx):
        self.extract(idx)
        self.compute_derivative()
        self.x = np.repeat(self.q2.reshape(self.n_q*2-1,1),self.n_p,axis=1)
        
        if self.type=="fixmassflux":
            Gfun = self.fun_gamma()
            P = (Patm+self.Q[idx]/2-Gfun-self.grav*self.h
                 -0.5*(1+np.power(self.hq,2))
                 /np.power(self.hp,2.0))
            
        elif self.type=="fixdepth":
            Gfun = self.p0[idx]*self.fun_gamma()
            P = (Patm+self.Q[idx]/2+Gfun-self.grav*self.d[idx]*(1+self.h)
                 -self.grav*self.d[idx]*self.p
                 -0.5*self.p0[idx]**2*(1/(self.d[idx]**2)+np.power(self.hq,2))
                 /np.power(self.hp+1,2.0))
        
        self.pres = np.zeros(self.y2.shape)
        self.pres[0:self.n_q,:] = P[-1::-1,:]
        self.pres[self.n_q:,:] = P[1:,:]
    
    def fun_gamma(self):
        if self.gamma.shape[0]==1:
            G = self.gamma[0]*self.p
        elif self.gamma.shape[0]==2:
            G = np.zeros(self.p.shape)
            idx = np.where(self.p==self.p1)[0]
            G[idx[0]:] = self.gamma[0]*self.p[idx[0]:]
            G[0:idx[0]] = self.gamma[1]*(self.p[0:idx[0]]-self.p1)+G[idx[0]]
        else:
            pass
        
        return G
    
    
    def prepare_path_computation(self, c_):
        self.c = c_
        
        if self.p[0]<0:
            c_u = np.concatenate((self.c_u[-1::-1,:],self.c_u[1:,:]),axis=0)
            # Velocity field and grid coordinate for one period
            v = np.concatenate((-self.v[-1::-1,:],self.v[1:,:]),axis=0)
            u = self.c - c_u
            x = self.q2
            y = np.concatenate((self.y[-1::-1,:],self.y[1:,:]),axis=0)
        else:
            c_u = np.concatenate((self.c_u[-1::-1,-1::-1],self.c_u[1:,-1::-1]),axis=0)
            # Velocity field and grid coordinate for one period
            v = np.concatenate((-self.v[-1::-1,-1::-1],self.v[1:,-1::-1]),axis=0)
            u = self.c - c_u
            x = self.q2
            y = np.concatenate((self.y[-1::-1,-1::-1],self.y[1:,-1::-1]),axis=0)
        
        
        # Velocity field and grid coordinate for multiple period
        self.xm = np.concatenate((x-8*np.pi,x[1:-1]-6*np.pi,x[1:-1]-4*np.pi,
                                  x[1:-1]-2*np.pi,x[1:-1],x[1:-1]+2*np.pi,x+4*np.pi),axis=0)
        self.ym = np.concatenate((y,y[1:-1,:],y[1:-1,:],y[1:-1,:],y[1:-1,:],y[1:-1,:],y),axis=0)
        self.um = np.concatenate((u,u[1:-1,:],u[1:-1,:],u[1:-1,:],u[1:-1,:],u[1:-1,:],u),axis=0)
        self.vm = np.concatenate((v,v[1:-1,:],v[1:-1,:],v[1:-1,:],v[1:-1,:],v[1:-1,:],v),axis=0)
        
    def path_ode(self,xy,t):
        x = xy[0]
        y = xy[1]
        
        # Interpolation to get the instantenous velocity
        xt = self.c*t+self.xm
        
        i = np.where((xt-x)>=0)[0][0]
        
        if self.ym[i-1,-1]<y:
            j = self.n_p-1
        else:
            j = np.where((self.ym[i-1,:]-y)>=0)[0][0]
            
        if self.ym[i,-1]<y:
            k = self.n_p-1
        else:
            k = np.where((self.ym[i,:]-y)>=0)[0][0]
        
        if j == 0 or j == self.n_p-1:
            j1 = j
        else:
            j1 = j-1
            
        if k == 0 or k == self.n_p-1:
            k1 = k
        else:
            k1 = k-1
        
        u1 = self.um[i-1,j1]
        u2 = self.um[i,k1]
        u3 = self.um[i-1,j]
        u4 = self.um[i,k]
        v1 = self.vm[i-1,j1]
        v2 = self.vm[i,k1]
        v3 = self.vm[i-1,j]
        v4 = self.vm[i,k]
        
        l1 = np.sqrt(np.power(xt[i-1]-x,2.)+np.power(self.ym[i-1,j1]-y,2.))
        l2 = np.sqrt(np.power(xt[i]-x,2.)+np.power(self.ym[i,k1]-y,2.))
        l3 = np.sqrt(np.power(xt[i-1]-x,2.)+np.power(self.ym[i-1,j]-y,2.))
        l4 = np.sqrt(np.power(xt[i]-x,2.)+np.power(self.ym[i-1,k]-y,2.))
        
        if np.abs(l1)<1e-6:
            u = u1
            v = v1
        elif np.abs(l2)<1e-6:
            u = u2
            v = v2
        elif np.abs(l3)<1e-6:
            u = u3
            v = v3
        elif np.abs(l4)<1e-6:
            u = u4
            v = v4
        else:
            l = 1./l1+1./l2+1./l3+1./l4;
            u = (u1/l1+u2/l2+u3/l3+u4/l4)/l
            v = (v1/l1+v2/l2+v3/l3+v4/l4)/l
            
        uv = [u,v]
        return uv
        
    
    def compute_path(self,xy0,t):
        return odeint(self.path_ode, xy0, t)
        
if __name__=="__main__":
    gamma = 0
    file_name_ext = 'gamma='+'{:.2f}'.format(gamma)
    
    flow1 = Flow("fix_p0.txt","FixMassFlux")
    flow1.load()
    
#    flow2 = Flow("fix_d.txt", "FixDepth")
#    flow2.load()