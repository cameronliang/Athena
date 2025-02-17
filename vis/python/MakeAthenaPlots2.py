import matplotlib
matplotlib.use('Agg')

import yt
import sys
sys.settrace
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import *
import struct
import array
import os
from scipy.interpolate import griddata

import h5py

from mpi4py import MPI

# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

crat=5694.76

rhomin=1.e-3
rhomax=0.05

Ermin=10.0
Ermax=2.e3

def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,time,outputname,vel=0):
    plots, axes = plt.subplots(figsize=(10.,8.4),dpi=300)
    plt.xlabel('$ x/r_g$', size = 30)
    plt.ylabel('$ y/r_g$', size = 30)
    plt.subplots_adjust(left=0.16,right=0.78,top=0.85,bottom=0.11)
    plt.title(time,size=25,y=1.12)
    
    if vel>0:
      speed=np.sqrt(vx_cart**2+vy_cart**2)
      if vlim2 < vlim1:
         vlim2=speed.max()
      speed=np.clip(speed,vlim1,vlim2)
      logspeed=np.log10(speed)
      vcolor=axes.imshow(speed,cmap='jet',norm=LogNorm(vmin=vlim1,vmax=vlim2),origin='lower', extent=[xmin,xmax,ymin,ymax])
      velaxes = plots.add_axes([0.24, 0.89, 0.5, 0.03])
      velcbar = plots.colorbar(vcolor, cax = velaxes, orientation = 'horizontal')
      velcbar.ax.tick_params(labelsize=20)
      if speed.max() > vlim1:
        velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1,arrowsize=3.0,color=logspeed)
      else:
        velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1,arrowsize=3.0)
        
      axes.set_ylim([ymin,ymax])

      #RdGy_r 
    im = axes.imshow(data_cart,cmap='RdYlBu_r', norm=LogNorm(vmin = minval, vmax=maxval), \
                     origin='lower', extent=[xmin,xmax,ymin,ymax])
                     
    cbaxes = plots.add_axes([0.8,0.15,0.03,0.6])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(label1, size=30)
    cbar.ax.tick_params(labelsize=25)
#    axes.set_xticks([-30.,-20.0, -10., 0., 10.,20.,30.0])
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    
    axes.set_aspect('auto')
    plt.savefig(outputname)
    plt.close(plots)


def makemovies(ni,no,minval,maxval,vlim1,vlim2,width=25,res=4,vel=0,var1='density',var2='vel1',var3='vel3',label1='$\\rho/\\rho_0$',label2='$v/c$'):
  
    for i in range(ni,no+1):
        print i
        filename='disk.out1.'+'{:05d}'.format(i)+'.athdf'
        data=yt.load(filename)

#get region slice
        reg=data.r[:,0.5*np.pi,:]
        rhoslice=reg[var1].d
        rcoord=reg['r'].d
        phi=reg['phi'].d

        vel1=reg[var2].d
        vel2=reg[var3].d

        vx=vel1*np.cos(phi)-vel2*np.sin(phi)
        vy=vel1*np.sin(phi)+vel2*np.cos(phi)

# convert to cartesian coordinate
        xcoord=rcoord*np.cos(phi)
        ycoord=rcoord*np.sin(phi)

#create grid for plotting
        rootgrid=data.domain_dimensions
        nr=rootgrid[0]*res
        ntheta=rootgrid[1]*res
        nx=nr
        ny=2*nr
        xmin=-width
# xmax=np.max(rcoord)
        xmax=width
        ymin=-xmax
        ymax=xmax
        xgrid=np.linspace(xmin, xmax, nx)
        ygrid=np.linspace(ymin, ymax, ny)

        xmesh,ymesh=np.meshgrid(xgrid,ygrid)

        rmesh=(xmesh**2.0+ymesh**2.0)**0.5

        rindix=rmesh < 2.0

        rho_cart=griddata(np.c_[xcoord,ycoord],rhoslice,(xmesh,ymesh),method='nearest')

        vx_cart=griddata(np.c_[xcoord,ycoord],vx,(xmesh,ymesh),method='nearest')
        vy_cart=griddata(np.c_[xcoord,ycoord],vy,(xmesh,ymesh),method='nearest')

        rho_cart[rindix] = minval
        vx_cart[rindix] = 0.0
        vy_cart[rindix] = 0.0

        xgrid = xgrid*2
        ygrid = ygrid*2
        xmin = xmin*2
        xmax = xmax*2
        ymin = ymin*2  
        ymax = ymax*2

# only scale to speed of light for velocity
        if var2=='vel1':
           vx_cart=vx_cart/crat
           vy_cart=vy_cart/crat

        outputname='disk_z.'+'{:05d}'.format(i)+'_'+var1+'.png'

        labelname='$\\rho/\\rho_0$'

        time='${\\rm time}='+"%4.2f"%(data.current_time*crat)+'{\ r_s/c}$'

        MakeRhoVSlice(rho_cart, vx_cart, vy_cart, minval,maxval, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,time,outputname,vel)

#################################################


# The main program

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

ni=10988
no=14567

for i in range(ni,no+1,nprocs):
  fi=i+rank
  print fi, rank
  makemovies(fi,fi,rhomin,rhomax,1.e-1,1.0,vel=1)
  makemovies(fi,fi,Ermin,Ermax,1.e-3,10.0,vel=1,var1='Er',var2='B1',var3='B2',label1='$E_r/a_rT_0^4$',label2='B')


