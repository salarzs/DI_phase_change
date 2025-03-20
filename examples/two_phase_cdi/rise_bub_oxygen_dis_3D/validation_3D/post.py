#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 00:14:59 2023
@author: salarzamanisalimi
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy as sc
from scipy import special
import glob
import pickle
import os
import sys as sys
import math
from numpy import diff
from numpy import array
from numpy.linalg import norm
from pylab import *
from matplotlib import rc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


sns.set()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
sys.setrecursionlimit(10000)


mpl.rcParams['xtick.major.size']  = 8
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['xtick.minor.size']  = 4
mpl.rcParams['xtick.minor.width'] = 3
#
mpl.rcParams['ytick.major.size']  = 8
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['ytick.minor.size']  = 4
mpl.rcParams['ytick.minor.width'] = 3

#
plt.rcParams['xtick.labelsize']   = 18
plt.rcParams['ytick.labelsize']   = 18
#
#
mpl.rcParams['axes.linewidth']    = 3
mpl.rcParams['axes.labelsize']    = 18


#############PLEASE SPECIFY THE DIRECTORY OF THE BINARY OUTPUTS HERE

direct     =   '/cluster/work/users/salarzs/DI_incompressible/incom_conserv/data_3D/'

##################


def loadfile(fn,n):
   f   = open(fn,'rb')
   fld = np.fromfile(f,dtype='f8',count=np.product(n))
   data = np.reshape(fld,(n[0],n[1],n[2]),order='f')
   return data

###
###
dx = 0.008/192

location = np.linspace(1,192,192)*dx/1e-3 ##convert to milimeter

N = [96,192,96]
#############################################################
####################################################
adress_p = str(direct)+"psi_fld_"

fig, ax_main = plt.subplots(figsize=(9, 7))
time = 0
dt   = 1e-6
for x in range(0,180000,30000):
    
   x   = '{:09d}'.format(x)
   
   adr_p = adress_p+str(x)+".bin"


   psi     = loadfile(adr_p,N)
   
   
   sns.lineplot(x=location,y=psi[48,:,48],marker='',linestyle='-',markersize=0,linewidth
                 ='2',ax=ax_main,label='t='+str(time*1e-6))
   time     = (time + 30000.0)
   
ax_main.set_title('Time history of bubble position', fontsize=16)
ax_main.set_xlabel('vertial-location [mm]', fontsize=18)
ax_main.set_ylabel(r'$C(t)$' , fontsize=18)

#ax_main.legend(loc='lower left')
ax_main.grid(True)

plt.savefig('conc_water_his.png',dpi=900)
plt.show()


##############################################################
######################################################
###  HISTORY OF OXYGEN IN THE WATER
######################################################
adress_c2 = str(direct)+"c_2_fld_"

fig, ax_main = plt.subplots(figsize=(9, 7))
time = 0
dt   = 1e-6
for x in range(0,180000,30000):
    
   x   = '{:09d}'.format(x)
   
   adr_c2 = adress_c2+str(x)+".bin"


   c2     = loadfile(adr_c2,N)
   
   
   sns.lineplot(x=location,y=c2[48,:,48],marker='',linestyle='-',markersize=0,linewidth
                 ='2',ax=ax_main,label='t='+str(time*1e-6))
   time     = (time + 30000.0)
   
ax_main.set_title('Time history of oxygen profile in the water', fontsize=16)
ax_main.set_xlabel('vertial-location [mm]', fontsize=18)
ax_main.set_ylabel(r'$C(t)$' , fontsize=18)
#ax_main.legend(loc='lower left')
ax_main.grid(True)

plt.savefig('conc_water_his.png',dpi=900)
plt.show()

######################################################
###  HISTORY OF OXYGEN IN THE BUBBLE
######################################################
adress_c1 = str(direct)+"c_1_fld_"

fig, ax_main = plt.subplots(figsize=(9, 7))
time = 0
dt   = 1e-6
for x in range(0,180000,30000):
    
   x   = '{:09d}'.format(x)
   
   adr_c1 = adress_c1+str(x)+".bin"


   c1     = loadfile(adr_c1,N)
   
   
   sns.lineplot(x=location,y=c1[48,:,48],marker='',linestyle='-',markersize=0,linewidth
                 ='2',ax=ax_main,label='t='+str(time*1e-6))
   time     = (time + 30000.0)
   
ax_main.set_title('Time history of oxygen profile inside the bubble', fontsize=16)
ax_main.set_xlabel('vertial-location [mm]', fontsize=18)
ax_main.set_ylabel(r'$C(t)$' , fontsize=18)

#ax_main.legend(loc='lower left')
ax_main.grid(True)

plt.savefig('conc_bub_his.png',dpi=900)
plt.show()


################################################
###  RISING VELOCITY 

adress_c1 = str(direct)+"c_1_fld_"
adress_v  = str(direct)+"vey_fld_"
adress_p  = str(direct)+"psi_fld_"


fig, ax_main = plt.subplots(figsize=(9, 7))

time = 0

for x in range(0,150000,10000):
    
   x   = '{:09d}'.format(x)
   
   adr_c1   = adress_c1+str(x)+".bin"
   adr_p    = adress_p +str(x)+".bin"
   adr_v    = adress_v +str(x)+".bin"

   psi = loadfile(adr_p,N)
   vey = loadfile(adr_v,N)
   c1  = loadfile(adr_c1,N)

   
   int_nom = sum(psi*vey)
   vol = sum(psi)
   
   rise_vel = int_nom/vol
   dt       = 1e-6
   tau      = 0.01899
   
   sns.scatterplot(x=[time*dt/tau],y=[rise_vel],marker='o',color='black',s=50,ax=ax_main)
   
   time     = (time + 10000.0)
   

ax_main.set_xlabel(r'$t/\tau$', fontsize=18)
ax_main.set_ylabel(r'$V$' , fontsize=18)

ax_main.grid(True)

plt.savefig('rise_vel.png',dpi=900)
plt.show()

##################################################
######################CONCENTRATION IN tHE BUBBLE
#Bothe data
#bothe = np.loadtxt('bothe.in') 

#X_val = bothe[:,0]
#Y_val = bothe[:,1]


adress_c1 = str(direct)+"c_1_fld_"

adress_p  = str(direct)+"psi_fld_"


fig, ax_main = plt.subplots(figsize=(9, 7))

time = 0

for x in range(0,100000,10000):

   x   = '{:09d}'.format(x)
   
   adr_c1   = adress_c1+str(x)+".bin"
   adr_p    = adress_p +str(x)+".bin"


   psi = loadfile(adr_p,N)

   c1  = loadfile(adr_c1,N)

   
   int_nom = sum(psi*c1)
   vol = sum(psi)
   
   concent = int_nom/vol
   
   c0      = 0.7606
   dt      = 1e-6
   tau     = 0.01899
   
   sns.scatterplot(x=[time*dt/tau],y=[concent/c0],marker='o',color='black',s=100,ax=ax_main)
   
   
   time     = (time + 10000.0)


ax_main.set_xlabel(r'$t/\tau$', fontsize=18)
ax_main.set_ylabel(r'$C/C(0)$' , fontsize=18)

legend_labels = [r"DNS-192$\times384$", r"Bothe and Fleckenstein(2013)-$512\times256$"]

legend_handles = [
    plt.Line2D([0], [0], marker='o', linestyle='' ,color='black', label=legend_labels[0], markersize=10)
]

legend=plt.legend(handles=legend_handles, loc='lower left')

for label in legend.get_texts():
    label.set_fontsize(16)  # Adjust the font size as neede
    
    
plt.ylim(0.8,1)
plt.xlim(0,5)
plt.savefig('bubble_concent.png',dpi=900)
plt.show()

##############################################

