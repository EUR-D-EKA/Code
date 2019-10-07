#! /usr/bin/python
#/* Copyright (c) 2018 Siddhartha Shelton */

import os
import subprocess
from subprocess import call
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math as mt
import matplotlib.patches as mpatches



def orbital_speeds(): #creates a plot of the orbits of the planets in the solar system
    speeds = [47.4, 35.0, 29.8, 24.1,13.1, 9.7, 6.8,5.4, 4.7]#km/s
    radii  = [0.39, 0.723, 1, 1.524, 5.203, 9.539, 19.18, 30.06, 39.53] #AU
    planet_dia = [4879., 12104., 12756., 6792., 142984., 120536., 51118., 49528., 2370.] #km
    colors = ['gray', 'green', 'mediumblue', 'red', 'orange', 'yellow', 'cyan', 'navy', 'black' ]
    plt.figure(figsize=(20, 10))
    plt.tick_params(axis='y', which='major', labelsize=24)
    plt.tick_params(axis='x', which='major', labelsize=24)
    plt.ylabel('Orbital Velocity (km/s)', fontsize=26)
    plt.xlabel('Distance From Sun (AU)', fontsize=26)
    plt.plot(radii, speeds, linestyle = '-', linewidth = 3, color ='k', alpha = 1)
    
    for i in range(0, len(radii)):
        sz = mt.log10(planet_dia[i]**4.0 / 16.0)
        
        print(sz)
        plt.scatter(radii[i], speeds[i], marker = 'o', s = 0.001 * sz**5, color = colors[i])
    plt.savefig('orbital_vels.png', format='png')
    
    
    

def bulge_pot(r):
    d = 0.7 #kpc
    M = 1.52954402e5# kpc^3/gy^2  
    a = r * M / (r * (r + d)**2.0)
    return a#kpc/gy^2


def disk_pot(r):
    M = 4.45865888e5  # kpc^3/gy^2  
    b = 6.5 #kpc
    c = 0.26 #kpc
    a = r * M / (r * r + (b + c)**2.)**(3./2.)
    return a#kpc/gy^2

def halo_pot(r):
    v0 = 73.0 / (0.977813107) #km/s # factor for converting the v0's km/s to kpc/gy
    d = 12.0 #kpc
    a = 2. * v0 * v0 * r / (r * r + d * d)
    return a #kpc/gy^2
    
def mw_galaxy_orbital_speeds():
    w_halo = []
    no_halo = []
    bulge = []
    disk = []
    halo = []
    rs = [] 
    r = 0.01
    esc = 1.0 #mt.sqrt(2.0 * (bulge_pot(8.) + disk_pot(8.) + halo_pot(8.)) )
    for i in range(1000):
        bulge_a = bulge_pot(r)
        disk_a  = disk_pot(r)
        halo_a  = halo_pot(r)
        
        #each component separately
        bul = mt.sqrt(r * (bulge_a )) # in kpc/gy. Multiply by 0.977813107  to get km/s
        dsk = mt.sqrt(r * (disk_a)) 
        hlo = mt.sqrt(r * (halo_a)) 
        
        nh  =  mt.sqrt(r * (bulge_a + disk_a)) #MW with no halo
        h   =  mt.sqrt(r * (bulge_a + disk_a + halo_a)) #MW with halo
        
        bulge.append(bul / esc)
        disk.append(dsk / esc)
        halo.append(hlo / esc)
        w_halo.append(h / esc)
        no_halo.append(nh / esc)
        rs.append(r)
        
        r += 0.1
        
    params = {'legend.fontsize': 22,
            'legend.handlelength': 1}
    plt.rcParams.update(params)    
    plt.figure(figsize=(13, 10))
    plt.tick_params(axis='y', which='major', labelsize=24)
    plt.tick_params(axis='x', which='major', labelsize=24)
    plt.ylabel('Orbital Velocity (kpc/Gyr)', fontsize=26)
    plt.xlabel('Distance From GC (kpc)', fontsize=26)
    plt.xlim(0,100 )
    plt.plot(rs, bulge, linestyle = ':', linewidth = 1, color ='green', alpha = 1, label = 'Central Bulge')
    plt.plot(rs, disk, linestyle = '-.', linewidth = 1, color ='cyan', alpha = 1, label = 'Disk')
    plt.plot(rs, halo, linestyle = '--', linewidth = 1, color ='red', alpha = 1, label = 'Dark Matter halo')
    
    plt.plot(rs, w_halo, linestyle = '-', linewidth = 3, color ='k', alpha = 1, label = 'with halo')
    plt.plot(rs, no_halo, linestyle = '-', linewidth = 3, color ='b', alpha = 1, label = 'without halo')
    plt.legend()
    plt.show()
        

def main():
    #orbital_speeds()
    mw_galaxy_orbital_speeds()
    
main()
