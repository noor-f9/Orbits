# -*- coding: utf-8 -*-
"""
Created : Tue Dec 20 17:56:25 2022
Modified: Sun Mar 31 09:22:30 2024

@author: Noor Alhasani (University of Southampton)

This is a 3d orbit visualisation tool for 
parabolic orbits (may include elliptic & hyperbolic orbits
in the future if there is the demand)

Please contact me at alnoor587@gmail.com for any 
suggestions/mistakes/questions about the code!

Enjoy!!!

"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import cross,sqrt,sin,cos,\
arccos,arctan2,radians,degrees,dot,pi,matmul
from numpy import linalg
import sys
from math import isnan

"Type in ' %matplotlib qt ' into the console to plot externally"
" recommended as the 3d plots are manipulable then! "

mu_body      = 398600 #km^3/s^2
body_radius  = 6378   #km
no_of_points = 1e4    #for orbit plot

input_mode = 'keplerian_elements'  #'state_vectors' or 'keplerian_elements'

#state vectors (position vector, velocity vector)
R_sc = [0, 0, 0]     # KM
V_sc = [0,-9,-3]   # KM/S

# keplerian orbital elements [a,e,i,Ω,ω,θ] below in KM & degrees
keplerian_elements = [26600,0.74,63.4,50,270,150]


plt.close()
    
def energy(R,V):
    r = linalg.norm(R)
    v = linalg.norm(V)
    ε = ((v**2)/2) - (mu_body/r)
    return ε

def kepler_to_state_vectors(element_list):
    #input a,e,i,Ω,ω,θ in km and degrees
    a,e,i,Ω,ω,θ = element_list
    i_rads,Ω_rads = radians(i),radians(Ω)
    ω_rads,θ_rads = radians(ω),radians(θ)
    
    h = sqrt(a*mu_body*(1-e**2)) #scalar ang. mom km^2/s
    r = (a*(1-e**2))/(1+e*cos(θ_rads))
    
    
    rotation_matrix = [
        [cos(Ω_rads)*cos(ω_rads)-(cos(i_rads)*sin(Ω_rads)*sin(ω_rads)),
         -(cos(Ω_rads)*sin(ω_rads))-(cos(i_rads)*cos(ω_rads)*sin(Ω_rads)),
         sin(Ω_rads)*sin(i_rads) ],
        [cos(ω_rads)*sin(Ω_rads) + cos(Ω_rads)*cos(i_rads)*sin(ω_rads),
         cos(Ω_rads)*cos(i_rads)*cos(ω_rads) - (sin(Ω_rads)*sin(ω_rads)),   
        -cos(Ω_rads)*sin(i_rads) ],
        [sin(i_rads)*sin(ω_rads),
         cos(ω_rads)*sin(i_rads),
         cos(i_rads)]  
        ]
    Q = rotation_matrix #matrix
    
    a1 = [[r*cos(θ_rads)],[r*sin(θ_rads)],[0]] #auxiliary matrix
    a2 = [[-sin(θ_rads)] ,[e+cos(θ_rads)],[0]] #auxiliary matrix
    
    r = dot(Q,a1)
    v = (mu_body/h) * dot(Q,a2)
    
    r = [float(r[0]),float(r[1]),float(r[2])] #turns from np array to list (for np cross function)
    v = [float(v[0]),float(v[1]),float(v[2])]
    
    return r,v #returns state_vectors
        
def state_vectors_to_kepler(R,V,print_elements=0,output='default'):
    
    r = linalg.norm(R) #instantaenous scalar radius
    #v = linalg.norm(V) #instantaneous scalar velocity
    
    H = cross(R,V) #vector
    h = linalg.norm(H) #returns magnitude of ang. mom
    
    try:
        E = (cross(V,H)/mu_body) - (R/r) #points to perigee
    except:
        print("Error computing eccentricity vector!\n")
        print("V = {} Km  H = {} Km^2/s".format(V,H))
        sys.exit()
        
    e = linalg.norm(E)
    
    
    if (e<0 or e>1):
        print("eccentricity exceeds preset limits!")
        sys.exit()
        
    a  = (h**2)/(mu_body*(1-e**2))
    
    #not strictly needed but may be useful
    ra = a*(1+e) #scalar in KM
    rp = a*(1-e) #scalar in KM
    RA = (-ra/e)*E #location of apogee
    RP = ( rp/e)*E
    
    if output=='rp,ra': #outputs perigee and apogee positions in 3d space
        return RP,RA
    
    
    try:
        print("\nnalpha\n")
        i_rads  = arccos(H[2]/h) #inclination 
    except:
        print("Error computing inclination!\n\n")
        print("H = [] km^2/s \n")
        print("h = {} km")
        raise Exception("??")

        #sys.exit()
        
    #as theres a pole for ω 
    i_rads  = i_rads + radians(1e-6)
    i_degs  = degrees(i_rads)
    if (i_degs < 0 or i_degs > 180):
        if (H[0] == 0) and (H[1]==0):#no inclination,orbit in x-y plane
            i_rads = 1e-6 #due to ω pole
            i_degs = degrees(i_rads) 
        else:
            print("Inclination is out of bounds! {}degs".format(i_degs))
            sys.exit()
        
    Ω_rads = arctan2(H[0],-H[1]) #no right ascention when i = 0
    if i_rads < 1e-5:
        Ω_rads = 0 #define RAAN as zero when inclination = 0
    Ω_degs = degrees(Ω_rads)
    if (Ω_degs < 0 or Ω_degs > 360):
        print("RAAN is out of bounds! {}degs".format(Ω_degs))
        sys.exit()
        
    arg1 = E[2]/sin(i_rads)
    arg2 = (E[1]*sin(Ω_rads)) + (E[0]*cos(Ω_rads))
    ω_rads = arctan2(arg1,arg2)
    #ω_rads = arccos(RP[1]*H[0] - RP[0]*H[1])
    ω_degs = degrees(ω_rads)
    if (ω_degs < 0 or ω_degs > 360): #attempts to fix omega before it exits
        if ω_degs < 0:
            ω_rads = ω_rads + 2*pi
            ω_degs = degrees(ω_rads)
        elif  ω_degs > 360:
            ω_rads = ω_rads - 2*pi
            ω_degs = degrees(ω_rads)
        
    if (ω_degs < 0 or ω_degs > 360): #if issue persists (shouldn't be more than 2 pi out)
        print("Arg. of Perigee is out of bounds! {}degs".format (ω_degs))
        sys.exit()
        
        
    if dot(R,V)>=0: #rising towards apogee
        θ_rads = arccos((dot(R,E))/(r*e))
    elif dot(R,V)<0: #descending towards perigee
        θ_rads = 2*pi - arccos((dot(R,E))/(r*e))
    θ_degs = degrees(θ_rads)
    if θ_degs < 0 or θ_degs > 360:
        print("True Anomaly is out of bounds!")
        sys.exit()
    if (R[0],R[1],R[2])==(RP[0],RP[1],RP[2]): #at perigee
        θ_rads = 1e-6
        θ_degs = degrees(θ_rads)
    if (R[0],R[1],R[2])==(RA[0],RA[1],RA[2]): #at apogee
        θ_rads = pi + 1e-6
        θ_degs = degrees(θ_rads)
    
        
    if print_elements == 1:
        
        print("\na = {} km".format(round(a,5)))
        print("e = {}".format(round(e,5)))
        print("i = {} degrees".format(round(i_degs,4)))
        print("Ω = {} degrees".format(Ω_degs))
        print("ω = {} degrees".format(ω_degs))
        print("θ = {} degrees\n".format(θ_degs))
        print("Perigee = {}km".format(RP))
        print("Apogee = {}km".format(RA))
        
    if output=='default':
        return a,e,i_rads,Ω_rads,ω_rads,θ_rads #in km and radians
    
    else:
        print("Select a valid output mode!")
        sys.exit()
    
def check_inputs():
    if input_mode == 'state_vectors':
        orbital_energy = energy(R_sc,V_sc)
        if orbital_energy > 0: #i.e. e not < 1
            print("Orbital energy is too high for \
    a circular/elliptic orbit! e = {}J/kg".format(round(orbital_energy,2)))
            sys.exit()
        else:   
            elements = state_vectors_to_kepler(R_sc,V_sc)
                   
    elif input_mode == 'keplerian_elements':
        i,Ω,ω,θ = keplerian_elements[2:]
        elements = keplerian_elements[0:2]+[radians(i),
        radians(Ω),radians(ω),radians(θ)]
        RR,VV = kepler_to_state_vectors(elements) #vectors
        rr,vv = linalg.norm(RR),linalg.norm(VV)   #scalars
        orbital_energy = energy(rr,vv)
        if orbital_energy >= 0: #i.e. e => 1
            print("Orbital energy is too high for \
    a circular/elliptic orbit! e = {}J/kg".format(round(orbital_energy,2)))
            sys.exit()
           
    else:
        print("Select 'state_vectors'\
     or 'keplerian_elements as an input mode!!")
        sys.exit()
    
    
def ellipse(a,e,n=int(no_of_points),plot='0',returnn='1'):
    b = a*sqrt(1-e**2)
    rp = a*(1-e)
    ra = a*(1+e)
    X = [i for i in np.linspace(-ra,rp,n)] #array
    Y_plus = [ b*sqrt( (1-((i+(a-rp))/a)**2)) for i in X]
    Y_minus =[-b*sqrt( (1-((i+(a-rp))/a)**2)) for i in X]
    if plot=='1':
        plt.plot(X,Y_plus,'k.-',label='S/C Orbit')
        plt.plot(X,Y_minus,'k.-')
    xx = X+X
    yy = Y_plus+Y_minus
    if returnn=='1':
        return xx,yy #return x,y points in perifocal frame
        
def perifocal_to_real_frame(xp,yp,zp,elements): #rotates frame 
    i_rads,Ω_rads,ω_rads,θ_rads = elements[2:] 
    R_perif = [[xp],[yp],[zp]] #matrix
    rotation_matrix = [
        [cos(Ω_rads)*cos(ω_rads)-(cos(i_rads)*sin(Ω_rads)*sin(ω_rads)),
         -(cos(Ω_rads)*sin(ω_rads))-(cos(i_rads)*cos(ω_rads)*sin(Ω_rads)),
         sin(Ω_rads)*sin(i_rads) ],
        
        [cos(ω_rads)*sin(Ω_rads) + cos(Ω_rads)*cos(i_rads)*sin(ω_rads),
         cos(Ω_rads)*cos(i_rads)*cos(ω_rads) - (sin(Ω_rads)*sin(ω_rads)),   
        -cos(Ω_rads)*sin(i_rads) ],
        
        [sin(i_rads)*sin(ω_rads),
         cos(ω_rads)*sin(i_rads),
         cos(i_rads)]  
        ]
    Q = rotation_matrix #matrix
    R = matmul(Q,R_perif)
    return R #returns x,y,z in global frame
    
def compute_inertial_points(input_data,input_mode='keplerian_elements'):
    
    inertial_x,inertial_y,inertial_z = [],[],[]
    
    if input_mode == 'keplerian_elements':
        elements = input_data
        
    elif input_mode == 'state_vectors':
        R,V = input_data[0:3],input_data[3:6]
        elements = state_vectors_to_kepler(R,V,print_elements=0,output='default')
        
    a,e = elements[0:2]
    
    #converts to radians for perif to real fcn
    elements[2] = radians(elements[2]) 
    elements[3] = radians(elements[3])
    elements[4] = radians(elements[4])
    elements[5] = radians(elements[5])
    
    perifocal_x,perifocal_y = ellipse(a,e)[0],ellipse(a,e)[1]
    perifocal_z = [0 for i in perifocal_x] #zero in 2d frame   
     
    for x,y,z in zip(perifocal_x,perifocal_y,perifocal_z):
        
        xi,yi,zi = perifocal_to_real_frame(x,y,z,elements)
        
        if (isnan(xi) == True or isnan(yi)==True or isnan(zi)==True):
            None
        else:
            inertial_x.append(float(xi))
            inertial_y.append(float(yi))
            inertial_z.append(float(zi))
 
    return inertial_x,inertial_y,inertial_z


#Plotting Orbit(s)

scale = 3

xmin,xmax = -scale*body_radius,scale*body_radius
ymin,ymax = -scale*body_radius,scale*body_radius
zmin,zmax = -scale*body_radius,scale*body_radius

plot_lims = 1.1*min([xmin,ymin,zmin]),1.1*max(([xmax,ymax,zmax]))

fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(projection='3d')

ax.set_xlim(plot_lims)
ax.set_ylim(plot_lims)
ax.set_zlim(plot_lims)

ax.set_xlabel('X (Km)'),ax.set_ylabel('Y (Km)'),ax.set_zlabel('Z (Km)')

#Body
rr = body_radius
density = 20j
u, v = np.mgrid[0:2*np.pi:density, 0:np.pi:density]
x = rr*np.cos(u) * np.sin(v)
y = rr*np.sin(u) * np.sin(v)
z = rr*np.cos(v)
ax.plot_wireframe(x,y,z,alpha=0.5)

#[a,e,i,Ω,ω,θ] below in KM & degrees

orbits = [[26600, 0.2, 20, 0, 0, 0],\
          [26600, 0.2, 40, 0, 0, 0]]
colours = ['r','g','b']
n = 1
for orbit in orbits:
     
    inertial_x,inertial_y,inertial_z = \
    compute_inertial_points(orbit,'keplerian_elements')
    
    ax.plot(inertial_x,inertial_y,inertial_z,\
    color=colours[n-1],markersize=1,label='Orbit {}'.format(n))
        
    n = n + 1
    

ax.legend()

#Inertial Axes

text_size = 15
length = 15e3 #size of inertial earth axes

ax.quiver(0,0,0,1,0,0,length=length,color='r')
ax.text(length,0,0,'$X_e$',color='k',fontsize=text_size)

ax.quiver(0,0,0,0,1,0,length=length,color='g')
ax.text(0,length,0,'$Y_e$',color='k',fontsize=text_size)

ax.quiver(0,0,0,0,0,1,length=length,color='b')
ax.text(0,0,length,'$Z_e$',color='k',fontsize=text_size)
