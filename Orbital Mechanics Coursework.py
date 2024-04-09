# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 07:05:27 2024

@author: Noor Alhasani
"""

from numpy import radians,sqrt,sin,pi,sin,cos,tan,arctan,e,arange,\
linspace,arange,degrees
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

plt.close()

#Earth Gravitational Parameter [Km^3/s^2]
mu = 398600 

#Earth nominal radius [Km]
re = 6378

#Initial parking orbit radius [Km]
r_initial = re+450 

#Geostationary Orbit Radius [Km]
r_geo = 42164

#Orbital inclination change [rad]
delta_i = radians(58.5107)

#Initial spacecraft mass [Kg]
m_initial = 1500 

#Spacecraft Thruster Exhaust Velocity [Km/s]
vx = (250 * 9.81) / 1000

#------------------------------------------------------------------------------------

def compute_resultant_dv(dv1,dv2):
    return sqrt((dv1**2)+(dv2**2))

def compute_ang_momentum(r1,r2):
    return sqrt(2*mu*((r1*r2)/(r1+r2)))

def compute_ecc(r1,r2):
    if r1 > r2:
        ra,rp = r1,r2
    else:
        rp,ra = r1,r2
    return (ra-rp)/(ra+rp)
    
def compute_inclination_change_dv(delta,r1,r2,location='apogee'):
    if r1 > r2:
        ra,rp = r1,r2
    else:
        rp,ra = r1,r2
    h  = compute_ang_momentum(rp,ra)
    if location == 'apogee':    #minimises inclination cost
        vt = h / ra
    elif location == 'perigee': #sub-optimal (dV wise)
        vt = h / rp
    else:
        raise Exception('Input perigee or apogee for burn location!')
    return 2*vt*sin(delta/2)

def compute_semi_major_axis(r1,r2):
    return (r1+r2)/2
    
def compute_full_orbital_period(r1,r2):
    a = compute_semi_major_axis(r1,r2)
    return 2*pi*sqrt((a**3)/(mu))

def compute_circular_velocity(r_const):
    return sqrt(mu/r_const)

def compute_perigee_velocity(r1,r2):
    if r1 > r2:
        ra,rp = r1,r2
    else:
        rp,ra = r1,r2
    return sqrt(((2*mu)/rp)-((2*mu)/(ra+rp)))
    
def compute_apogee_velocity(r1,r2):
    if r1 > r2:
        ra,rp = r1,r2
    else:
        rp,ra = r1,r2
    return sqrt(((2*mu)/ra)-((2*mu)/(ra+rp)))

def compute_fuel_mass_percentage(dv):
    return (1 - (e**(-dv/vx)))*100

def compute_burn_dv(v1,i_1,v2,i_2=0): 
    i_1,i_2 = radians(i_1),radians(i_2)
    dvx = v2*cos(i_2) - v1*cos(i_1)
    dvz = v2*sin(i_2) - v1*sin(i_1)
    return compute_resultant_dv(dvx,dvz)
#-----------------------------------------------------------------------------
def compute_total_bielliptic_dv(r_transfer):
    
    #Boost Apogee far out beyond GEO
    dv1      = compute_perigee_velocity(r_initial, r_transfer) - \
               compute_perigee_velocity(r_initial, r_initial)
    
    #Rotate plane and boost slightly
    v1 = compute_apogee_velocity(r_transfer, r_initial)
    v2 = compute_apogee_velocity(r_transfer, r_geo)
    dv2 = compute_burn_dv(v1,degrees(delta_i),v2,i_2=0)
    
    
    #Circularise at perigee to enter final GEO orbit
    dv3      = compute_perigee_velocity(r_transfer, r_geo) - \
               compute_perigee_velocity(r_geo,      r_geo)
          
    return dv1 + dv2 + dv3

def compute_total_transfer_time(r_transfer):
    
    #first rising to apogee of the bi-elliptic orbit
    t1 = 0.5*compute_full_orbital_period(r_initial,r_transfer)
    
    #now falling towards perigee (which is now at geo height)
    t2 = 0.5*compute_full_orbital_period(r_geo,r_transfer)
    
    return t1 + t2

radii              = arange(r_geo,800e3,1e3)
dvs                = [compute_total_bielliptic_dv(r) for r in radii]
transfer_time      = [compute_total_transfer_time(r) for r in radii]
transfer_time_days = [i/(3600*24) for i in transfer_time]

plt.style.use('default')

fig, ax1 = plt.subplots(figsize=(10, 4))  
plt.style.use('default')

def with_commas(x, pos):
    return '{:,}'.format(int(x))

ax1.xaxis.set_major_formatter(FuncFormatter(with_commas))

color = 'tab:blue'
ax1.set_xlabel('Bi-Elliptic Transfer Apogee [km]')
ax1.set_ylabel('Total ΔV\n[km/s]', rotation=0, labelpad=30, color=color)
ax1.plot(radii, dvs, color=color)  
ax1.tick_params(axis='y', labelcolor=color)
plt.xticks(linspace(0,800e3,9))
ax1.grid()

ax1.plot(42164,5.006,'ro')
ax1.text(50000,5,'Hohmann\nTransfer',color='r')

ax2 = ax1.twinx()
color = 'tab:green'
ax2.set_ylabel('Total\nTransfer\nTime\n[Days]', rotation=0, labelpad=30, color=color)
ax2.plot(radii, transfer_time_days, color=color)
ax2.tick_params(axis='y', labelcolor=color)

ax1_y_ticks = linspace(4.5, 5.1, 7)  
ax2_y_ticks = arange(1.83, 34, 5.0167) 

ax1.set_yticks(ax1_y_ticks)
ax2.set_yticks(ax2_y_ticks)


plt.title('Total ΔV & Total Transfer time v.s. Bi-elliptic Apogee Radius')
plt.show()

def true_anom_to_time(rp,ra,theta1,theta2):  #return answer in hours
    theta1,theta2 = radians(theta1),radians(theta2)
    e  = compute_ecc(rp,ra)
    h  = compute_ang_momentum(rp,ra)
    #Auxiliary Terms----------------------------------
    K  = (h**3)/((mu**2)*(1-e**2)**1.5)
    A1 = 2*arctan( sqrt((1-e)/(1+e))*(tan(theta1/2)) )
    B1 = (e*sqrt(1-e**2)*sin(theta1))/(1+e*cos(theta1))
    A2 = 2*arctan( sqrt((1-e)/(1+e))*(tan(theta2/2)) )
    B2 = (e*sqrt(1-e**2)*sin(theta2))/(1+e*cos(theta2))
    #-------------------------------------------------
    t1 = K*(A1-B1)
    t2 = K*(A2-B2)
    return (t2-t1)/3600