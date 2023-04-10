import numpy as np
from numpy import sin,cos,sqrt,arctan2,pi,floor
import matplotlib.pyplot as plt
    
#  %matplotlib qt <- put this into console in spyder to plot externally

plt.close()
        
G = 6.67430e-11 #Universal Gravitational Constant

M1 = 6e24   #Mass of body 1 in kg
M2 = 7e22   #Mass of body 2 in kg
M3 = 2e30  #Mass of body 3 in kg


Time_of_simulation = 24*366 #IN HOURS 
deltaT = 60 #Time physics delta in seconds (smaller = more accurate and more computationall expensive)

deltaT_hours = deltaT / 3600

#initial positions and velocities of 3 bodies
R1 = [[1.5e11,0,0]]
V1 = [[0,30e3,0]]               

R2 = [[1.504e11,0,0]]
V2 = [[0,31e3,0]]             

R3 = [[0,0,0]] 
V3 = [[0,0,0]]                

threshold = 1 * 1e3 #if bodies come closer than this value, terminate the simulation

iterations = Time_of_simulation * 3600 / deltaT #max no. iterations to be performed

#below block is arrays of x,y,z distance differences between bodies 1,2,3. Eg: delX_12 is x distances between bodies 1 & 2
delX_12, delY_12, delZ_12 = [R2[0][0]-R1[0][0]] , [R2[0][1]-R1[0][1]] , [R2[0][2]-R1[0][2]]
delX_13, delY_13, delZ_13 = [R3[0][0]-R1[0][0]] , [R3[0][1]-R1[0][1]] , [R3[0][2]-R1[0][2]] 
delX_23, delY_23, delZ_23 = [R3[0][0]-R2[0][0]] , [R3[0][1]-R2[0][1]] , [R3[0][2]-R2[0][2]]

#absolute distances (3d pythagoras basically)
distances_12 = []
distances_13 = []
distances_23 = []

distances_12.append(sqrt(delX_12[-1]**2 + delY_12[-1]**2 + delZ_12[-1]**2))
distances_13.append(sqrt(delX_13[-1]**2 + delY_13[-1]**2 + delZ_13[-1]**2))
distances_23.append(sqrt(delX_23[-1]**2 + delY_23[-1]**2 + delZ_23[-1]**2))

#copies x,y,z difference arrays as acrtan2 function seems to modify input arrays (I don't wanna modify orig. array)
delX_12_copy, delY_12_copy, delZ_12_copy = [i for i in delX_12],[i for i in delY_12],[i for i in delZ_12]
delX_13_copy, delY_13_copy, delZ_13_copy = [i for i in delX_13],[i for i in delY_13],[i for i in delZ_13]
delX_23_copy, delY_23_copy, delZ_23_copy = [i for i in delX_23],[i for i in delY_23],[i for i in delZ_23]

#In 3D, a distance and two angles are needed to fully describe relative positions of bodies.
#theta denotes azimuth angle (sideways, +ve counter-clockwise from positive x-axis)
#phi denotes elevation angle (upwards, +ve upwards from positive x-axis to positive z-axis)
#subscript 12 means looking from body 1 towards body 2
thetas12,thetas21,thetas13,thetas31,thetas23,thetas32 = [],[],[],[],[],[]
phis12,phis21,phis13,phis31,phis23,phis32 = [],[],[],[],[],[]

theta12 = float(arctan2(delY_12_copy,delX_12_copy)) 
phi12 = float(arctan2(delZ_12_copy, np.sqrt(delX_12_copy[-1]**2+delY_12_copy[-1]**2))) 
theta21 = float(pi + theta12) 
phi21 = float(-1 * phi12) 

thetas12.append(theta12)
phis12.append(phi12)
thetas21.append(theta21)
phis21.append(phi21)

theta13 = float(arctan2(delY_13_copy,delX_13_copy))
phi13 = float(arctan2(delZ_13_copy, np.sqrt(delX_13_copy[-1]**2+delY_13_copy[-1]**2)))
theta31 = float(pi + theta13)
phi31 = float(-1 * phi13)

thetas13.append(theta13)
phis13.append(phi13)
thetas31.append(theta31)
phis31.append(phi31)

theta23 = float(arctan2(delY_23_copy,delX_23_copy))
phi23 = float(arctan2(delZ_23_copy, np.sqrt(delX_23_copy[-1]**2+delY_23_copy[-1]**2)))
theta32 = float(pi + theta23)
phi32 = float(-1 * phi23)

thetas23.append(theta23)
phis23.append(phi23)
thetas32.append(theta32)
phis32.append(phi32)

#Initialising force vectors for each body
F12 = [(G*M1*M2)/(delX_12[-1]**2+delY_12[-1]**2+delZ_12[-1]**2)]
F13 = [(G*M1*M3)/(delX_13[-1]**2+delY_13[-1]**2+delZ_13[-1]**2)]
F23 = [(G*M2*M3)/(delX_23[-1]**2+delY_23[-1]**2+delZ_23[-1]**2)]

#Initialising accelerations vectors for body 1
A12 = [[(F12[-1]/M1)*cos(phi12)*cos(theta12),(F12[-1]/M1)*cos(phi12)*sin(theta12),(F12[-1]/M1)*sin(phi12)]]#Acceleration on body 1 due to body 2
A13 = [[(F13[-1]/M1)*cos(phi13)*cos(theta13),(F13[-1]/M1)*cos(phi13)*sin(theta13),(F13[-1]/M1)*sin(phi13)]]#Acceleration on body 1 due to body 3,
A1  = [[A12[-1][0]+A13[-1][0], A12[-1][1]+A13[-1][1], A12[-1][2]+A13[-1][2]]] #Resultant acceleration on body 1 due to body 2 & body 3\

#Initialising accelerations vectors for body 2
A21 = [[(F12[-1]/M2)*cos(phi21)*cos(theta21),(F12[-1]/M2)*cos(phi21)*sin(theta21),(F12[-1]/M2)*sin(phi21)]]
A23 = [[(F23[-1]/M2)*cos(phi23)*cos(theta23),(F23[-1]/M2)*cos(phi23)*sin(theta23),(F23[-1]/M2)*sin(phi23)]]
A2  = [[A21[-1][0]+A23[-1][0], A21[-1][1]+A23[-1][1], A21[-1][2]+A23[-1][2]]] #Resultant acceleration of body 2 due to body 1 and body 3

#Initialising accelerations vectors for body 3
A31 = [[(F13[-1]/M3)*cos(phi31)*cos(theta31),(F13[-1]/M3)*cos(phi31)*sin(theta31),(F13[-1]/M3)*sin(phi31)]]
A32 = [[(F23[-1]/M3)*cos(phi32)*cos(theta32),(F23[-1]/M3)*cos(phi32)*sin(theta32),(F23[-1]/M3)*sin(phi32)]]
A3  = [[A31[-1][0]+A32[-1][0], A31[-1][1]+A32[-1][1], A31[-1][2]+A32[-1][2]]] #Resultant acceleration of body 2 due to body 1 and body 3

V1.append( [ V1[-1][0]+A1[-1][0]*deltaT, V1[-1][1]+A1[-1][1]*deltaT, V1[-1][2]+A1[-1][2]*deltaT ] ) 
V2.append( [ V2[-1][0]+A2[-1][0]*deltaT, V2[-1][1]+A2[-1][1]*deltaT, V2[-1][2]+A2[-1][2]*deltaT ] )
V3.append( [ V3[-1][0]+A3[-1][0]*deltaT, V3[-1][1]+A3[-1][1]*deltaT, V3[-1][2]+A3[-1][2]*deltaT ] )

R1.append( [ R1[-1][0]+V1[-1][0]*deltaT, R1[-1][1]+V1[-1][1]*deltaT, R1[-1][2]+V1[-1][2]*deltaT ] )
R2.append( [ R2[-1][0]+V2[-1][0]*deltaT, R2[-1][1]+V2[-1][1]*deltaT, R2[-1][2]+V2[-1][2]*deltaT ] )
R3.append( [ R3[-1][0]+V3[-1][0]*deltaT, R3[-1][1]+V3[-1][1]*deltaT, R3[-1][2]+V3[-1][2]*deltaT ] )

simulation_termination_loop = 0
simulation_terminated = 0
counter = 0

while counter < iterations: #runs the iterative computational algorithm

    delX_12.append( float(R2[-1][0] - R1[-1][0]) ) 
    delY_12.append( float(R2[-1][1] - R1[-1][1]) )
    delZ_12.append( float(R2[-1][2] - R1[-1][2]) )

    delX_13.append( float(R3[-1][0] - R1[-1][0]) )
    delY_13.append( float(R3[-1][1] - R1[-1][1]) )
    delZ_13.append( float(R3[-1][2] - R1[-1][2]) )

    delX_23.append( float(R3[-1][0] - R2[-1][0]) )
    delY_23.append( float(R3[-1][1] - R2[-1][1]) )
    delZ_23.append( float(R3[-1][2] - R2[-1][2]) ) 

    distances_12.append(sqrt(delX_12[-1]**2 + delY_12[-1]**2 + delZ_12[-1]**2))
    distances_13.append(sqrt(delX_13[-1]**2 + delY_13[-1]**2 + delZ_13[-1]**2))
    distances_23.append(sqrt(delX_23[-1]**2 + delY_23[-1]**2 + delZ_23[-1]**2))

    if (distances_12[-1] < threshold):
        print("\nSimulation has been stopped!!!!\n")
        print("\nClosest point of approach is between bodies 1 & 2, distance is {} km\n".format(round((distances_12[-1]/1000),4)))
        simulation_terminated = 1
        simulation_termination_loop = counter
        counter = iterations #basically means this is the last loop, effectively terminates the program

    elif (distances_13[-1] < threshold):
        print("\nSimulation has been stopped!!!!\n")
        print("\nClosest point of approach is between bodies 1 & 3, distance is {} km\n".format(round((distances_13[-1]/1000),4)))
        simulation_terminated = 2
        simulation_termination_loop = counter
        counter = iterations #basically means this is the last loop, effectively terminates the program

    elif (distances_23[-1] < threshold):
        print("\nSimulation has been stopped!!!!\n")
        print("\nClosest point of approach is between bodies 2 & 3, distance is {} km\n".format(round((distances_23[-1]/1000),4)))
        simulation_terminated = 3
        simulation_termination_loop = counter
        counter = iterations #basically means this is the last loop, effectively terminates the program

    delX_12_copy.append(delX_12[-1]) 
    delY_12_copy.append(delY_12[-1])
    delZ_12_copy.append(delZ_12[-1])

    delX_13_copy.append(delX_13[-1])
    delY_13_copy.append(delY_13[-1])
    delZ_13_copy.append(delZ_13[-1])

    delX_23_copy.append(delX_23[-1])
    delY_23_copy.append(delY_23[-1])
    delZ_23_copy.append(delZ_23[-1])

    theta12 = float(arctan2(delY_12_copy[-1],delX_12_copy[-1])) 
    phi12 = float(arctan2(delZ_12_copy[-1], sqrt(delX_12_copy[-1]**2+delY_12_copy[-1]**2))) 
    theta21 = float(pi + theta12) 
    phi21 = float(-1 * phi12) 

    thetas12.append(theta12)
    phis12.append(phi12)
    thetas21.append(theta21)
    phis21.append(phi21)

    theta13 = float(arctan2(delY_13_copy[-1],delX_13_copy[-1]))
    phi13 = float(arctan2(delZ_13_copy[-1], sqrt(delX_13_copy[-1]**2+delY_13_copy[-1]**2)))
    theta31 = float(pi + theta13)
    phi31 = float(-1 * phi13)

    thetas13.append(theta13)
    phis13.append(phi13)
    thetas31.append(theta31)
    phis31.append(phi31)

    theta23 = float(arctan2(delY_23_copy[-1],delX_23_copy[-1]))
    phi23 = float(arctan2(delZ_23_copy[-1], sqrt(delX_23_copy[-1]**2+delY_23_copy[-1]**2)))
    theta32 = float(pi + theta23)
    phi32 = float(-1 * phi23)

    thetas23.append(theta23)
    phis23.append(phi23)
    thetas32.append(theta32)
    phis32.append(phi32)

    #below forces are magnitudes of gravitational forces, these need to later be resolved into X,Y,Z
    F12.append((G*M1*M2)/(delX_12[-1]**2+delY_12[-1]**2+delZ_12[-1]**2))
    F13.append((G*M1*M3)/(delX_13[-1]**2+delY_13[-1]**2+delZ_13[-1]**2))
    F23.append((G*M2*M3)/(delX_23[-1]**2+delY_23[-1]**2+delZ_23[-1]**2))

    #A12 is the acceleration on body 1 due to body 2, A13 is the acceleration on body 1 due to body 3
    A12.append([(F12[-1]/M1)*cos(phi12)*cos(theta12),(F12[-1]/M1)*cos(phi12)*sin(theta12),(F12[-1]/M1)*sin(phi12)])
    A13.append([(F13[-1]/M1)*cos(phi13)*cos(theta13),(F13[-1]/M1)*cos(phi13)*sin(theta13),(F13[-1]/M1)*sin(phi13)])
    A1.append( [A12[-1][0]+A13[-1][0], A12[-1][1]+A13[-1][1], A12[-1][2]+A13[-1][2]] ) #Resultant acceleration on body 1 due to body 2 & body 3

    A21.append([(F12[-1]/M2)*cos(phi21)*cos(theta21),(F12[-1]/M2)*cos(phi21)*sin(theta21),(F12[-1]/M2)*sin(phi21)])
    A23.append([(F23[-1]/M2)*cos(phi23)*cos(theta23),(F23[-1]/M2)*cos(phi23)*sin(theta23),(F23[-1]/M2)*sin(phi23)])
    A2.append( [A21[-1][0]+A23[-1][0], A21[-1][1]+A23[-1][1], A21[-1][2]+A23[-1][2]]) #Resultant acceleration of body 2 due to body 1 and body 3

    A31.append([(F13[-1]/M3)*cos(phi31)*cos(theta31),(F13[-1]/M3)*cos(phi31)*sin(theta31),(F13[-1]/M3)*sin(phi31)])
    A32.append([(F23[-1]/M3)*cos(phi32)*cos(theta32),(F23[-1]/M3)*cos(phi32)*sin(theta32),(F23[-1]/M3)*sin(phi32)])
    A3.append( [A31[-1][0]+A32[-1][0], A31[-1][1]+A32[-1][1], A31[-1][2]+A32[-1][2]] )

    V1.append( [ V1[-1][0]+A1[-1][0]*deltaT, V1[-1][1]+A1[-1][1]*deltaT, V1[-1][2]+A1[-1][2]*deltaT ] ) #Velocity vector for body 1
    V2.append( [ V2[-1][0]+A2[-1][0]*deltaT, V2[-1][1]+A2[-1][1]*deltaT, V2[-1][2]+A2[-1][2]*deltaT ] )
    V3.append( [ V3[-1][0]+A3[-1][0]*deltaT, V3[-1][1]+A3[-1][1]*deltaT, V3[-1][2]+A3[-1][2]*deltaT ] )

    R1.append( [ R1[-1][0]+V1[-1][0]*deltaT, R1[-1][1]+V1[-1][1]*deltaT, R1[-1][2]+V1[-1][2]*deltaT ] ) #Position vector for body 1
    R2.append( [ R2[-1][0]+V2[-1][0]*deltaT, R2[-1][1]+V2[-1][1]*deltaT, R2[-1][2]+V2[-1][2]*deltaT ] )
    R3.append( [ R3[-1][0]+V3[-1][0]*deltaT, R3[-1][1]+V3[-1][1]*deltaT, R3[-1][2]+V3[-1][2]*deltaT ] )

    counter = counter + 1
    
#computes closest distance between each pair
closest_12 = min(distances_12)
closest_13 = min(distances_13)
closest_23 = min(distances_23)

#divides lists by 1E3 to get numbers from M into KM
x1,y1,z1 = [i[0]/1000 for i in R1], [i[1]/1000 for i in R1], [i[2]/1000 for i in R1]
x2,y2,z2 = [i[0]/1000 for i in R2], [i[1]/1000 for i in R2], [i[2]/1000 for i in R2]
x3,y3,z3 = [i[0]/1000 for i in R3], [i[1]/1000 for i in R3], [i[2]/1000 for i in R3]

scale_down = 1e3 #animates 1 in n positions
fps = 100 #max is 100

if (len(x1) < scale_down): #n is the scale down factor (animates in 1 in n frames)
    n = 1
else:
    n = int(scale_down)

#Shortens lists by interpolating between values to make animation easier, also converts to numpy lists
x1,y1,z1 = x1[::n],y1[::n],z1[::n]
x2,y2,z2 = x2[::n],y2[::n],z2[::n]
x3,y3,z3 = x3[::n],y3[::n],z3[::n]
    
width = 12 #how wide to make plot window
height = 12 #how high to make plot window

line_transparency = 1 #transparency of lines in plot
line_thickness = 2 #thickness of line in plots

visual_size_of_body1 = 10 #how fat the dots in the animation are
visual_size_of_body2 = 5
visual_size_of_body3 = 20

if simulation_terminated != 0:
    Time_of_simulation_actual = (simulation_termination_loop / iterations) * Time_of_simulation
    if simulation_terminated == 1:
        print("The simulation has been terminated prematurely after {} hours as body 1 & 2 come within {} km of eachother."
              .format(round(Time_of_simulation_actual,2),round(distances_12[-1]/1000,3)))
    if simulation_terminated == 2:
        print("The simulation has been terminated prematurely after {} hours as body 1 & 3 come within {} km of eachother."
              .format(round(Time_of_simulation_actual,2),round(distances_13[-1]/1000,3)))
    if simulation_terminated == 3:
        print("The simulation has been terminated prematurely after {} hours as body 2 & 3 come within {} km of eachother."
              .format(round(Time_of_simulation_actual,2),round(distances_23[-1]/1000,3)))
else: 
    Time_of_simulation_actual = Time_of_simulation
    
simulation_time_per_frame = deltaT_hours*n #n is animation factor - hours

fig = plt.figure(0)
plt.style.use("dark_background") #to give it a bit of space snazz
ax = fig.add_subplot(projection='3d')
fig.set_size_inches(width,height)
ax.set_title("3-Body Numerical Simulation ({} hours) - Inertial Reference Frame".format(round(Time_of_simulation_actual,2)))
ax.set_xlabel('X (Km)'),ax.set_ylabel('Y (Km)'),ax.set_zlabel('Z (km)')

autoscale = 'no'

if autoscale == 'yes':
    None #graph automatically scales during animation
else:                  #manually set physical size of animation space
    ax.set_xlim(-3e11,3e11)
    ax.set_ylim(-3e11,3e11)
    ax.set_zlim(-3e11,3e11)

while(1): #loops animation

    i = 0
    
    path_1x, path_1y, path_1z= [], [], []
    path_2x, path_2y, path_2z= [], [], []
    path_3x, path_3y, path_3z= [], [], []
    
    for data in x1:
        
        path_1x.append(x1[int(i)]), path_1y.append(y1[int(i)]), path_1z.append(z1[int(i)])
        path_2x.append(x2[int(i)]), path_2y.append(y2[int(i)]), path_2z.append(z2[int(i)])
        path_3x.append(x3[int(i)]), path_3y.append(y3[int(i)]), path_3z.append(z3[int(i)])
        
        point1, = ax.plot(x1[i], y1[i], z1[i], 'bo',label='Body 1 ({} kg)'.format(round(M1,5)),markersize=visual_size_of_body1)
        point2, = ax.plot(x2[i], y2[i], z2[i], 'ro',label='Body 2 ({} kg)'.format(round(M2,5)),markersize=visual_size_of_body2)
        point3, = ax.plot(x3[i], y3[i], z3[i], 'yo',label='Body 3 ({} kg)'.format(round(M3,5)),markersize=visual_size_of_body3)
        
        line1, = ax.plot(path_1x, path_1y, path_1z, 'b-',linewidth=line_thickness,alpha=line_transparency) 
        line2, = ax.plot(path_2x, path_2y, path_2z, 'r-',linewidth=line_thickness,alpha=line_transparency)
        line3, = ax.plot(path_3x, path_3y, path_3z, 'y-',linewidth=line_thickness,alpha=line_transparency) 
        
        ax.legend(bbox_to_anchor=(0.1, 1))
        text = ax.text(1,1,1e3,'Time = {} hours ({} days)'\
                    .format(round(i*simulation_time_per_frame,2),\
                    round(i*(simulation_time_per_frame)/24,2)),color='k')
        
        plt.pause(1/fps)
        
        point1.remove(),line1.remove()
        point2.remove(),line2.remove()
        point3.remove(),line3.remove()
        text.remove()
    
        i = i + 1