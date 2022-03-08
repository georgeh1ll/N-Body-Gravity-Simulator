
# =============================================================================
# Final Code for N Body 
# =============================================================================

#Importing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import LSODA #not directly used, but can be used to compare how fast LSODA solves compared to RK methods

#%%


G = 6.67430e-11 #Gravitational constant 

## Sun inital conditions ##

x_sun_inital=0 #Sun x coord
y_sun_inital=0 #Sun y coord
z_sun_inital=0 #Sun z coord
vx_sun_inital=0  #Sun velocity in x-direction
vy_sun_inital=0  #Sun velocity in y-direction
vz_sun_inital=0 #Sun velocity in z-direction
M_s=1.989e30 #Sun mass in kg

## Earth inital conditions ##

x_earth_inital= 1.496*10**11 #Earth x coord - 1AU initally 
y_earth_inital=0 #Earth y coord
z_earth_inital=0 #Earth z coord
vx_earth_inital=0 #Earth velocity in x-direction 
vy_earth_inital=np.sqrt((G*M_s)/x_earth_inital) #Earth velocity in y-direction 
vz_earth_inital=0 #Earth velocity in z-direction 
M_e=5.972*10**24 #Earth mass in kg

## Time the System evolves over ##

year = 3.154*10**7 #Year in seconds
ti=0 #Inital time
tf=5*year #Solves up to  5 years
t=np.arange(ti,tf,10)


#Defining 2D system of Earth and Sun
def solving_system_earth(System_Earth,t):
    
    #Defining a 2D system of all variables to solve at any time t
    x_earth,y_earth,x_sun,y_sun,vx_earth,vy_earth,vx_sun,vy_sun  = System_Earth 
    r_se=np.sqrt((x_sun-x_earth)**2 +(y_sun-y_earth)**2) #Radius vector Sun - Earth
    
    return [vx_earth,
            vy_earth,
            vx_sun,
            vy_sun,
            (G*M_s/r_se**3) *(x_sun-x_earth),
            (G*M_s/r_se**3) *(y_sun-y_earth), 
            (G*M_e/r_se**3) * (x_earth-x_sun),
            (G*M_e/r_se**3) *(y_earth-y_sun)]


#Solving 2D System of Earth and Sun 
Solution_2D_Earth = odeint(solving_system_earth, y0=[x_earth_inital, y_earth_inital, x_sun_inital, y_sun_inital,
                                         vx_earth_inital,vy_earth_inital,
                                         vx_sun_inital,vy_sun_inital],
                                         t=t)
Solution_2D_Earth = Solution_2D_Earth/1.496e11 #Converting solution into AU
t1=Solution_2D_Earth.T[0] #time


#%%
# Plotting distance from sun against time (test plot)

fig1=plt.figure(1,figsize=(10,10))
axsec=plt.gca() #gets current axis
axsec.plot((Solution_2D_Earth.T[0]))
axsec.tick_params(labelsize=15) #Increasing tick size

plt.xlabel("Time (Seconds)",fontsize=18)
plt.ylabel("Distance from the Sun in AU",fontsize=18)
plt.title("$x$⨁ against time over 5 years",fontsize=24,x=0.5,y=1.1)

#Adding year axis
axyears=axsec.twiny()
axyears.set_xticks([0,1,2,3,4,5])
axyears.set_xlabel("Time (Years)",fontsize=18)
axyears.tick_params(labelsize=15) #making ticks readable size

plt.show()

#%%
# Plotting full orbit view (test plot 2)

fig2=plt.figure(2,figsize=(12,12))
x_earth_sol= Solution_2D_Earth[:,0] #x coord of Earth
y_earth_sol= Solution_2D_Earth[:,1] #y coord of Earth
x_sun_sol= Solution_2D_Earth[:,2] #x coord of the Sun
y_sun_sol= Solution_2D_Earth[:,3] #y coord of the Sun

plt.plot(x_earth_sol,y_earth_sol,'b') #Plotting Earth's orbit
plt.plot(x_sun_sol,y_sun_sol,'orange',linewidth=5) #Plotting the Sun's orbit
plt.title("Earth's Orbit around the Sun",fontsize=24)
plt.xlabel('$x$' r'$\bigoplus$',fontsize=18)
plt.ylabel('$y$' r'$\bigoplus$',fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.show()

#%%
## 3D Plotting of Earth around the Sun
fig3= plt.figure(3,figsize=(10,10))
ax3=plt.axes(projection='3d') #3d axis setup
plt.plot(x_earth_sol,y_earth_sol,0,linewidth=5) #Plotting Earth Sun orbit with no z components.
plt.plot(x_sun_sol,y_sun_sol,0,linewidth=5)
plt.title("Earth Orbit around Sun 3D Axis",fontsize=20)
plt.xlabel('$x$' r'$\bigoplus$',fontsize=16)
plt.ylabel('$y$' r'$\bigoplus$',fontsize=16)
ax3.set_zlabel('$z$' r'$\bigoplus$',fontsize=16)

ax3.locator_params(nbins=6) #6 ticks on each axis for no overlapping
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax3.zaxis.set_tick_params(labelsize=14)
ax3.set_aspect('auto') #auto selects best aspect ratio to display 
plt.show()
#%%
## Attempting with Mars ##

## Mars Inital Conditions ##

x_mars_inital= 1.5*1.496e11 #x coord of Mars in AU
y_mars_inital=0 #y coord of Mars
z_mars_inital=0 #Z coord of Mars
vx_mars_inital= 0 #Velocity of Mars in x component
vy_mars_inital= np.sqrt((G*M_s)/x_mars_inital) #Velocity of Mars in y component 
vz_mars_inital=0 #Velocity of Mars in z component
M_m= 6.39e23 #Mar's mass in kg


##Defining Mars Sun Problem ##

def evolving_system_mars(System_Mars,t):

    #Defining a 2D system of all variables to solve at any time tm
    x_mars,y_mars,x_sun,y_sun,vx_mars,vy_mars,vx_sun,vy_sun = System_Mars
    r_ms= np.sqrt((x_sun-x_mars)**2 +(y_sun-y_mars)**2)
    
    
    return [vx_mars,
            vy_mars,
            vx_sun,
            vy_sun,
            (G*M_m/r_ms**3)*(x_sun-x_mars),
            (G*M_m/r_ms**3) *(y_sun-y_mars), 
            (G*M_m/r_ms**3) * (x_mars-x_sun),
            (G*M_m/r_ms**3) *(y_mars-y_sun)]


#Solving Mars Sun problem 

Solution_Mars = odeint(evolving_system_mars, y0=[x_mars_inital, y_mars_inital ,
                                                 x_sun_inital,y_sun_inital,
                                         vx_mars_inital,vy_mars_inital,
                                         vx_sun_inital,vy_sun_inital,],
                                         t=t)
Solution_Mars = Solution_Mars/1.496e11 #Converting solution into AU

x_mars_sol= Solution_Mars[:,0] #x coord of Mars
y_mars_sol= Solution_Mars[:,1] #y coord of Mars



#Solving Mars 2D system
def evolving_system_mars(System_Mars,t):
    
    #Defining a 2D system of all variables to solve at any time t
    x_mars,y_mars,x_sun,y_sun,vx_mars,vy_mars,vx_sun,vy_sun  = System_Mars
    rm=np.sqrt((x_sun-x_mars)**2 +(y_sun-y_mars)**2) #Radius vector
    
    return [vx_mars,
            vy_mars,
            vx_sun,
            vy_sun,
            (G*M_s/rm**3) *(x_sun-x_mars),
            (G*M_s/rm**3) *(y_sun-y_mars), 
            (G*M_m/rm**3) * (x_mars-x_sun),
            (G*M_m/rm**3) *(y_mars-y_sun)]



Solution_2D_Mars = odeint(evolving_system_mars, y0=[x_mars_inital, y_mars_inital ,
                                                    x_sun_inital,y_sun_inital,
                                         vx_mars_inital,vy_mars_inital,
                                         vx_sun_inital,vy_sun_inital],
                                         t=t)
Solution_2D_Mars = Solution_2D_Mars/1.496e11 #Converting solution into AU

x_mars_sol= Solution_2D_Mars[:,0] #x coord of Earth
y_mars_sol= Solution_2D_Mars[:,1] #y coord of Earth
x_sun_sol= Solution_2D_Mars[:,2] #x coord of the Sun
y_sun_sol= Solution_2D_Mars[:,3] #y coord of the Sun



## 3D Plotting of Earth, Mars, Sun orbit.
fig4= plt.figure(4,figsize=(10,10))
ax4=plt.axes(projection='3d')
plt.plot(x_mars_sol,y_mars_sol,0,label="Mars Orbit",color='Red') #plots x,y coords of mars
plt.title("Earth and Mars Orbit 3D",fontsize=20)
plt.plot(x_earth_sol,y_earth_sol,color='blue',label="Earth Orbit")
plt.plot(x_sun_sol,y_sun_sol,0,label="Sun Orbit",color='orange',linewidth=4) #Plotting Mars Sun orbit with no z components.
plt.xlabel('$x$' r'$\bigoplus$',fontsize=16)
plt.ylabel('$y$' r'$\bigoplus$',fontsize=16)
ax4.set_zlabel('$z$' r'$\bigoplus$',fontsize=16)
plt.show()


#%%
# =============================================================================
# 2 Heavy Stars and 1 Smaller Mass 
# =============================================================================

#Setting inital conditions

#Inital masses
M_e=5.972e24 
M_Star1=1e50
M_Star2=1e35
M_Planet=1e20
G=6.6743e-11

#Inital positions
x_star1_inital = 1e10
y_star1_inital = 0
z_star1_inital = 0
x_star2_inital=2e10
y_star2_inital = 1e10
z_star2_inital =0
x_planet_inital =-2e10
y_planet_inital =-2e10
z_planet_inital = 0

#Inital radius vectors
r_s1_s2= np.sqrt((x_star2_inital-x_star1_inital)**2)
r_s1_p3= np.sqrt((x_planet_inital-x_star2_inital)**2)
r_s2_p3 = np.sqrt((x_planet_inital-x_star1_inital)**2)

#Inital velocites
vx_star1_inital =0
vy_star1_inital = np.sqrt(G*M_Star2/np.abs(r_s1_s2))+np.sqrt(G*M_Planet/np.abs(r_s1_p3))
vz_star1_inital = 0
vx_star2_inital = 0
vy_star2_inital = np.sqrt(G*M_Star1/np.abs(r_s1_s2))+np.sqrt(G*M_Planet/np.abs(r_s2_p3))
vz_star2_inital=0
vx_planet_inital = 0
vy_planet_inital = np.sqrt(G*M_Star1/np.abs(r_s1_p3))+np.sqrt(G*M_Star2/np.abs(r_s2_p3))
vz_planet_inital = 0


#Defining three body systems with 2 stars, 1 planet 
def three_body_2stars(t, System_2stars):
    x_star1,y_star1,z_star1,x_star2,y_star2,z_star2,x_planet,y_planet,z_planet,vx_star1,vy_star1,vz_star1,vx_star2, vy_star2,vz_star2,vx_planet,vy_planet,vz_planet = System_2stars
    r_s1_s2 = np.sqrt((x_star2-x_star1)**2 + (y_star2-y_star1)**2 + (z_star2-z_star1)**2)
    r_s1_p3 = np.sqrt((x_planet-x_star1)**2 + (y_planet-y_star1)**2 +(z_planet-z_star1)**2)
    r_s2_p3 = np.sqrt((x_star2-x_planet)**2 + (y_star2-y_planet)**2 + (z_star2-z_planet)**2)
    return [ vx_star1,
            vy_star1,
            vz_star1,
            vx_star2,
            vy_star2,
            vz_star2,
            vx_planet,
            vy_planet,
            vz_planet,
            G*M_Star2/r_s1_s2**3 * (x_star2-x_star1) + M_Planet/r_s1_p3**3 * (x_planet-x_star1), #Star1
            G*M_Star2/r_s1_s2**3 * (y_star2-y_star1) + M_Planet/r_s1_p3**3 * (y_planet-y_star1),
            G*M_Star2/r_s1_s2**3 * (z_star2-z_star1)+ M_Planet/r_s1_p3**3 *(z_planet-z_star1),
            G*M_Star1/r_s1_s2**3 * (x_star1-x_star2) + M_Planet/r_s2_p3**3 * (x_planet-x_star2), #Star2
            G*M_Star1/r_s1_s2**3 * (y_star1-y_star2) + M_Planet/r_s2_p3**3 * (y_planet-y_star2),
            G*M_Star1/r_s1_s2**3 * (z_star1-z_star2) +M_Planet/r_s2_p3**3 * (z_planet-z_star2),
            G*M_Star1/r_s1_p3**3 * (x_star1-x_planet) + M_Star2/r_s2_p3**3 * (x_star2-x_planet), #Planet
            G*M_Star1/r_s1_p3**3 * (y_star1-y_planet) + M_Star2/r_s2_p3**3 * (y_star2-y_planet),
            G*M_Star1/r_s1_p3**3 *(z_star1-z_planet) + M_Star2/r_s2_p3**3 *(z_star2-z_planet)]


#time system runs over
t_min=0
t_max=1000
t = np.linspace(t_min, t_max, 100000)


#Solving three body system of 2 stars, 1 planet 

Solution_3_Body_2_Stars= solve_ivp(three_body_2stars,y0=[x_star1_inital,
                                                         y_star1_inital,
                                                         z_star1_inital,
                                          x_star2_inital, y_star2_inital ,
                                          z_star2_inital,
                                          x_planet_inital, y_planet_inital,
                                          z_planet_inital,
                       vx_star1_inital, vy_star1_inital,vz_star1_inital, 
                       vx_star2_inital, vy_star2_inital,vz_star2_inital,
                       vx_planet_inital, vy_planet_inital,vz_planet_inital],
                     method='RK23', t_span=(0,1000))

#coordinates of each object over time 
x_star1_sol = Solution_3_Body_2_Stars.y[0]
y_star1_sol = Solution_3_Body_2_Stars.y[1]
z_star1_sol = Solution_3_Body_2_Stars.y[2]
x_star2_sol = Solution_3_Body_2_Stars.y[3]
y_star2_sol = Solution_3_Body_2_Stars.y[4]
z_star2_sol = Solution_3_Body_2_Stars.y[5]
x_planet_sol = Solution_3_Body_2_Stars.y[6]
y_planet_sol = Solution_3_Body_2_Stars.y[7]
z_planet_sol = Solution_3_Body_2_Stars.y[8]

t = Solution_3_Body_2_Stars.t

#Animates the three body system by plotting positions to line objects 
def animate_2stars_1planet(i):
    line1.set_data([x_star1_sol[i]], [y_star1_sol[i]])
    line2.set_data([x_star2_sol[i],y_star2_sol[i]])
    line3.set_data([x_planet_sol[i],y_planet_sol[i]])
    
fig5=plt.figure(figsize=(12,12))
ax5=plt.axes()
ax5.set_facecolor('black') #background black for space theme
plt.grid() #adds grid to plot background 
#Plotting positions
line1, = plt.plot([], [],'r*', lw=3, markersize=20,label="Star1")
line2, =plt.plot([],[],'b*',lw=3,label="Star2",markersize=20)
line3, = plt.plot([],[],'go',label="Planet",markersize=10)

#Axis labelling
plt.xlabel("$x$(metres)",fontsize=18)
plt.ylabel("$y$(metres)",fontsize=18)
plt.xlim(-10e10,10e10)
plt.ylim(-10e10,10e10)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend()
plt.title("2 Stars and a Planet Orbit",fontsize=22)
#blit false for three body systems 
ani1 = animation.FuncAnimation(fig5, animate_2stars_1planet,
                              frames=1000, interval=1,blit=False)
plt.show()


#%%
#3D plotting coordinates over time for 2 star, one planet 
fig9=plt.figure(figsize=(22,14))
plt.axis('off')
plt.title("Coordinates Plotted in 3D over Time",fontsize=26)
#Setting up 3 subplots wit 3D axes
ax9=fig9.add_subplot(1,3,1,projection='3d')
ax10=fig9.add_subplot(1,3,2,projection='3d')
ax11=fig9.add_subplot(1,3,3,projection='3d')
plt.subplots_adjust(hspace=0,wspace=0.3,left=0,right=None)

#Plotting star 1 coords
#labelpad used so axes ticks and axes labels do not overlap 
ax9.plot(x_star1_sol,y_star1_sol,z_star1_sol,color='r')
ax9.set_xlabel("   X Coordinate (10^10 metres)",fontsize=18,labelpad=30)
ax9.set_ylabel("   Y Coordinate (10^10 metres)",fontsize=18,labelpad=30)
ax9.set_zlabel("Z Coordinate (metres)",fontsize=18,labelpad=30)
ax9.set_title("Coordinates of Star 1",fontsize=22)
ax9.tick_params(axis='both',labelsize=16,pad=10)

#Plotting star 2 coords
ax10.plot(x_star2_sol,y_star2_sol,z_star2_sol,color='b')
ax10.set_xlabel("   X Coordinate (10^13 metres)",fontsize=18,labelpad=30)
ax10.set_ylabel("   Y Coordinate (10^13 metres)",fontsize=18,labelpad=30)
ax10.set_zlabel("Z Coordinate (metres)",fontsize=18,labelpad=30)
ax10.set_title("Coordinates of Star 2",fontsize=22)
ax10.tick_params(axis='both',labelsize=14,pad=10)

#Plotting planet coords 
ax11.plot(x_planet_sol,y_planet_sol,z_planet_sol,color='g')
ax11.tick_params(axis='both',labelsize=14,pad=10)
ax11.set_xlabel("   X Coordinate (10^13 metres)",fontsize=18,labelpad=30)
ax11.set_ylabel("   Y Coordinate (10^13 metres)",fontsize=18,labelpad=30)
ax11.set_zlabel("Z Coordinate (metres)",fontsize=18,labelpad=30)
ax11.set_title("Coordinates of Planet",fontsize=22)
plt.show()

#%%

# =============================================================================
# 2 Planets and 1 Star
# =============================================================================

#Setting inital conditions

#masses inital conditions
M_Star1=1e50
M_Planet1=1e20
M_Planet2=1e20
G=6.6743e-11

#positions inital conditions
x_star_inital = 1e10
y_star_inital = 0
z_star_inital = 0
x_planet1_inital=10e10
y_planet1_inital = 10e10
z_planet1_inital =0
x_planet2_inital =-10e10
y_planet2_inital =-10e10
z_planet2_inital = 0

#inital radius vectors 
r_p1_s= np.sqrt((x_planet1_inital-x_star_inital)**2)
r_p2_p1= np.sqrt((x_planet2_inital-x_planet1_inital)**2)
r_p2_s = np.sqrt((x_planet2_inital-x_star_inital)**2)

#inital velocities 
vx_star_inital =0
vy_star_inital = np.sqrt(G*M_Planet2/np.abs(r_p1_s))+np.sqrt(G*M_Planet1/np.abs(r_p2_p1))
vz_star_inital = 0
vx_planet1_inital = 0
vy_planet1_inital = np.sqrt(G*M_Star1/np.abs(r_p1_s))+np.sqrt(G*M_Planet1/np.abs(r_p2_s))
vz_planet1_inital=0
vx_planet2_inital = 0
vy_planet2_inital = np.sqrt(G*M_Star1/np.abs(r_p2_p1))+np.sqrt(G*M_Planet2/np.abs(r_p2_s))
vz_planet2_inital = 0


#defining three body system with 1 star, 2 planets 
def three_body_1star_2planets(t, System_1star_2planets):
    x_star,y_star,z_star,x_planet1,y_planet1,z_planet1,x_planet2,y_planet2,z_planet2,vx_star,vy_star,vz_star,vx_planet1,vy_planet1,vz_planet1,vx_planet2,vy_planet2,vz_planet2 = System_1star_2planets
    r_p1_s = np.sqrt((x_planet1-x_star)**2 + (y_planet1-y_star)**2 + (z_planet1-z_star)**2)
    r_p2_p1 = np.sqrt((x_planet2-x_star)**2 + (y_planet2-y_star)**2 +(z_planet2-z_star)**2)
    r_p2_s = np.sqrt((x_planet1-x_planet2)**2 + (y_planet1-y_planet2)**2 + (z_planet1-z_planet2)**2)
    return [ vx_star,
            vy_star,
            vz_star,
            vx_planet1,
            vy_planet1,
            vz_planet1,
            vx_planet2,
            vy_planet2,
            vz_planet2,
            G*M_Planet2/r_p1_s**3 * (x_planet1-x_star) + M_Planet1/r_p2_p1**3 * (x_planet2-x_star), #Star
            G*M_Planet2/r_p1_s**3 * (y_planet1-y_star) + M_Planet1/r_p2_p1**3 * (y_planet2-y_star),
            G*M_Planet2/r_p1_s**3 * (z_planet1-z_star)+ M_Planet1/r_p2_p1**3 *(z_planet2-z_star),
            G*M_Star1/r_p1_s**3 * (x_star-x_planet1) + M_Planet1/r_p2_s**3 * (x_planet2-x_planet1), #Planet 1
            G*M_Star1/r_p1_s**3 * (y_star-y_planet1) + M_Planet1/r_p2_s**3 * (y_planet2-y_planet1),
            G*M_Star1/r_p1_s**3 * (z_star-z_planet1) +M_Planet1/r_p2_s**3 * (z_planet2-z_planet1),
            G*M_Star1/r_p2_p1**3 * (x_star-x_planet2) + M_Planet2/r_p2_s**3 * (x_planet1-x_planet2), #Planet 2
            G*M_Star1/r_p2_p1**3 * (y_star-y_planet2) + M_Planet2/r_p2_s**3 * (y_planet1-y_planet2),
            G*M_Star1/r_p2_p1**3 *(z_star-z_planet2) + M_Planet2/r_p2_s**3 *(z_planet1-z_planet2)]

#time to evolve over 
t_min=0
t_max=1000
t = np.linspace(t_min, t_max, 100000)


#solving three body system of 1 star, 2 planets 
Solution_3_Body_1_Star_2_Planets= solve_ivp(three_body_1star_2planets,
                                   y0=[x_star_inital, y_star_inital,
                                       z_star_inital,
                                          x_planet1_inital, y_planet1_inital ,
                                          z_planet1_inital,
                                          x_planet2_inital, y_planet2_inital,
                                          z_planet2_inital,
                       vx_star_inital, vy_star_inital,vz_star_inital, 
                       vx_planet1_inital, vy_planet1_inital,vz_planet1_inital,
                       vx_planet2_inital, vy_planet2_inital,vz_planet2_inital],
                     method='RK45', t_span=(0,1000))


#inidividual positions component solutions for each object
x_star_sol = Solution_3_Body_1_Star_2_Planets.y[0]
y_star_sol = Solution_3_Body_1_Star_2_Planets.y[1]
z_star_sol = Solution_3_Body_1_Star_2_Planets.y[2]
x_planet1_sol = Solution_3_Body_1_Star_2_Planets.y[3]
y_planet1_sol = Solution_3_Body_1_Star_2_Planets.y[4]
z_planet1_sol = Solution_3_Body_1_Star_2_Planets.y[5]
x_planet2_sol = Solution_3_Body_1_Star_2_Planets.y[6]
y_planet2_sol = Solution_3_Body_1_Star_2_Planets.y[7]
z_planet2_sol = Solution_3_Body_1_Star_2_Planets.y[8]

t = Solution_3_Body_1_Star_2_Planets.t

#animinating 3 lines with positions being plotted 
def animate_1star_2planets(i):
    line4.set_data([x_star_sol[i]], [y_star_sol[i]])
    line5.set_data([x_planet1_sol[i],y_planet1_sol[i]])
    line6.set_data([x_planet2_sol[i],y_planet2_sol[i]])
    
fig5=plt.figure(figsize=(12,12))
ax5=plt.axes()
ax5.set_facecolor('black')
plt.grid()
#plot the line data 
line4, = plt.plot([], [],'r*', lw=3, markersize=20,label="Star")
line5, =plt.plot([],[],'bo',lw=3,label="Planet 1",markersize=10)
line6, = plt.plot([],[],'go',lw=3,label="Planet 2",markersize=10)

plt.xlabel("$x$(metres)",fontsize=16)
plt.ylabel("$y$(metres)",fontsize=16)
plt.xlim(-50e10,50e10)
plt.ylim(-50e10,50e10)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend()
plt.title("Star and 2 Planets Orbit",fontsize=22)
#blit = false for three body problems 
ani2 = animation.FuncAnimation(fig5, animate_1star_2planets,
                              frames=1000, interval=1,blit=False)
plt.show()


#%%

# =============================================================================
# N Body Functions
# =============================================================================

#function to define general forces and break into components 
def Forces_Function(i,k):
    r_vector=np.sqrt(((x_pos[i]-x_pos[k])**2)+(y_pos[i]+y_pos[k]**2)+
              ((z_pos[i]+z_pos[k])**2)) 
    #r seperation disatnce for each point
    force_general= -G*masses[i]*masses[k]/(r_vector**2) #Newton's Law for gravitational force
    
    #Splitting force into three components of motion
    force_y=force_general*(y_pos[i]-y_pos[k])/r_vector
    force_x=force_general*(x_pos[i]-x_pos[k])/r_vector
    force_z=force_general*(z_pos[i]-z_pos[k])/r_vector
    
    return force_x,force_y,force_z #return x,y,z components of force 

#function to get the resultant force for each object based on all other objects. Then seperates into components.
def Resultant_Forces_Function(i):
    #setting all inital resultant forces to zero 
    force_x_res,force_y_res,force_z_res=0,0,0  
    for k in range(n):
        if i!=k:
            #Setting forces to function above
            force_x,force_y,force_z = Forces_Function(i,k) 
            #adding incrimental forces 
            force_x_res +=force_x
            force_y_res +=force_y
            force_z_res +=force_z
            
            #returns all resultant force components
        return force_x_res,force_y_res,force_z_res 
        

#get accelerations from resultant force using Newton's 2nd law 
        
def Accelerations_Function():
    #Initally setting all components of acceleration to zero in an n sized array. 
    acc_x=np.zeros(n)
    acc_y=np.zeros(n)
    acc_z=np.zeros(n)
    
    for i in range(n):
        force_x,force_y,force_z=Resultant_Forces_Function(i) 
        #Updating acceleration using Newton's second law
        acc_x[i]=force_x/masses[i]
        acc_y[i]=force_y/masses[i]
        acc_z[i]=force_z/masses[i]
        #Returns acceleration components 
    return acc_x,acc_y,acc_z

#animation as before but over n bodies
def Animation_N_Body(i):
    global vel_x,vel_y,vel_z,x_pos,y_pos,z_pos
    
    #Gets acceleration of bodies from function.
    acc_x,acc_y,acc_z=Accelerations_Function() 
    
    #Update velocity using simple SUVATs.
    vel_x += acc_x*time_int
    vel_y += acc_y*time_int
    vel_z += acc_z*time_int
    
    #Update coords of bodies using simple SUVATs. 
    x_pos += vel_x*time_int
    y_pos += vel_y*time_int
    z_pos += vel_z*time_int
    
    #Setting x,y coords to line object that can be animated. 
    line7.set_xdata(x_pos)
    line7.set_ydata(y_pos)
    
    return line7, #Returns an x,y line to be animated.


# =============================================================================
# Modelling Star Cluster Data
# =============================================================================

G=6.67*10**-11 #Gravitational constant 
time_int=0.5e12 #time interval

n=1000 #Number of bodies/stars in cluster 

M_s=1.989e30 #Sun mass in kg

#Seed so each run can be compared, good seeds: 1,2,5,8,
np.random.seed(2)

#Stars have a random mass between 1/1000 and 1000x the mass of the Sun 
masses=np.random.uniform(low=1/1000 *M_s, high=1000*M_s, size=(n,))

#After formation all stars have approximatley the same velcoities
vel_x= np.random.uniform(low=-25, high=25, size=(n,))
vel_y= np.random.uniform(low=-25, high=25, size=(n,))
vel_z= np.random.uniform(low=-25, high=25, size=(n,))

#Stars start at a random location in x and y between -1ly and +1ly
x_pos= np.random.uniform(low=-9.461e+15, high=9.461e+15, size=(n,))
y_pos= np.random.uniform(low=-9.461e+15, high=9.461e+15, size=(n,))
z_pos= np.random.uniform(low=-9.461e+15, high=9.461e+15, size=(n,))


# =============================================================================
# Modelling Star Cluster Animation
# =============================================================================

#Opening figure window and a set of axes.
fig6=plt.figure(figsize=(14,14))
ax6=plt.axes()
ax6.set_facecolor('black') #Setting plot background colour to black to fit theme.

#Plotting a title that changes N= when n is changed in code.
plt.title("Star Cluster Formation with N=%i Stars"%n,fontsize=24)

line7, = ax6.plot([],[],'r*') #Plotting all coords as red dots.

#Setting axes limits to 4 light years
plt.xlim(-4*9.461e+15,4*9.461e+15)
plt.ylim(-4*9.461e+15,4*9.461e+15)


#Changing axes from m to ly 
xlabels=[item.get_text() for item in ax6.get_xticklabels()]
xlabels[2]="-2ly"
xlabels[4]="0"
xlabels[6]="2ly"
ax6.set_xticklabels(xlabels)
ylabels=[item.get_text() for item in ax6.get_yticklabels()]
ylabels[2]="-2ly"
ylabels[4]="0"
ylabels[6]="2ly"
ax6.set_yticklabels(ylabels)
#Make tick sizes larger so more readable
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

#Plot axes labels.
ax6.set_xlabel("$X$ Distance from (0ly,0ly)",fontsize=18)
ax6.set_ylabel("$Y$ Distance from (0ly,0ly)",fontsize=18)

#Animating and displaying. Using blit True as renders faster for so many stars. 
ani3=animation.FuncAnimation(fig6,Animation_N_Body,interval=10,blit=True)
plt.show()

