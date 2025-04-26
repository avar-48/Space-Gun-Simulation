#Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#User Inputs

"All inputs should be in SI units unless otherwise noted"

projectile_mass = 1000
initial_velocity = 5500
launch_angle_degrees = 17.75
launch_angle_radians = (np.pi / 180) * (launch_angle_degrees)

C_d = 0.15
A = np.pi * (0.15**2) / 4
empty_mass = 457


Isp = 300   #in seconds
max_accel = 100 * 9.81

#Contstants and Initializations
G= 6.67 * 10**(-11)
M_e = 5.972 * 10**24
R_e = 6.378 * 10**6
g0 = 9.81
v_sin = initial_velocity * np.sin(launch_angle_radians)
v_cos = initial_velocity * np.cos(launch_angle_radians)
engine_on_time = 0
engine_off_time = 0
is_engine_on = False

#Functions

def gravity(h_center, x, y):
    g = -G*M_e/(h_center)**2
    gx = g * x / h_center
    gy = g * y / h_center
    return gx, gy

def propulsion(current_atltiude, current_mass, Isp, empty_mass, t):
    global engine_on_time,engine_off_time, max_accel, is_engine_on

    if is_engine_on == False:
        "Activate rocket engine at set altitude"
        activation_condition = current_atltiude >= 300000
    else:
        "Keep engine activated until fuel runs out"
        activation_condition = current_mass >= empty_mass

    if activation_condition:
        is_engine_on = True
        ta = max_accel
        mdot = -ta*current_mass/(Isp*9.81)
        if engine_on_time==0:
            engine_on_time = t
        else:
            engine_off_time = t
    else:
        ta = 0
        mdot = 0
    return ta, mdot

def dVdt(S, t):
    [xdot,ydot,x,y,m] = S
    h_center = np.sqrt(x**2 + y**2)
    altitude = h_center-R_e

    #Gravitational Accelartion#
    [gx, gy] = gravity(h_center, x, y)

    #Thrust Acceleration#
    [ta, mdot] = propulsion(altitude, m, Isp, 457,t)

    yaw = -20*np.pi/180#np.arctan(ydot/xdot)

    tax = ta*np.cos(yaw)
    tay = ta*np.sin(yaw)

    #Drag Acceleration#
    d = 0

    xdotdot = gx + tax

    ydotdot = gy + tay

    if altitude <=0:
        [xdotdot, ydotdot, xdot, ydot, mdot] = [0,0,0,0,0]

    return xdotdot, ydotdot, xdot, ydot, mdot

#Main Code

t = np.linspace(0,10000,5000)
S0 = (v_cos+460, v_sin, 0, R_e+0.001, projectile_mass)
Sol = odeint(dVdt, S0, t, hmax=1)

[xdot_rocket, ydot_rocket, x_rocket, y_rocket, mass] = Sol.T

rocket_fire_time = (t >= engine_on_time) & (t <= engine_off_time)
powered_y = y_rocket[rocket_fire_time]
powered_x = x_rocket[rocket_fire_time]


#Plots
plot1 = plt.subplot2grid((3, 3), (0, 1), rowspan = 3, colspan=3)
plot2 = plt.subplot2grid((3, 3), (0, 0), rowspan = 1, colspan=1)
plot3 = plt.subplot2grid((3, 3), (1, 0), rowspan = 1, colspan=1)
plot4 = plt.subplot2grid((3, 3), (2, 0), rowspan = 1, colspan=1)

#Rocket Trajectory
theta = np.linspace(0,2*np.pi, 1000)
x_earth = R_e*np.cos(theta)
y_earth = R_e*np.sin(theta)
plot1.plot(x_earth, y_earth, 'b')
plot1.plot(x_rocket, y_rocket, 'y')
plot1.plot(powered_x, powered_y, 'r')
plot1.set_title("Trajectory")
plot1.legend(['Earth', 'Unpowered Flight', 'Powered Flight'])
plot1.axes.get_xaxis().set_ticks([])
plot1.axes.get_yaxis().set_ticks([])

#Mass w.r.t time
plot2.plot(t, mass)
plot2.set (ylabel = "Mass (Kg)")

#Altitude
altitude = np.sqrt(x_rocket**2+y_rocket**2)-R_e
plot3.plot(t, altitude)
plot3.set(ylabel = "Altitude (m)")

#Velocity
vel = np.sqrt(xdot_rocket**2 + ydot_rocket**2)
plot4.plot(t, vel)
plot4.set(xlabel="Time (s)", ylabel = "Orbital Velocity (m/s)")

plt.show()
