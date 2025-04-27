#Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#User Inputs

"All inputs should be in SI units unless otherwise noted"

projectile_mass = 1000
initial_velocity = 5500
launch_angle_degrees = 19.4
launch_angle_radians = (np.pi / 180) * (launch_angle_degrees)

C_d = 0.15
Rocket_CS_Area = np.pi * (1 ** 2) / 4
empty_mass = 420  #Mass of rocket after fuel is spent
yaw_degrees = -18
yaw_radians = yaw_degrees*np.pi/180


Isp = 300   #in seconds
max_accel = 100 * 9.81

#Contstants and Initializations
G= 6.67 * 10**(-11)
M_e = 5.972 * 10**24
R_e = 6.378 * 10**6
g0 = 9.81
v_sin = initial_velocity * np.sin(launch_angle_radians)
v_cos = initial_velocity * np.cos(launch_angle_radians)
is_engine_on = False

#Air density data and curve fitting
altitude_data = (10 ** 3) * np.array([
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    15, 20, 25, 30, 40, 50, 60, 70, 80
])

measured_density_data = np.array([
    1.225, 1.112, 1.007, 0.9093, 0.8194,
    0.7364, 0.6601, 0.59, 0.5258, 0.4671,
    0.4135, 0.1948, 0.08891, 0.04008,
    0.01841, 0.003996, 0.001027, 0.0003097,
    0.0008283, 0.00001846
])

coefficients = np.polyfit(altitude_data, measured_density_data, 6)

#Functions
air_density = np.poly1d(coefficients)

def gravity(h_center, x, y):
    global G, M_e
    g = -G * M_e / (h_center)**2
    gx = g * x / h_center
    gy = g * y / h_center
    return gx, gy

def propulsion(current_atltiude, current_mass, Isp, empty_mass, yaw):
    global max_accel, is_engine_on

    if is_engine_on == False:
        "Activate rocket engine at set altitude"
        activation_condition = current_atltiude >= 300000
    else:
        "Keep engine activated until fuel runs out"
        activation_condition = current_mass >= empty_mass

    if activation_condition:
        is_engine_on = True
        ta = max_accel
        tax = ta*np.cos(yaw)
        tay = ta*np.sin(yaw)
        mdot = -ta * current_mass / (Isp * 9.81)

    else:
        tax = 0
        tay = 0
        mdot = 0

    return tax, tay, mdot

def drag(altitude, xdot, ydot, current_mass):
    global C_d, Rocket_CS_Area
    velocity = np.sqrt(xdot**2 + ydot**2)
    inst_theta = np.arctan(ydot/xdot)
    inst_density = air_density(altitude)
    if inst_density<0:
        inst_density = 0
    d = -0.5 * inst_density * (velocity ** 2) * C_d * Rocket_CS_Area / current_mass
    dx = d * np.cos(inst_theta)
    dy = d * np.sin(inst_theta)
    return dx, dy


def dVdt(S, t):
    global empty_mass
    [xdot,ydot,x,y,m] = S
    h_center = np.sqrt(x**2 + y**2)
    altitude = h_center-R_e

    #Gravitational Accelartion#
    [gx, gy] = gravity(h_center, x, y)

    #Thrust Acceleration#
    [tax, tay, mdot] = propulsion(altitude, m, Isp, empty_mass,yaw_radians)

    #Drag Acceleration#
    [dx, dy] = drag(altitude, xdot, ydot, m)

    xdotdot = gx + tax + dx

    ydotdot = gy + tay + dy

    if altitude <=0:
        [xdotdot, ydotdot, xdot, ydot, mdot] = [0,0,0,0,0]

    return xdotdot, ydotdot, xdot, ydot, mdot

#Main Code

t = np.linspace(0,10000,5000)
S0 = (v_cos+460, v_sin, 0, R_e+0.001, projectile_mass)
Sol = odeint(dVdt, S0, t)

[xdot_rocket, ydot_rocket, x_rocket, y_rocket, mass] = Sol.T

rocket_fire_time_mask = (mass >= empty_mass) & (mass < projectile_mass)
powered_y = y_rocket[rocket_fire_time_mask]
powered_x = x_rocket[rocket_fire_time_mask]


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