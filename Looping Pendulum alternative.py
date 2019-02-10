"""
Looping Pendulum mechanics
This model both big and small masses accounting for Capstan axle friction.  
When T2max exceeds M's weight, M stops, and this program allows m to loop alone. 
"""
import math
import matplotlib.pyplot as plt
pi = math.pi
import numpy

g = 9.81
h = .00001        # small timestep
count = 560000       # counter
# Material constants
m = 0.2           # light mass (looping)
M = 0.5           # heavy mass
mu = 0.03         # axle friction coefficient
k = 0.0001        # air drag coefficient for -kv force model
R = 0.01          # radius of axle
length = 1        # total length of string -> a constant. MEASURE!

#Initial release conditions for state array
theta0 = 0        # starting angle (from top)
theta_dot0 = 5    # starting angular velocity (from rest)
y2 = 0.2          # starting vertical string length to M
y2_dot0 = 0       # starting vertical speed. downwards is negative!
r = length - y2 - R*(pi/2 + theta0)  # geometry, starting string length to m
r_dot0 = 0        # starting radial velocity

# State vector order: theta, theta dot, r, r_dot, y2, y2_dot, y2_acc, T, T2
state0 = numpy.array([theta0, theta_dot0, r, r_dot0, y2_dot0])
print state0
print "Initial state: theta, theta dot, r, r_dot0, y2_dot0"
check_length0 = y2 + R*(state0[0] + pi/2) + state0[2]    
print "Expected length:", length, "Calculated langth:", check_length0

list_theta = []          # initialise list for angle data
list_thetadot = []
time = []
list_y2_dot = []         # check stoppage
list_x = []              # for trajectory visualisation
list_y = []
list_check_length = []   # check total string length

def bothmove(state):
    """Input: state array as defined above containing information about the system
    currently, at step n: [theta, theta_dot, r, r_dot, y2_dot]
    Propagates motion for one time step (h seconds) forward, assuming both masses
    are moving (i.e. M has not stopped yet)    
    Does NOT know whether M stops
    Output: array for next timestep, n + 1: [theta, theta_dot, r, r_dot, y2_dot]"""
    theta1 = state[0]      # current values at n step. read from array
    theta_dot1 = state[1]
    r1 = state[2]
    r_dot1 = state[3]
    y2_dot1 = state[4]
       
    # Step 1: Get next string r distance
    r2 = r1 + r_dot1*h                     
    # Step 2: Solve for angular velocity using tangential acceleration formula    
    theta_dot2 = (g*math.cos(theta1)*h + r1*theta_dot1)/(r2*math.cos(theta_dot1*h))
    # Step 3: Get next angle
    theta2 = theta1 + theta_dot2*h                
    # Step 4: Solve for tension to m by centripetal acceleration
    T12 = m*r2*theta_dot2**2 + m*g*math.sin(theta2)  
    # Step 5: Solve for tension to M by Capstan equation. This is max tension
    T22 = T12*math.exp(mu*(theta2 + pi/2))
    # Step 6: Obtain M acceleration by Newton's law
    y2_acc2 = (T22 - M*g)/M
    # Step 7: Get M vertical speed
    y2_dot2 = y2_dot1 + y2_acc2*h
    # Step 8: Get radial velocity of string
    r_dot2 = -R*theta_dot2 - y2_dot2
       
    state2 = numpy.array([theta2, theta_dot2, r2, r_dot2, y2_dot2])
    return state2
    
def onlymloops(state):
    """Input: state array as defined above containing information about the system
    currently, at step n: [theta, theta dot, r, r_dot, y2_dot]
    Given M has stopped, this propagates motion for one time step (h seconds) 
    forward. Only m will move. *Does* evaluate whether T2 is too small for M. 

    If max friction >= Mg, output y2_dot = 0
    elif max friction < Mg, output y2_dot by Newton Law so as to start bothmove
    Output: state array for next timestep, n + 1: [theta, theta dot, r, r_dot, y2_dot]"""
    
    theta1 = state[0]      # current values at n step. read from array
    theta_dot1 = state[1]
    r1 = state[2]
    r_dot1 = state[3]      # we don't need state[4], y2_dot, to compute here
    
    # Step 1: conservation of string if M stops
    r2 = r1 - R*theta_dot1*h   
    # Step 2: tangential acceleration
    theta_dot2 = (g*math.cos(theta1)*h + r1*theta_dot1)/(r2*math.cos(theta_dot1*h)) 
    # Step 3: next step angle    
    theta2 = theta1 + theta_dot2*h
    # Step 4: centripetal acceleration for T12
    T12 = m*r2*theta_dot2**2 + m*g*math.sin(theta2)    # centripetal acc for T
    # Step 5: string speed conservation    
    r_dot2 = -R*theta_dot2    
    
    # Test whether M will slip, with theoretical max limit for T2 by Capstan equation    
    T22_max = T12*math.exp(mu*(theta2 + pi/2))    
    if T22_max >= M*g:                # M remains stationary, y2_dot2 = 0
        state2 = numpy.array([theta2, theta_dot2, r2, r_dot2, 0])
        return state2
    elif T22_max < M*g:               # M slips!
        y2_acc2 = (T22_max - M*g)/M   # By vertical force balance, assign max friction
        y2_dot2 = y2_acc2*h           # accelerate from rest
        state2 = numpy.array([theta2, theta_dot2, r2, r_dot2, y2_dot2])
        return state2                 # once y2_dot2 <0, bothmove shall run

print "=====Start simulation====="
Mstop = False # a flag for whether M is moving. 
theta = theta0 
theta_dot = theta_dot0
n = 0
state_now = state0
for n in range(count):
    if Mstop == False:                     # M is falling
        state_new = bothmove(state_now)    # propagate both M, m forward by h
    elif Mstop == True:                    # if Mstop is true, M not falling
        state_new = onlymloops(state_now)  # only loop m forward, M stays put
    
    if state_new[4] >= 0:    # if y2_dot (M speed) stalls to 0, switch 
        Mstop = True
      
    x = state_new[2]*math.cos(state_new[0]) + R*math.sin(state_new[0]) 
    y = - state_new[2]*math.sin(state_new[0]) + R*math.cos(state_new[0])         # check sign lol
    y2 = y2 + state_new[4]*h           # to check length of M string
    check_length = y2 + R*(state_new[0] + pi/2) + state_new[2]    
    
    time.append(n*h)                   # save time in sec, to list
    list_theta.append(state_new[0])    # save angle to list
    list_x.append(x)                   # save m trajectory for analysis
    list_y.append(y)        
    list_y2_dot.append(state_new[4])   # append M speed
    list_check_length.append(round(check_length,4))  # string check, don't need so many dp
    state_now = state_new

print "Done! Timestep:", h
#Plot quantities
plt.figure(figsize=(7,7))              # plot theta (angle) across time
plt.plot(time, list_theta, 'b')
plt.ylabel('Theta Position (rad)')
plt.xlabel('Time (sec)')
plt.show() 

plt.figure(figsize=(7,7))              # plot theta (angle) across time
plt.plot(list_x, list_y, 'r')
plt.ylabel('m y position (m)')
plt.xlabel('m x position (m)')
plt.show() 

plt.figure(figsize=(7,12))             # plot M velocity across time
plt.plot(time, list_y2_dot, 'g')
plt.ylabel('M velocity (m/s)')
plt.xlabel('Time (sec)')
plt.show() 

print "Check total string length. See axis label carefully!"
plt.figure(figsize=(7,7))              # plot string length across time
plt.plot(time, list_check_length, 'r')
plt.ylabel('Total string length check (m)')
plt.xlabel('Time (sec)')
plt.show()                             # string length decreases a little actually