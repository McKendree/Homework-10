import numpy as np
import matplotlib.pyplot as plt
import time

G = 6.6738e-11 #Newton's constant
M = 1.9891e30 #mass of sun in kg
AU = 1.496e11 #1 astronomical unit in meters

def f(r):
    """derivative function to pass to rk4
    pass state vector r = [x, y, vx, vy]"""

    x, y, vx, vy = r
    rcubed = np.sqrt(x**2 + y**2)**3

    fx = vx #return 1st parameter
    fy = vy #return 2nd parameter
    fvx = -G*M*x/rcubed #return 3rd parameter
    fvy = -G*M*y/rcubed #return last parameter
    return np.array([fx, fy, fvx, fvy])

def rk4_step(r=None, h=None, f=None):
    '''r = [x, y, vx, vy]
    h = step size
    f = derivative function'''
    k1 = h*f(r)
    k2 = h*f(r+0.5*k1)
    k3 = h*f(r+0.5*k2)
    k4 = h*f(r+k3)
    return (1/6)*(k1+2*k2+2*k3+k4)

def run_rk4_fixed(initial_state=None, initial_h=None, t0=None, tmax=None, f=None):
    '''initial_state = array of [x, y, vx, vy]
    initial_h = initial step size
    t0 = initial time
    tmax = final time
    f = derivative function'''
    r = initial_state
    h = initial_h
    t = t0
    xpoints = []
    ypoints = []
    while t<tmax:
        t = t+h
        r = r + rk4_step(r=r, h=h, f=f)
        xpoints.append(r[0])
        ypoints.append(r[1])
    return np.array([xpoints, ypoints])

def run_rk4_adaptive(initial_state=None, initial_h=None, t0=None, tmax=None, f=None):
    '''initial_state = array of [x, y, vx, vy]
    initial_h = initial step size
    t0 = initial time
    tmax = final time
    f = derivative function'''
    r = initial_state
    h = initial_h
    t = t0
    xpoints = []
    ypoints = []
    while t<tmax:
        #Do one large step
        r1 = r + rk4_step(r=r, h=2*h, f=f)
        #Do two small steps
        r2 = r + rk4_step(r=r, h=h, f=f)
        r2 = r2 + rk4_step(r=r2,h=h, f=f)

        #calculate value of rho
        ex = (1/30)*(r1[0]-r2[0])
        ey = (1/30)*(r1[1]-r2[1])
        rho = (h*delta)/(ex**2+ey**2)**0.5

        #calculate new values of t, h, r
        #update points if appropriate
        if rho>=1.0:
            h = h*rho**(1/4)
            t = t+h
            r = r + rk4_step(r=r, h=h, f=f)
            xpoints.append(r[0])
            ypoints.append(r[1])
        else:
            h = h*rho**(1/4)
    return np.array([xpoints, ypoints])

if __name__ == "__main__":
    h0 = 1.0e4 #initial step size
    tmax = 2.0e9 #total time
    t0 = 0 #initial time

    delta = 1e3/(365.25*24*3600) #meters accuracy per second

    x0, y0 = 4e12, 0 #starting pos, 4 billion kilometers
    vx0, vy0 = 0, 500 #starting velocity, m/s
    r0 = np.array([x0, y0, vx0, vy0])

    #runs the fixed rk4 method and records time taken to run
    startTime = time.time()
    xpos, ypos = run_rk4_fixed(initial_state=r0,initial_h=h0,t0=t0,tmax=tmax,f=f)
    endTime = time.time()
    print('Time to run fixed rk4 method:', endTime-startTime, 'seconds')
    print('Step size required for accurate orbit is approx.: <=1e4 meters')
    
    #plots fixed rk4
    plt.plot(xpos/AU, ypos/AU)
    plt.title('Fixed RK4 Method Orbit')
    plt.xlabel('AU')
    plt.ylabel('AU')
    plt.savefig('fixed_rk4_orbit.png')
    plt.show()
    plt.clf()

    h0 = 1.0e5 #initial step size
    tmax = 1.54e9 #total time
    t0 = 0 #initial time

    delta = 1e3/(365.25*24*3600) #meters accuracy per second

    x0, y0 = 4e12, 0 #starting pos, 4 billion kilometers
    vx0, vy0 = 0, 500 #starting velocity, m/s
    r0 = np.array([x0, y0, vx0, vy0])

    #runs the adaptive rk4 method and records time taken to run
    startTime = time.time()
    xpos, ypos = run_rk4_adaptive(initial_state=r0,initial_h=h0,t0=t0,tmax=tmax,f=f)
    endTime = time.time()
    print('\nTime to run adaptive rk4 method:', endTime-startTime, 'seconds')
    print('Initial step size is 1e5 meters')

    #plots adaptive rk4
    plt.plot(xpos/AU, ypos/AU, alpha = 0.5)
    plt.plot(xpos/AU, ypos/AU, 'k.')
    plt.title('Adaptive RK4 Method Orbit')
    plt.xlabel('AU')
    plt.ylabel('AU')
    plt.savefig('adaptive_rk4_orbit.png')
    plt.show()
    plt.clf()
