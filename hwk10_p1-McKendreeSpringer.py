import numpy as np
import matplotlib.pyplot as plt

def f(x,v):
    return v**2-x-5

def leapfrog_method(initial_x=None,
                    initial_v=None,
                    h=None,
                    t0=None,
                    tmax=None,
                    f=None):
    '''initial_x = initial x value
    initial_v = initial v value
    h = time step
    t0 = initial time
    tmax = maximum time
    f = function to use leapfrog method on'''
    x = initial_x
    v = initial_v
    t = t0
    xVals = []
    while t<tmax:
        t = t+h
        x1 = x+v*h+0.5*f(x,v)*h**2
        v1 = v+0.5*(f(x,v)+f(x1,v))*h
        x = x1
        v = v1
        xVals.append(x)
    return np.array(xVals)

if __name__ == "__main__":
    t0 = 0 #initial time
    tmax = 50 #final time
    h = 0.001 #timestep
    x0 = 1 #initial x
    v0 = 0 #initial v

    #runs leapfrog function
    xVals = leapfrog_method(initial_x=x0,initial_v=v0,h=h,t0=t0,tmax=tmax,f=f)

    #plots leapfrog function results
    plt.plot(np.arange(t0,tmax+h,h), xVals, '.')
    plt.title('Leapfrog Method')
    plt.savefig('Leapfrog_Method.png')
    plt.show()
