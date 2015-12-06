import numpy as np

#RK4 solver
def rungekutta(accel,m,r,h,v):
    """
    alternate K calculations for position (r) and velocity (v)
    velocity Ks call prior position K's (Kr) and position calls prior velocity K's (Kv)
    this makes sense when you consider the units of Ks, v, and r below
    """ 
    k1v = accel(m,r) 
    k1r = v 
    k2v = accel(m,r + h*0.5*k1r) 
    k2r = v+k1v*h*0.5 
    k3v = accel(m,r + h*0.5*k2r) 
    k3r = v+k2v*h*0.5
    k4v = accel(m,r + h*k3r) 
    k4r = v+k3v*h
    new_v = v + h*(k1v + 2*k2v + 2*k3v + k4v)/float(6)
    new_r = r + h*(k1r + 2*k2r + 2*k3r + k4r)/float(6)
    return new_v,new_r

#runge-kutta-fehlberg adaptive solver
def rkf45(accel,m,r,h,v,recur,emin=10**-12,emax=10**-8,hmax=.001,hmin=.000001,recurmax=100):
    """
    alternate K calculations for position (r) and velocity (v)
    velocity Ks call prior position K's (Kr) and position calls prior velocity K's (Kv)
    this makes sense when you consider the units of Ks, v, and r below 
    """
    k1v = accel(m,r)
    k1r = v
    k2v = accel(m,r + 0.25*k1r)
    k2r = v + 0.25*k1v
    k3v = accel(m,r + 3/32.*k1r + 9/32.*k2r)
    k3r = v + 3/32.*k1v + 9/32.*k2v
    k4v = accel(m,r + 1932/2197.*k1r - 7200/2197.*k2r + 7296/2197.*k3r)
    k4r = v + 1932/2197.*k1v - 7200/2197.*k2v + 7296/2197.*k3v
    k5v = accel(m,r + 439/216.*k1r - 8*k2r + 3680/513.*k3r - 845/4104*k4r)
    k5r = v + 439/216.*k1v - 8*k2v + 3680/513.*k3v - 845/4104.*k4v
    k6v = accel(m,r - 8/27.*k1r + 2*k2r - 3544/2565.*k3r + 1859/4104.*k4r - 11/40.*k5r)
    k6r = v - 8/27.*k1v + 2*k2v - 3544/2565.*k3v + 1859/4104.*k4v - 11/40.*k5v

    #4th order calculation
    new_v4 = v + h*(25/216.*k1v + 1408/2565.*k3v + 2197/4104.*k4v - 1/5.*k5v)
    new_r4 = r + h*(25/216.*k1r + 1408/2565.*k3r + 2197/4104.*k4r - 1/5.*k5r)
    
    #5th order calculation
    new_v5 = v + h*(16/135.*k1v + 6656/12825.*k3v+28561/56430.*k4v - 9/50.*k5v + 2/55.*k6v) 
    new_r5 = r + h*(16/135.*k1r + 6656/12825.*k3r+28561/56430.*k4r - 9/50.*k5r + 2/55.*k6r) 

    #calculate truncation % error between 5th and 4th order
    eps = (new_r5 - new_r4)/new_r4

    #compare eps to emin and emax and update h accordingly
    if np.max(eps) < emin:
        if h*2.0 < hmax:
            h *= 2.0
        new_v = new_v5
        new_r = new_r5        
    
    if np.max(eps) > emax:
        if h/2.0 > hmin:
            h /= 2.0
            #error too large, call rkf45 again with smaller h
            if recur < recurmax:
                recur += 1
                rkf45(accel,m,r,h,v,recur)
        new_v = new_v5
        new_r = new_r5
    
    else:
        new_v = new_v5
        new_r = new_r5
    
    return new_v, new_r, h

#leapfrog solver
def leapfrog(accel,m,r,h,v):
    #velocity Verlet (2nd order)
    vhalf = v + 0.5*h*accel(m,r)
    new_r = r + h*vhalf
    new_v = vhalf + 0.5*h*accel(m,r)
    return new_v, new_r

#Forest-Ruth leapfrog solver 
def leapfrogFR(accel,m,r,h,v):
    #Forest-Ruth algorithm (4th order) 
    #http://arxiv.org/pdf/cond-mat/0110585v1.pdf
    eps =  0.1786178958448091 
    lam = -0.2123418310626054
    chi = -0.06626458266981849
    r += v*eps*h
    v += accel(m,r)*(1-2*lam)*h*0.5
    r += v*chi*h
    v += accel(m,r)*lam*h
    r += v*(1-2*(chi+eps))/h
    v += accel(m,r)*lam*h
    r += v*chi*h
    new_v = v + accel(m,r)*(1-2*lam)*h*0.5
    new_r = r + new_v*eps*h
    return new_v, new_r









