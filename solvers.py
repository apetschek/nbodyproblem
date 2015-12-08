import numpy as np

# Runge-Kutta solver (fixed time step)
def rk4(accel,m,r,h,v):
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

# Runge-Kutta-Fehlberg (Adaptive timestep) 
def rk5(accel,m,r,h,v):
    """
    alternate K calculations for position (r) and velocity (v)
    velocity Ks call prior position K's (Kr) and position calls prior velocity K's (Kv)
    this makes sense when you consider the units of Ks, v, and r below 
    """
    k1v = accel(m,r)
    k1r = v
    k2v = accel(m,r + 0.25*k1r*h)
    k2r = v + (0.25*k1v)*h
    k3v = accel(m,r + (3/32.*k1r + 9/32.*k2r)*h)
    k3r = v + (3/32.*k1v + 9/32.*k2v)*h
    k4v = accel(m,r + (1932/2197.*k1r - 7200/2197.*k2r + 7296/2197.*k3r)*h)
    k4r = v + (1932/2197.*k1v - 7200/2197.*k2v + 7296/2197.*k3v)*h
    k5v = accel(m,r + (439/216.*k1r - 8*k2r + 3680/513.*k3r - 845/4104.*k4r)*h)
    k5r = v + (439/216.*k1v - 8*k2v + 3680/513.*k3v - 845/4104.*k4v)*h
    k6v = accel(m,r - (8/27.*k1r + 2*k2r - 3544/2565.*k3r + 1859/4104.*k4r - 11/40.*k5r)*h)
    k6r = v - (8/27.*k1v + 2*k2v - 3544/2565.*k3v + 1859/4104.*k4v - 11/40.*k5v)*h

    # 5th order calculation
    new_v5 = v + h*(16/135.*k1v + 6656/12825.*k3v+28561/56430.*k4v - 9/50.*k5v + 2/55.*k6v) 
    new_r5 = r + h*(16/135.*k1r + 6656/12825.*k3r+28561/56430.*k4r - 9/50.*k5r + 2/55.*k6r) 
    
    return new_v5, new_r5

# Runge-Kutta solver (fixed time step)
def rk8(accel,m,r,h,v):
   """
   alternate K calculations for position (r) and velocity (v)
   velocity Ks call prior position K's (Kr) and position calls prior velocity K's (Kv)
   this makes sense when you consider the units of Ks, v, and r below
   """ 
   k1v = accel(m,r)
   k1r = v
   k2v = accel(m,r + 0.25*k1r*h)
   k2r = v + (0.25*k1v)*h
   k3v = accel(m,r + (5/72.*k1r + 1/72.*k2r)*h)
   k3r = v + (5/72.*k1v + 1/72.*k2v)*h
   k4v = accel(m,r + (1/32.*k1r +3/32.*k3r)*h)
   k4r = v + (1/32.*k1v +3/32.*k3v)*h
   k5v = accel(m,r + (106/125.*k1r- 408/125.*k3r + 352/125.*k4r)*h)
   k5r = v + (106/125.*k1v- 408/125.*k3v + 352/125.*k4v)*h
   k6v = accel(m,r + (1/48.*k1r+ 8/33.*k4r - 125/528.*k5r)*h)
   k6r = v + (1/48.*k1v+ 8/33.*k4v - 125/528.*k5v)*h
   k7v = accel(m,r + (-13893*k1r+ 39936*k4r -64125*k5r+ 60720*k6r)*h/26411.)
   k7r = v +(-13893*k1v+ 39936*k4v -64125*k5v+ 60720*k6v)*h/26411.
   k8v = accel(m,r + (37/392.*k1r+ 1625/9408.*k5r -2/15.*k6r+ 61/6720*k7r)*h)
   k8r = v + (37/392.*k1v+ 1625/9408.*k5v -2/15.*k6v+ 61/6720*k7v)*h
   k9v = accel(m,r +(17176/25515.*k1r - 47104/25515.*k4r + 1325/504.*k5r - 41792/25515.*k6r + 20237/145800.*k7r + 4312/6075.*k8r)*h)
   k9r = v + (17176/25515.*k1v - 47104/25515.*k4v + 1325/504.*k5v - 41792/25515.*k6v + 20237/145800.*k7v + 4312/6075.*k8v)*h
   k10v = accel(m,r + ( -23834/180075.*k1r - 77824/1980825.*k4r- 636635/633864.*k5r + 254048/300125.*k6r - 183/7000.*k7r + 8/11.*k8r - 324/3773.*k9r)*h)
   k10r = v + ( -23834/180075.*k1v - 77824/1980825.*k4v- 636635/633864.*k5v + 254048/300125.*k6v - 183/7000.*k7v + 8/11.*k8v - 324/3773.*k9v)*h
   k11v= accel(m,r + (12733/7600.*k1r - 20032/5225.*k4r + 456485/80256.*k5r - 42599/7125.*k6r + 339227/912000.*k7r - 1029/4108.*k8r + 1701/1408.*k9r + 5145/2432.*k10r)*h)
   k11r = v + (12733/7600.*k1v - 20032/5225.*k4v + 456485/80256.*k5v - 42599/7125.*k6v + 339227/912000.*k7v - 1029/4108.*k8v + 1701/1408.*k9v + 5145/2432.*k10v)*h
   k12v = accel(m,r + h*(-27061/204120.*k1r + 40448/280665.*k4r -1353775/1197504.*k5r + 17662/25515.*k6r - 71687/1166400.*k7r + 98/225.*k8r + 1/16.*k9r + 3773/11664.*k10r))
   k12r = v + h*(-27061/204120.*k1v + 40448/280665.*k4v -1353775/1197504.*k5v + 17662/25515.*k6v - 71687/1166400.*k7v + 98/225.*k8v + 1/16.*k9v + 3773/11664.*k10v)
   k13v = accel(m,r + h*(11203/8680.*k1r - 38144/11935.*k4r + 2354425/458304.*k5r - 84046/16275.*k6r + 673309/1636800.*k7r + 4704/8525.*k8r + 9477/10912.*k9r - 1029/992.*k10r + 19/341.*k12r))
   k13r = v + h*(11203/8680.*k1v - 38144/11935.*k4v + 2354425/458304.*k5v - 84046/16275.*k6v + 673309/1636800.*k7v + 4704/8525.*k8v + 9477/10912.*k9v - 1029/992.*k10v + 19/341.*k12v)


   new_v8 = v + h*(13/288.*k1v +32/125.*k6v + 31213/144000.*k7v + 2401/12375.*k8v + 1701/14080.*k9v + 2401/19200.*k10v + 19/450.*k11v) 
   new_r8 = r + h*(13/288.*k1r +32/125.*k6r + 31213/144000.*k7r + 2401/12375.*k8r + 1701/14080.*k9r + 2401/19200.*k10r + 19/450.*k11r) 
   
   return new_v8,new_r8

# Runge-Kutta-Fehlberg (Adaptive timestep) 
def rk_adaptive(accel,m,r,h,v,recur,emin=10**-12,emax=10**-8,hmax=.1,hmin=.01,recurmax=100):
    """
    alternate K calculations for position (r) and velocity (v)
    velocity Ks call prior position K's (Kr) and position calls prior velocity K's (Kv)
    this makes sense when you consider the units of Ks, v, and r below 
    """
    k1v = accel(m,r)
    k1r = v
    k2v = accel(m,r + 0.25*k1r*h)
    k2r = v + (0.25*k1v)*h
    k3v = accel(m,r + (3/32.*k1r + 9/32.*k2r)*h)
    k3r = v + (3/32.*k1v + 9/32.*k2v)*h
    k4v = accel(m,r + (1932/2197.*k1r - 7200/2197.*k2r + 7296/2197.*k3r)*h)
    k4r = v + (1932/2197.*k1v - 7200/2197.*k2v + 7296/2197.*k3v)*h
    k5v = accel(m,r + (439/216.*k1r - 8*k2r + 3680/513.*k3r - 845/4104.*k4r)*h)
    k5r = v + (439/216.*k1v - 8*k2v + 3680/513.*k3v - 845/4104.*k4v)*h
    k6v = accel(m,r - (8/27.*k1r + 2*k2r - 3544/2565.*k3r + 1859/4104.*k4r - 11/40.*k5r)*h)
    k6r = v - (8/27.*k1v + 2*k2v - 3544/2565.*k3v + 1859/4104.*k4v - 11/40.*k5v)*h

    # 4th order calculation
    new_v4 = v + h*(25/216.*k1v + 1408/2565.*k3v + 2197/4104.*k4v - 1/5.*k5v)
    new_r4 = r + h*(25/216.*k1r + 1408/2565.*k3r + 2197/4104.*k4r - 1/5.*k5r)
    
    # 5th order calculation
    new_v5 = v + h*(16/135.*k1v + 6656/12825.*k3v+28561/56430.*k4v - 9/50.*k5v + 2/55.*k6v) 
    new_r5 = r + h*(16/135.*k1r + 6656/12825.*k3r+28561/56430.*k4r - 9/50.*k5r + 2/55.*k6r) 

    # Calculate truncation error between 5th and 4th order
    eps = np.abs( (np.max(np.abs(new_r5)) - np.max(np.abs(new_r4))) / np.max(np.abs(new_r4)))
    
    # Compare eps to emin and emax and update h accordingly
    if np.max(eps) < emin:
        if h*2.0 < hmax:
            h *= 2.0
        new_v = new_v5
        new_r = new_r5        
    
    if np.max(eps) > emax:
        if h/2.0 > hmin:
            h /= 2.0
            print h
            # Error too large, call rk_adaptive again with smaller h
            if recur < recurmax:
                recur += 1
                rk_adaptive(accel,m,r,h,v,recur)
        new_v = new_v5
        new_r = new_r5
    
    else:
        new_v = new_v5
        new_r = new_r5
    
    return new_v, new_r, h

# Symplectic solvers:

# Stromer Verlet
def leapfrog2(accel,m,r,h,v):
    # 2nd order
    vhalf = v + 0.5*h*accel(m,r)
    new_r = r + h*vhalf
    new_v = vhalf + 0.5*h*accel(m,r)
    return new_v, new_r

# Forest-Ruth leapfrog solver http://arxiv.org/pdf/cond-mat/0110585v1.pdf
def leapfrog4(accel,m,r,h,v):
    # 4th order 
    eps =  .1644986515575760
    lam = -.02094333910398989
    chi =  .1235692651138917
    v += accel(m,r)*eps*h
    r += v*(1-2*lam)*h*0.5
    v += accel(m,r)*chi*h
    r += v*lam*h
    v += accel(m,r)*(1-2*(chi+eps))*h
    r += v*lam*h
    v += accel(m,r)*chi*h
    new_r = r + v*(1-2*lam)*h*0.5
    new_v= v + accel(m,new_r)*eps*h
    return new_v, new_r

def adamsbash(accel, m, r, h, v, step, new_v0, new_r0, new_v1, new_r1, new_v2, new_r2):
    if step == 0:
        new_v0, new_r0 = rungekutta(accel, m, r, h, v)
        return new_v0, new_r0,new_v1, new_r1, new_v2, new_r2
    elif step ==1:
        new_v1, new_r1 = rungekutta(accel, m, r, h, v)
        return new_v0, new_r0, new_v1, new_r1, new_v2, new_r2
    elif step ==2:
        new_v2, new_r2 = rungekutta(accel, m, r, h, v)
        new_v0 = 0; new_r0=0; new_v1=0; new_r1=0
        return new_v0, new_r0, new_v1, new_r1, new_v2, new_r2
    else:
        new_v, new_r = rungekutta(accel, m, r, h, v) + h*(1/12.)*(23*rungekutta(accel, m, new_r))
