from NewHydrogenic import *
import FancyIntegrations as fi
import numpy as np
import scipy.integrate as integ
import scipy.special as scp
import sympy.physics.wigner as wg
import scipy.constants as cte
from scipy.optimize import curve_fit as cf


def cg(l1,l2,l3,m1,m2,m3):
    return float(wg.clebsch_gordan(l1,l2,l3,m1,m2,m3).n(10))
#A linear model to extrapolate
def lin(x,a,b):
    return a*x+b

#Radial dipole integral between Regular function at energy EcontAU and angular momentum l, and the hydrogenic orbital with principal quantum number n and orbital angular momentum ln.
def lcont_RadIntRombRegHyd(EcontAU,l,n,ln):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*n**2.)//dBW),int((200)//dBW))+2
    Nromb = 2**(int(np.log2(100*NdBW)))+1
    rf = max(200,5*n**2)
    mesh = np.linspace(1e-5,rf,Nromb,endpoint=True)
    y =  [Ucl(EcontAU,l,xx)*xx*Unl(n,ln,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint

#Radial dipole integral between irregular function at energy EcontAU and angular momentum l, and the hydrogenic orbital with principal quantum number n and orbital angular momentum ln.
def lcont_RadIntRombIregHyd(EcontAU,l,n,ln):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    print(dBW)
    NdBW = max(int((5.*n**2.)//dBW),int((200)//dBW))+2
    Nromb = 2**(int(np.log2(100*NdBW)))+1
    rf = (200,5*n**2)
    mesh = np.linspace(1e-5,rf,Nromb,endpoint=True)
    y =  [Ccl(EcontAU,l,xx)*xx*Unl(n,ln,xx)*SB(max(l,ln),xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint

#Int1 and int2 define the two integrals present in the integral equation of the radial dalgarno Lewis state

#Integrands for positive energy for direct integration
def int1func(r,e,l,n,ln):
    return Ucl(e,l,r)*r*Unl(n,ln,r)

def int2func(r,e,l,n,ln):
    return Ccl(e,l,r)*r*Unl(n,ln,r)

#Functions that recast the integrals as a first order differential equation
def int1funcode(y,r,e,l,n,ln):
    return Ucl(e,l,r)*r*Unl(n,ln,r)

def int2funcode(y,r,e,l,n,ln):
    return -Ccl(e,l,r)*r*Unl(n,ln,r)

#Integrand for negative energies.
def int1funcn(x,e,l,n,ln,rl):
    nu = nuf(e)
    if(x<rl):
        return Ucln(e,l,x)*x*Unl(n,ln,x)
    else:
        return LRUcln(e,l,x)*x*Unl(n,ln,x)

def int2funcn(x,e,l,n,ln,rl):
    nu = np.sqrt(1/(-2*e))
    if(x<rl):
        return W(nu,l,x)*x*Unl(n,ln,x)
    else:
        return LRW(nu,l,x)*x*Unl(n,ln,x)

#Integrand for negative energies using the analytical Green's function.
def int1funcn_ana(x,e,l,n,ln,rl):
    nu = nuf(e)
    if(x<rl):
        return Uclno(e,l,x)*x*Unl(n,ln,x)
    else:
        return LRUclno(e,l,x)*x*Unl(n,ln,x)

def int2funcn_ana(x,e,l,n,ln,rl):
    nu = np.sqrt(1/(-2*e))
    if(x<rl):
        return W_scal(nu,l,x)*x*Unl(n,ln,x)
    else:
        return LRW_scal(nu,l,x)*x*Unl(n,ln,x)


#define the two green's functions for positive and negative energies, and the analytic Green's function.

def G_pos(E, l, r, rp):
    rl = min(rp,r)
    rg = max(rp,r)
    return -np.sqrt(2)*np.pi*Ucl(E, l, rl)*(1/np.sqrt(2) * (-Ccl(E,l,rg)+1j*Ucl(E,l,rg)))

def G_neg(E,l,r,rp):
    nu = nuf(E)
    bet = np.pi*(nuf(E)-l)
    rl =  min(rp,r)
    rg = max(rp,r)
    return -np.pi/np.sin(bet) * Ucln(E,l,rl) * W(nu,l,rg)

def G_ana(E,l,r,rp):
    nu = nuf(E)
    bet = np.pi*(nuf(E)-l)
    rl =  min(rp,r)
    rg = max(rp,r)
    w = 1
    if((bet/np.pi)%1 == 0):
        if(bet>0):
            w = 1e-18
        else:
            w = -2*nu**l / scp.factorial(l-nu)
    else:
        w = -np.sin(bet) * nu**l * scp.gamma(nu-l)*2./np.pi

    return 2/w * Uclno(E,l,rl) * W_scal(nu,l,rg)

#Definition of the radial functions of the Dalgarno-Lewis state. 
#This functions return the mesh where they are calculated and the function evaluated in this mesh.
#It it useful to test if the function fullfils the differential equation.

#Negative energy
def DLradfuncNE(n,ln,k,ldl):
    En = -0.5/n**2
    rf = max(5*n**2,100)
    rmesh = np.linspace(1e-7,rf,2000)
    rl = min(300,50/kappa(En-k))
    dmesh = rmesh[1]-rmesh[0]
    #Integrate the regular function integral
    Integrand = [ int1funcn(rr,En-k,ldl,n,ln,rl) for rr in rmesh]
    I12t = fi.cumulative_composite(Integrand,dmesh,order=7)
    #Multiply by the whittaker function, we do this later to avoid problems with the divergence of W near zero.
    I12 = [W(nuf(En-k),ldl,rmesh[i+1])*I12t[i] if(rmesh[i+1]<rl) else LRW(nuf(En-k),ldl,rmesh[i+1])*I12t[i] for i in range(len(I12t))]
    I12=np.insert(I12,0,0)
    #Integrate the whittaker function integral
    #in this case, since we neeg to ensure that the vanishes at large r, we use the cumulative integration method to chnage variables:
    # We want to get \int_r^{RL} f(r')dr' we do it by computing changing variables to r'=RL-x so the integral goes to -\int_{Rl-r}^0 f(Rl-x)dx.
    # Then we sample the function as f(RL-ri) integrate from ri=0 to Rl, we then integrate this cumulative giving the values of the desrired integral in inverted order.
    irtest = np.array(list(reversed(rmesh)))
    Integrand = [int2funcn(rr,En-k,ldl,n,ln,rl) for rr in irtest]
    I22t = fi.cumulative_composite(Integrand,dmesh,order=7)
    #Now multiply by the regular function and invert to get the right order.
    I22 = np.array(list(reversed([Ucln(En-k,ldl,irtest[i+1])*(I22t[i]) if(irtest[i+1]<rl) else LRUcln(En-k,ldl,irtest[i+1])*(I22t[i]) for i in range(len(I22t))])))
    I22=np.insert(I22,-1,0)
    bet = np.pi*(nuf(En-k)-ldl)
    w = -2/np.pi * np.sin(bet)
    Dlr = 2.0/w*I12+2.0/w*I22
    return [rmesh,Dlr]

#Negative energy with analytic Green's function
def DLradfuncNE_ana(n,ln,k,ldl):
    En = -0.5/(n**2)
    rf = max(5*n**2,100)
    rmesh = np.linspace(1e-7,rf,2000)
    rl = min(300,50/kappa(En-k))
    dmesh = rmesh[1]-rmesh[0]
    #Integrate the regular function integral
    Integrand = [ int1funcn_ana(rr,En-k,ldl,n,ln,rl) for rr in rmesh]
    I12t = fi.cumulative_composite(Integrand,dmesh,order=7)
    #Multiply by the whittaker function, we do this later to avoid problems with the divergence of W near zero.
    I12 = [W_scal(nuf(En-k),ldl,rmesh[i+1])*I12t[i] if(rmesh[i+1]<rl) else LRW_scal(nuf(En-k),ldl,rmesh[i+1])*I12t[i] for i in range(len(I12t))]
    I12=np.insert(I12,0,0)
    #Integrate the whittaker function integral
    #in this case, since we neeg to ensure that the vanishes at large r, we use the cumulative integration method to chnage variables:
    # We want to get \int_r^{RL} f(r')dr' we do it by computing changing variables to r'=RL-x so the integral goes to -\int_{Rl-r}^0 f(Rl-x)dx.
    # Then we sample the function as f(RL-ri) integrate from ri=0 to Rl, we then integrate this cumulative giving the values of the desrired integral in inverted order.
    irtest = np.array(list(reversed(rmesh)))
    Integrand = [int2funcn_ana(rr,En-k,ldl,n,ln,rl) for rr in irtest]
    I22t = fi.cumulative_composite(Integrand,dmesh,order=7)
    #Now multiply by the regular function and invert to get the right order.
    I22 = np.array(list(reversed([Uclno(En-k,ldl,irtest[i+1])*(I22t[i]) if(irtest[i+1]<rl) else LRUclno(En-k,ldl,irtest[i+1])*(I22t[i]) for i in range(len(I22t))])))
    I22=np.insert(I22,-1,0)
    nu = nuf(En-k)
    bet = np.pi*(nuf(En-k)-ldl)
    w = 1
    if((bet/np.pi)%1 == 0):
        if(bet>0):
            print(bet)
            w = 1e-18
        else:
            print("Beta is equal or smaller than zero")
            w = -2*nu**ldl / scp.factorial(ldl-nu)
    else:
        w = -1* (nu**(ldl)) * scp.gamma(nu-ldl) * np.sin(bet) * 2./np.pi
    Dlr = 2.0/w*I12+2.0/w*I22
    return [rmesh,Dlr] 

#Positive energy solved by direct integration (complex valued function)
def DLradfuncPE(n,ln,k,ldl):
    #Generate the function the integrals in the DL in a mesh.
    En = -0.5/n**2
    rf = max(200,5*n**2)
    rmesh = np.linspace(1e-7,rf,2000)
    dmesh = rmesh[1]-rmesh[0]
    #Integrate the regular function integral
    Integrand = [ int1func(xx,En+k,ldl,n,ln) for xx in rmesh]
    I12t = fi.cumulative_composite(Integrand,dmesh,order=7)
    I12t = np.insert(I12t,0,0)
    I12 = np.array([Ccl(En+k,ldl,rmesh[xx])*I12t[xx] for xx in range(len(I12t))])
    #Initial value of for the I2
    irtest = np.array(list(reversed(rmesh)))
    Integrand = [int2func(rr,En+k,ldl,n,ln) for rr in irtest]
    I22t = fi.cumulative_composite(Integrand,dmesh,order=7)
    I22t = np.insert(I22t,0,0)
    #Now multiply by the regular function and invert to get the right order.
    I22 = np.array(list(reversed([Ucl(En+k,ldl,irtest[i])*I22t[i]  for i in range(len(I22t))])))
    #The value of I1 at infinity
    I1inf = lcont_RadIntRombRegHyd(En+k,ldl,n,ln)
    #Use simpsons algorithm to do the radial integral
    yint2 = [(-np.pi*1j*I1inf*Ucl(En+k,ldl,rmesh[xx])+np.pi*I12[xx]+np.pi*I22[xx]) for xx in range(len(rmesh))]
    return [rmesh,yint2]

#Positive energy solved by integration of first order ODE (less accurate)
def DLradfuncPEodeint(n,ln,k,ldl):
    #Generate the function the integrals in the DL in a mesh.
    rf = max(5*n**2,100)
    rmesh = np.linspace(1e-7,rf,2000)
    dmesh = rmesh[1]-rmesh[0]
    #Integrate the regular function integral
    I14 = integ.odeint(int1func,0,rmesh,(En+k,ldl,n,ln)).flatten()
    #Initial value of for the I2
    I2o = lcont_RadIntRombIregHyd(En+k,ldl,n,ln)
    I24 = integ.odeint(int2func,I2o,rmesh,(En+k,ldl,n,ln)).flatten()
    #The value of I1 at infinity
    I1inf = lcont_RadIntRombRegHyd(En+k,ldl,n,ln)
    #Use simpsons algorithm to do the radial integral
    Dlr = [-np.pi*1j*I1inf*Ucl(En+k,ldl,rmesh[xx])+np.pi*Ccl(En+k,ldl,rmesh[xx])*I14[xx]+np.pi*Ucl(En+k,ldl,rmesh[xx])*I24[xx] for xx in range(len(rmesh))]
    return [rmesh,Dlr]


#Radial integral, of the dipole operator, between: 
# --- orbital with principal number f and angular momentum lf
# --- Dalgarno-Lewis radial function with angular momentum ldl at energy E=-0.5/n**2 + k. 
# --- This DL radial function solves the inhomogenous differential equation for the orbital with principal number n and angular momentum ln

#Positive Energy (Solved with direct integration)
def DLmatElemRadintsPEl(f,lf,n,ln,k,ldl):
    nl = max(f,n)
    En = -0.5/n**2
    rf = max(200,5*nl**2)
    rmesh = np.linspace(1e-7,rf,2000)
    dmesh = rmesh[1]-rmesh[0]
    #Integrate the regular function integral
    Integrand = [ int1func(xx,En+k,ldl,n,ln) for xx in rmesh]
    I12t = fi.cumulative_composite(Integrand,dmesh,order=7)
    I12t = np.insert(I12t,0,0)
    I12 = np.array([Ccl(En+k,ldl,rmesh[xx])*I12t[xx] for xx in range(len(I12t))])
    #Initial value of for the I2
    irtest = np.array(list(reversed(rmesh)))
    Integrand = [int2func(rr,En+k,ldl,n,ln) for rr in irtest]
    I22t = fi.cumulative_composite(Integrand,dmesh,order=7)
    I22t = np.insert(I22t,0,0)
    #Now multiply by the regular function and invert to get the right order.
    I22 = np.array(list(reversed([Ucl(En+k,ldl,irtest[i])*I22t[i]  for i in range(len(I22t))])))
    #The value of I1 at infinity
    I1inf = lcont_RadIntRombRegHyd(En+k,ldl,n,ln)
    #Use simpsons algorithm to do the radial integral
    yint2 = [Unl(f,lf,rmesh[xx])*rmesh[xx]*(-np.pi*1j*I1inf*Ucl(En+k,ldl,rmesh[xx])+np.pi*I12[xx]+np.pi*I22[xx]) for xx in range(len(rmesh))]
    #tot = np.array([(-np.pi*1j*I1inf*Ucl(En+k,ldl,rmesh[xx])+np.pi*I12[xx]+np.pi*I22[xx]) for xx in range(len(rmesh))])
    radint2 = fi.simpson(yint2,dx=dmesh)
    return radint2

#Negative energy
def DLmatElemRadintsNEl(f,lf,n,ln,k,ldl):
    nl = max(f,n)
    En = -0.5/n**2
    rf = max(5*nl**2,200)
    rl = min(300,50/kappa(En-k))
    rmesh = np.linspace(1e-7,rf,2000)
    dmesh = rmesh[1]-rmesh[0]
    #Integrate the regular function integral
    Integrand = [ int1funcn(rr,En-k,ldl,n,ln,rl) for rr in rmesh]
    I12t = fi.cumulative_composite(Integrand,dmesh,order=7)
    #Multiply by the whittaker function, we do this later to avoid problems with the divergence of W near zero.
    I12 = [W(nuf(En-k),ldl,rmesh[i+1])*I12t[i] if(rmesh[i+1]<rl) else LRW(nuf(En-k),ldl,rmesh[i+1])*I12t[i] for i in range(len(I12t))]
    I12=np.insert(I12,0,0)
    #Integrate the whittaker function integral
    #in this case, since we neeg to ensure that the vanishes at large r, we use the cumulative integration method to chnage variables:
    # We want to get \int_r^{RL} f(r')dr' we do it by computing changing variables to r'=RL-x so the integral goes to -\int_{Rl-r}^0 f(Rl-x)dx.
    # Then we sample the function as f(RL-ri) integrate from ri=0 to Rl, we then integrate this cumulative giving the values of the desrired integral in inverted order.
    irtest = np.array(list(reversed(rmesh)))
    Integrand = [int2funcn(rr,En-k,ldl,n,ln,rl) for rr in irtest]
    I22t = fi.cumulative_composite(Integrand,dmesh,order=7)
    #Now multiply by the regular function and invert to get the right order.
    I22 = np.array(list(reversed([Ucln(En-k,ldl,irtest[i+1])*(I22t[i]) if(irtest[i+1]<rl) else LRUcln(En-k,ldl,irtest[i+1])*(I22t[i]) for i in range(len(I22t))])))
    I22=np.insert(I22,-1,0)
    bet = np.pi*(nuf(En-k)-ldl)
    w = -2/np.pi * np.sin(bet)
    Dlr = 2.0/w*I12+2.0/w*I22
    #The overlap with the final orbital
    Integrand = Unl(f,lf,rmesh)*rmesh*Dlr
    rint = fi.simpson(Integrand, dx=dmesh)
    return  rint

#Negative energy using the analytic Green's function.
def DLmatElemRadintsNEl_ana(f,lf,n,ln,k,ldl):
    nl = max(f,n)
    En = -0.5/(n**2)
    rf = max(5*nl**2,100)
    rmesh = np.linspace(1e-7,rf,2000)
    rl = min(300,50/kappa(En-k))
    dmesh = rmesh[1]-rmesh[0]
    #Integrate the regular function integral
    Integrand = [ int1funcn_ana(rr,En-k,ldl,n,ln,rl) for rr in rmesh]
    I12t = fi.cumulative_composite(Integrand,dmesh,order=7)
    #Multiply by the whittaker function, we do this later to avoid problems with the divergence of W near zero.
    I12 = [W_scal(nuf(En-k),ldl,rmesh[i+1])*I12t[i] if(rmesh[i+1]<rl) else LRW_scal(nuf(En-k),ldl,rmesh[i+1])*I12t[i] for i in range(len(I12t))]
    I12=np.insert(I12,0,0)
    #Integrate the whittaker function integral
    #in this case, since we neeg to ensure that the vanishes at large r, we use the cumulative integration method to chnage variables:
    # We want to get \int_r^{RL} f(r')dr' we do it by computing changing variables to r'=RL-x so the integral goes to -\int_{Rl-r}^0 f(Rl-x)dx.
    # Then we sample the function as f(RL-ri) integrate from ri=0 to Rl, we then integrate this cumulative giving the values of the desrired integral in inverted order.
    irtest = np.array(list(reversed(rmesh)))
    Integrand = [int2funcn_ana(rr,En-k,ldl,n,ln,rl) for rr in irtest]
    I22t = fi.cumulative_composite(Integrand,dmesh,order=7)
    #Now multiply by the regular function and invert to get the right order.
    I22 = np.array(list(reversed([Uclno(En-k,ldl,irtest[i+1])*(I22t[i]) if(irtest[i+1]<rl) else LRUclno(En-k,ldl,irtest[i+1])*(I22t[i]) for i in range(len(I22t))])))
    I22=np.insert(I22,-1,0)
    nu = nuf(En-k)
    bet = np.pi*(nuf(En-k)-ldl)
    w = 1
    if((bet/np.pi)%1 == 0):
        if(bet>0):
            w = 1e-18
        else:
            w = -2*nu**ldl / scp.factorial(ldl-nu)
    else:
        w = -1* (nu**(ldl)) * scp.gamma(nu-ldl) * np.sin(bet) * 2./np.pi
    Dlr = 2.0/w*I12+2.0/w*I22
    yint = Unl(f,lf,rmesh)*rmesh*Dlr
    rint = fi.simpson(yint, dx=dmesh)
    return rint


#Left Hand Side of the inhomogeneous radial diff equation: (E-H)Dl. Needs the Dalgarno-Lewis radial function sampled in the mesh rt, the angular momentum l, the energy E and the photon energy w.
def lhs(Dlh,rt,l,en,w):
    return (en-w)*Dlh-(-0.5*np.gradient(np.gradient(Dlh,rt),rt)+0.5*l*(l+1)/rt**2 * Dlh -1/rt * Dlh)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Here we define the functions for the matrix element in the Rayleigh-Raman cross section.

#start with the chi_a part of the matrix element. (see paper or notes for this)

#Start with the element that depend on the magnetic numbers and in the initial polarization

#This returns a vector with each spherical components of the position operator.

#For one extra unit of angular momentum
def Xa_ma_mb_e_lp(la,ma,lb,mb,eo):
    ep = -1.0/np.sqrt(2) * (eo[0]+1j*eo[1])
    e0 = eo[2]
    em = 1.0/np.sqrt(2) * (eo[0]-1j*eo[1])

    prefac_lp = np.sqrt((2*la+1)/(2*lb+1))*cg(1,la,la+1,0,0,0)*cg(1,la+1,lb,0,0,0)
    rp_lp = -ep*cg(1,la,la+1,1,ma,ma+1)*cg(1,la+1,lb,1,ma+1,mb)-em*cg(1,la,la+1,-1,ma,ma-1)*cg(1,la+1,lb,1,ma-1,mb)+e0*cg(1,la,la+1,0,ma,ma)*cg(1,la+1,lb,1,ma,mb)
    rm_lp = -ep*cg(1,la,la+1,1,ma,ma+1)*cg(1,la+1,lb,-1,ma+1,mb)-em*cg(1,la,la+1,-1,ma,ma-1)*cg(1,la+1,lb,-1,ma-1,mb)+e0*cg(1,la,la+1,0,ma,ma)*cg(1,la+1,lb,-1,ma,mb)
    r0_lp = -ep*cg(1,la,la+1,1,ma,ma+1)*cg(1,la+1,lb,0,ma+1,mb)-em*cg(1,la,la+1,-1,ma,ma-1)*cg(1,la+1,lb,0,ma-1,mb)+e0*cg(1,la,la+1,0,ma,ma)*cg(1,la+1,lb,0,ma,mb)
    
    return np.array([prefac_lp*rp_lp,prefac_lp*rm_lp,prefac_lp*r0_lp])

#One less
def Xa_ma_mb_e_lm(la,ma,lb,mb,eo):
    ep = -1.0/np.sqrt(2) * (eo[0]+1j*eo[1])
    e0 = eo[2]
    em = 1.0/np.sqrt(2) * (eo[0]-1j*eo[1])

    prefac_lm = np.sqrt((2*la+1)/(2*lb+1))*cg(1,la,la-1,0,0,0)*cg(1,la-1,lb,0,0,0)
    rp_lm = -ep*cg(1,la,la-1,1,ma,ma+1)*cg(1,la-1,lb,1,ma+1,mb)-em*cg(1,la,la-1,-1,ma,ma-1)*cg(1,la-1,lb,1,ma-1,mb)+e0*cg(1,la,la-1,0,ma,ma)*cg(1,la-1,lb,1,ma,mb)
    rm_lm = -ep*cg(1,la,la-1,1,ma,ma+1)*cg(1,la-1,lb,-1,ma+1,mb)-em*cg(1,la,la-1,-1,ma,ma-1)*cg(1,la-1,lb,-1,ma-1,mb)+e0*cg(1,la,la-1,0,ma,ma)*cg(1,la-1,lb,-1,ma,mb)
    r0_lm = -ep*cg(1,la,la-1,1,ma,ma+1)*cg(1,la-1,lb,0,ma+1,mb)-em*cg(1,la,la-1,-1,ma,ma-1)*cg(1,la-1,lb,0,ma-1,mb)+e0*cg(1,la,la-1,0,ma,ma)*cg(1,la-1,lb,0,ma,mb)

    return np.array([prefac_lm*rp_lm,prefac_lm*rm_lm,prefac_lm*r0_lm])


#Follow with the chi_b function

def Xb_ma_mb_e_lp(la,ma,lb,mb,eo):
    ep = -1.0/np.sqrt(2) * (eo[0]+1j*eo[1])
    e0 = eo[2]
    em = 1.0/np.sqrt(2) * (eo[0]-1j*eo[1])

    prefac_lp = np.sqrt((2*lb+1)*(2*la+1))/(2*(lb+1)+1) * cg(1,lb,lb+1,0,0,0)*cg(1,la,lb+1,0,0,0)
    rp_lp = -ep*cg(1,lb,lb+1,1,mb,mb+1)*cg(1,la,lb+1,1,ma,mb+1)-em*cg(1,lb,lb+1,-1,mb,mb-1)*cg(1,la,lb+1,1,ma,mb-1)+e0*cg(1,lb,lb+1,0,mb,mb)*cg(1,la,lb+1,1,ma,mb)
    rm_lp = -ep*cg(1,lb,lb+1,1,mb,mb+1)*cg(1,la,lb+1,-1,ma,mb+1)-em*cg(1,lb,lb+1,-1,mb,mb-1)*cg(1,la,lb+1,-1,ma,mb-1)+e0*cg(1,lb,lb+1,0,mb,mb)*cg(1,la,lb+1,-1,ma,mb)
    r0_lp = -ep*cg(1,lb,lb+1,1,mb,mb+1)*cg(1,la,lb+1,0,ma,mb+1)-em*cg(1,lb,lb+1,-1,mb,mb-1)*cg(1,la,lb+1,0,ma,mb-1)+e0*cg(1,lb,lb+1,0,mb,mb)*cg(1,la,lb+1,0,ma,mb)

    return np.array([prefac_lp*rp_lp,prefac_lp*rm_lp,prefac_lp*r0_lp])

def Xb_ma_mb_e_lm(la,ma,lb,mb,eo):
    ep = -1.0/np.sqrt(2) * (eo[0]+1j*eo[1])
    e0 = eo[2]
    em = 1.0/np.sqrt(2) * (eo[0]-1j*eo[1])

    prefac_lm = np.sqrt((2*lb+1)*(2*la+1))/(2*(lb-1)+1) * cg(1,lb,lb-1,0,0,0)*cg(1,la,lb-1,0,0,0)
    rp_lm = -ep*cg(1,lb,lb-1,1,mb,mb+1)*cg(1,la,lb-1,1,ma,mb+1)-em*cg(1,lb,lb-1,-1,mb,mb-1)*cg(1,la,lb-1,1,ma,mb-1)+e0*cg(1,lb,lb-1,0,mb,mb)*cg(1,la,lb-1,1,ma,mb)
    rm_lm = -ep*cg(1,lb,lb-1,1,mb,mb+1)*cg(1,la,lb-1,-1,ma,mb+1)-em*cg(1,lb,lb-1,-1,mb,mb-1)*cg(1,la,lb-1,-1,ma,mb-1)+e0*cg(1,lb,lb-1,0,mb,mb)*cg(1,la,lb-1,-1,ma,mb)
    r0_lm = -ep*cg(1,lb,lb-1,1,mb,mb+1)*cg(1,la,lb-1,0,ma,mb+1)-em*cg(1,lb,lb-1,-1,mb,mb-1)*cg(1,la,lb-1,0,ma,mb-1)+e0*cg(1,lb,lb-1,0,mb,mb)*cg(1,la,lb-1,0,ma,mb)

    return np.array([prefac_lm*rp_lm,prefac_lm*rm_lm,prefac_lm*r0_lm])


#Definition of the norm squared of the matrix element summed over final magnetic numbers and averaged over initial polarizations.
#Be careful since this will be a vector and there are two separate contributions.
def Msqrd(omega, na, la, nb, lb):
    Ena = -0.5/na**2
    Enb = -0.5/nb**2
    #We compute the radial part of the integral on each term of the sum. 

    #start with Xa whose energy is given by Ea+omega 
    if(Ena+omega<0):
        Xa_lp1 =  DLmatElemRadintsNEl_ana(nb,lb,na,la,-1*omega,la+1)
        if(la-1>=0):
            Xa_lm1 =  DLmatElemRadintsNEl_ana(nb,lb,na,la,-1*omega,la-1)
        else:
            Xa_lm1=0
    else:
        Xa_lp1 =  DLmatElemRadintsPEl(nb,lb,na,la,omega,la+1)
        if(la-1>=0):
            Xa_lm1 =  DLmatElemRadintsPEl(nb,lb,na,la,omega,la-1)
        else:
            Xa_lm1 = 0

    #Follow with Xb that has energy Eb-omega so always below threshold.
    Xb_lp1 = DLmatElemRadintsNEl_ana(na,la,nb,lb,omega,lb+1)
    if(lb-1>=0):
        Xb_lm1 = DLmatElemRadintsNEl_ana(na,la,nb,lb,omega,lb-1)
    else:
        Xb_lm1 = 0
    
    
    tot = 0
    for ma in range(-la,la+1):
        for mb in range(-lb,lb+1):
            #vtemp_x = Xa_lp1*Xa_ma_mb_e_lp(la,ma,lb,mb,[1,0,0])+Xa_lm1*Xa_ma_mb_e_lm(la,ma,lb,mb,[1,0,0])+Xb_lp1*Xb_ma_mb_e_lp(la,ma,lb,mb,[1,0,0])+Xb_lm1*Xb_ma_mb_e_lm(la,ma,lb,mb,[1,0,0])
            #vtemp_y = Xa_lp1*Xa_ma_mb_e_lp(la,ma,lb,mb,[0,1,0])+Xa_lm1*Xa_ma_mb_e_lm(la,ma,lb,mb,[0,1,0])+Xb_lp1*Xb_ma_mb_e_lp(la,ma,lb,mb,[0,1,0])+Xb_lm1*Xb_ma_mb_e_lm(la,ma,lb,mb,[0,1,0])
            vtemp_z = Xa_lp1*Xa_ma_mb_e_lp(la,ma,lb,mb,[0,0,1])+Xa_lm1*Xa_ma_mb_e_lm(la,ma,lb,mb,[0,0,1])+Xb_lp1*Xb_ma_mb_e_lp(la,ma,lb,mb,[0,0,1])+Xb_lm1*Xb_ma_mb_e_lm(la,ma,lb,mb,[0,0,1])

            tot += np.real( np.dot(np.conjugate(vtemp_z),vtemp_z) )

    return tot*1/(2*la+1)

#Cross section
def sigma(omega,na,la,nb,lb):
    omegap = -0.5/na**2+omega+0.5/nb**2
    if(omegap<=0):
        return 0
    else:
        ao = cte.physical_constants["Bohr radius"][0]*1e2
        M2 = Msqrd(omega,na,la,nb,lb)
        return 8./3.*np.pi/(137**4) * ao**2 * omega * omegap**3 * M2

def aufromnm(nm):
    return 1/(nm*1e-7) * 1/219474.6

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Test 0 verifies that the inhomogenous equation is fulfilled for the orbital n with angular momentum l and photon energy k. To get l+1 and l-1 the variable ldl has to be changed.
#This test also outputs the function to be verify boundary conditions of the positive and negative energies case.
def test0(n,l,k):
    En = -0.5/n**2
    ldl = l-1
    if(En+k>0):
        [rt,Dlp] =  DLradfuncPE(n,l,k,ldl)
        np.savetxt("TestOutput/DLp.dat", np.column_stack(np.array([np.real(rt).astype(float),np.real(Dlp).astype(float),np.imag(Dlp).astype(float)])))
        np.savetxt("TestOutput/Intp.dat", np.column_stack(np.array([np.real(rt).astype(float),np.real(Intp).astype(float),np.imag(Intp).astype(float)])))
    #[rtn,DLn] =  DLradfuncNE(n,l,k,ldl)
    [rtna, Dlna] = DLradfuncNE_ana(n,l,k,ldl)
    #Intp = [Unl(n,l,rt[rr])*rt[rr]*Dlp[rr] for rr in range(len(rt))]
    Intna = [Unl(n,l,rtna[rr])*rtna[rr]*Dlna[rr] for rr in range(len(rtna))]
    rhs = rtna*Unl(n,l,rtna)
    
    #np.savetxt("TestOutput/DLn.dat", np.column_stack(np.array([rtn,DLn,rhs,lhs(Dlna,rtna,ldl,En,k)])))
    np.savetxt("TestOutput/Intna.dat", np.column_stack(np.array([rtna,Intna])))
    np.savetxt("TestOutput/DLa.dat", np.column_stack(np.array([rtna,Dlna,rhs,lhs(Dlna,rtna,ldl,En,k)])))

#Test 01 gives the diagonal radial integrals in a mesh of photon energies with M uniformly distributed points from km to kl for the orbital n l.
def test01(M,km,kl,n,l):
    ks = np.linspace(km,kl,M)
    f = open("TestOutput/wdepofRdint_ana.dat",'w',buffering=1)
    for kk in ks:
        f.write("%.3f %.6f %.6f %.6f %.6f \n"%(kk,np.abs(DLmatElemRadintsPEl(n,l,n,l,kk,l-1)) if(l-1>=0) else 0,DLmatElemRadintsNEl_ana(n,l,n,l,kk,l-1) if(l-1>=0) else 0,np.abs(DLmatElemRadintsPEl(n,l,n,l,kk,l+1)),DLmatElemRadintsNEl_ana(n,l,n,l,kk,l+1)))
    f.close()

#Test 1 gives table one from Sadeghpour and Dalgarno, and the static polarizability of the Hydrogen ground state.
def test1():
    ksample = np.linspace(1e-6,1e-3,50)
    mlist=[DLmatElemRadintsNEl(1,0,1,0,k,1) for k in ksample]
    popt, cov = cf(lin,ksample,mlist)
    print("Static Polarizability of Hydrogen Ground state: %.8f"%(-2/3*popt[1]))
    omega = [103.2,103.8,113.6,102.5,101.8,101.2,100.6,99.5,95.5,94.3,93.5,92.5,91.5,91.3]
    Ry = [1.82e-23,4.25e-24,4.32e-24,3.57e-20,2.50e-23,7.87e-24,3.64e-24,9.97e-25,3.62e-26,1.04e-24,3.05e-24,4.93e-24,3.75e-25,1.60e-23]
    Ra = [3.66e-24,1.48e-24,4.81e-26,4.73e-21,1.83e-24,3.11e-25,5.59e-26,7.40e-27,1.21e-25,6.14e-28,6.64e-26,2.0e-25,2.57e-26,1.29e-24]
    print("Wavelength(nm)      Rayleigh     Rayleigh(paper)        Raman        Raman(Paper)")
    for i  in range(len(omega)):
        w=omega[i]
        print("%.3f            %.2e        %.2e          %.2e        %.2e   "%(w, sigma(aufromnm(w),1,0,1,0),Ry[i], sigma(aufromnm(w),1,0,2,0),Ra[i]))
    return 0

#Test gf plots the Green function as function of energy for set values of r and rp
def testgf(l):
    emesh = np.linspace(-0.65,0.5,300)
    r = 0.8
    rp = 1.5
    file =open("TestOutput/gantest.dat","w",buffering=1)
    file2
    for e in emesh:
        if(e>0):
            file.write("%.9f  %.9f %.9f \n"%(e,float(np.real(G_pos(e,l,r,rp))),float(np.imag(G_pos(e,l,r,rp)))))
        else:
            file.write("%.9f  %.9f %.9f \n"%(e,float(np.real(G_ana(e,l,r,rp))),float(np.imag(G_ana(e,l,r,rp)))))
    file.close()
    return 0
#Test 3 computes the Rayleigh cross section for the groundstate the 2p and 3s state as well as the Raman cross section fro 1s-2s and 3s-3d. 
def test3():
    omega = np.linspace(0.05,0.75,401)
    f = open("TestOutput/rayleight1.dat","w", buffering=1)
    for w in omega:
        if(nuf(-w)%1 != 0):
            f.write(("%.8f   %.8e   %.8e   %.8e \n")%(w, sigma(w,1,0,1,0),sigma(w,2,1,2,1),sigma(w,1,0,2,0)))
    f.close()

    omega = np.linspace(0.005,0.25,401)
    f = open("TestOutput/rayleight2.dat","w", buffering=1)
    for w in omega:
        if(nuf(-w)%1 != 0):
            f.write(("%.8f   %.8e    %.8e \n")%(w, sigma(w,3,0,3,0),sigma(w,3,0,3,2)))
    f.close()
    return 0

#Test 4 generates data to be compared with figure 2,3,4,5 and 6 of McNamara et al.
def test4():
    '''omega1s = np.linspace(1e-3,np.sqrt(0.45),801)
    file = open("Cross_Sections/1sRamanRay.dat","w")
    for w in np.flip(omega1s):
        ww = 0.5-w**2
        file.write(("%.8f   %.8e    %.8e   %.8e    %.8e \n")%(ww, sigma(ww,1,0,1,0),sigma(ww,1,0,2,0),sigma(ww,1,0,3,0),sigma(ww,1,0,3,2)))
    omega1s = np.linspace(1e-3,0.25,300)
    for w in omega1s:
        ww=0.5+w
        file.write(("%.8f   %.8e    %.8e   %.8e    %.8e \n")%(ww, sigma(ww,1,0,1,0),sigma(ww,1,0,2,0),sigma(ww,1,0,3,0),sigma(ww,1,0,3,2)))
    file.close()
    
    omega2p = np.linspace(1e-3,np.sqrt(0.5/4-0.01),401)
    file = open("Cross_Sections/2pRay.dat","w")
    for w in np.flip(omega2p):
        ww = 0.5/4-w**2
        file.write(("%.8f   %.8e \n")%(ww, sigma(ww,2,1,2,1)))
    omega2p = np.linspace(1e-3,0.75-0.5/4,300)
    for w in omega2p:
        ww = 0.5/4+w
        file.write(("%.8f   %.8e \n")%(ww, sigma(ww,2,1,2,1)))
    file.close()

    omega3p = np.linspace(1e-4,np.sqrt(0.5/9-0.0001),401)
    file = open("Cross_Sections/3pRaman.dat", "w")
    for w in np.flip(omega3p):
        ww = 0.5/9-w**2
        file.write(("%.8f   %.8e   %.8e   %.8e \n")%(ww, sigma(ww,3,1,2,1), sigma(ww,3,1,4,1), sigma(ww,3,1,4,3)))
    omega3p = np.linspace(1e-6,1-0.5/9,350)
    for w in omega3p:
        ww = 0.5/9+w
        file.write(("%.8f   %.8e   %.8e   %.8e \n")%(ww, sigma(ww,3,1,2,1), sigma(ww,3,1,4,1), sigma(ww,3,1,4,3)))
    file.close()'''

    omega3s = np.linspace(1e-3,np.sqrt(0.5/9-0.0001),401)
    file = open("Cross_Sections/3sRamanRay.dat", "w")
    for w in np.flip(omega3s):
        ww = 0.5/9-w**2
        file.write(("%.8f   %.8e   %.8e \n")%(ww, sigma(ww,3,0,3,0), sigma(ww,3,0,3,2)))
    omega3s = np.linspace(1e-6,0.25-0.5/9,300)
    for w in omega3s:
        ww = 0.5/9+w
        file.write(("%.8f   %.8e   %.8e \n")%(ww, sigma(ww,3,0,3,0), sigma(ww,3,0,3,2)))
    file.close()
    return 0

