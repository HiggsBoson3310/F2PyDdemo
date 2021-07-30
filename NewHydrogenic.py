import coul92 as cl
import couln as cln
import gensub as gs
import numpy as np
import scipy.special as spc
import math 
#import Hydrogenic_functions as Hf
import time

def k(e):
    return np.sqrt(2*e)
def kappa(e):
	return np.sqrt(-2*e)
def etaf(e,l):
	return 1/k(e) * np.log(2*k(e)) + np.angle(spc.gamma(l+1-1j/k(e)))-1/2. * l *np.pi
def nuf(e):
	return np.sqrt(1/(-2*e))
def D(e,l):
	if(spc.gamma(nuf(e)-l)>0):
		return np.sqrt(np.pi) * (2/nuf(e))**(nuf(e)) * 1/np.sqrt(spc.gamma(nuf(e)-l)*spc.gamma(l+1+nuf(e)))
	else:
		return np.sqrt(np.pi) * (2/nuf(e))**(nuf(e)) * 1/np.sqrt((-1*spc.gamma(nuf(e)-l))*spc.gamma(l+1+nuf(e)))

def F(l,eta,rho):
    [f,*args] = cl.coul90(rho,eta,0,l,0,0)
    return f[-1]

def Fp(l,eta,rho):
    [f,g,fp,*args] = cl.coul90(rho,eta,0,l,0,0)
    return fp[-1]

def G(l,eta,rho):
    [f,g,*args] = cl.coul90(rho,eta,0,l,0,0)
    return -1*g[-1]

def Gp(l,eta,rho):
    [*args,gp] = cl.coul90(rho,eta,0,l,0,0)
    return -1*gp[-1]



#Regular Coulomb function for an attractive potential and positive energy
def Ucl(e,l,r):
	if(k(e)*r > 1):
		return np.sqrt(2/(np.pi*k(e)))*F(l,-1/k(e),r*k(e))
	else:
		return gs.seaton(l,e*2.0,r,1.0)[0]

#Long Range behavior
def LRUcl(e,l,r):
	return np.sqrt(2/(np.pi*k(e))) * np.sin(k(e)*r+1/k(e) * np.log(r) + etaf(e,l)) 

#And its derivative
def Uclp(e,l,r):
	if(k(e)*r>1):
		return np.sqrt(2/(np.pi*k(e)))*Fp(l,-1/k(e),r*k(e))*k(e)
	else:
		return gs.seaton(l,e*2.0,r,1.0)[1]
 

#Regular Coulomb function for an attractive potential and negative energy
def Ucln(e,l,r):
	nu = nuf(e)
	if(nu>l):
		return gs.seaton(l,e*2.0,r,1.0)[0]
	else:
		return -1*gs.seaton1(l,e*2.0,r,1)[0]

#Analytic coulomb function for negative energies and nu<l
def Uclno(e,l,r):
	return gs.coulfg(l,e*2.0,r,1e-14)[0]
def Uclnop(e,l,r):
	return gs.coulfg(l,e*2.0,r,1e-14)[1]
def LRUclno(e,l,r):
	nu = nuf(e)
	return nu**(l+1+nu)/spc.gamma(l+1-nu) * (1/(2*r))**nu * np.exp(r/nu)


#Long Range behavior
def LRUcln(e,l,r):
	nu = nuf(e)
	return 1.0/np.sqrt(np.pi*kappa(e)) * np.sin(np.pi*(nuf(e)-l)) * 1.0/D(e,l) * r**(-1.0*nuf(e)) * np.exp(kappa(e)*r)
	

#And its derivative
def Uclpn(e,l,r):
	nu = nuf(e)
	if(nu>l):
		return gs.seaton(l,e*2.0,r,1.0)[1]
	else:
		return -1*gs.seaton1(l,e*2.0,r,1.0)[1]

#Irregular Coulomb function for an attractive potential and positive energy
def Ccl(e,l,r):
	if(k(e)*r>1):
		return np.sqrt(2/(np.pi*k(e)))*G(l,-1/k(e),r*k(e))
	else:
		return gs.seaton(l,e*2.0,r,1.0)[2]

#And its derivative.
def Cclp(e,l,r):
	if(k(e)*r>1): 
		return np.sqrt(2/(np.pi*k(e)))*Gp(l,-1/k(e),r*k(e))*k(e)
	else:
		return gs.seaton(l,e*2.0,r,1.0)[3]


#Long Range behavior
def LRCcl(e,l,r):
	return -np.sqrt(2/(np.pi*k(e))) * np.cos(k(e)*r+1/k(e) * np.log(r) + etaf(e,l))


#Regular Coulomb function for an attractive potential and negative energy
def Ccln(e,l,r):
	nu = nuf(e)
	if(nu>l):
		return gs.seaton(l,e*2.0,r,1.0)[2]
	else:
		return -1*gs.seaton1(l,e*2.0,r,1.0)[2]

#Long range behavior
def LRCcln(e,l,r):
	nu = nuf(e)
	return -1.0/np.sqrt(np.pi*kappa(e)) * np.cos(np.pi*(nuf(e)-l)) * 1.0/D(e,l) * r**(-1.0*nuf(e)) * np.exp(kappa(e)*r)

#And its derivative
def Cclpn(e,l,r):
	nu = nuf(e)
	if(nu>l):
		return gs.seaton(l,e*2.0,r,1.0)[3]
	else:
		return -1*gs.seaton1(l,e*2.0,r,1.0)[3]
#Hydrogenic orbitals
def Unl(n,l,r):
    return r*np.sqrt(math.factorial(n-l-1)/(2*n*math.factorial(n+l)))*spc.assoc_laguerre(2*r/n,n-l-1,2*l+1)*np.exp(-r/n)*(r**l)*((2/n)**(l+3./2.))
#Closed channel function normalized to unit norm at zero quantum defect
def W(nu,l,r):
	wa1 = np.zeros((20,10))
	wa2 = np.zeros((20,10))
	wa3 = np.zeros((20,10))
	wa4 = np.zeros((20,10))
	W = cln.couln(l,1,-1/nu**2,r,1e-12,wa1,wa2,wa3,wa4,0,70)
	gam = spc.gamma(nu-l)
	if gam>0:
		prefac = nu**(3/2.) * 1/np.sqrt(nu**2 * spc.gamma(nu+l+1) * gam)
	else:
		prefac = nu**(3/2.) * 1/np.sqrt(nu**2 * spc.gamma(nu+l+1) * (-1*gam))
	return prefac*W[0]

def LRW(nu,l,r):
	return 1*np.sqrt(nu/np.pi) * D(-0.5/nu**2,l) * r**(nu) * np.exp(-r/nu)

#Scaled Whittaker function

def W_scal(nu,l,r):
	wa1 = np.zeros((20,10))
	wa2 = np.zeros((20,10))
	wa3 = np.zeros((20,10))
	wa4 = np.zeros((20,10))
	W = cln.couln(l,1,-1/nu**2,r,1e-12,wa1,wa2,wa3,wa4,0,70)
	return W[0]

def LRW_scal(nu,l,r):
	return np.exp(-r/nu) * (2*r/nu)**nu

def Wp(nu,l,r):
	wa1 = np.zeros((20,10))
	wa2 = np.zeros((20,10))
	wa3 = np.zeros((20,10))
	wa4 = np.zeros((20,10))
	W = cln.couln(l,1,-1/nu**2,r,1e-12,wa1,wa2,wa3,wa4,0,70)
	return nu**(3/2.) * 1/np.sqrt(nu**2 * spc.gamma(nu+l+1) * spc.gamma(nu-l))*W[1]
def W_scalp(nu,l,r):
	wa1 = np.zeros((20,10))
	wa2 = np.zeros((20,10))
	wa3 = np.zeros((20,10))
	wa4 = np.zeros((20,10))
	W = cln.couln(l,1,-1/nu**2,r,1e-12,wa1,wa2,wa3,wa4,0,70)
	return W[1]
#Seaton-Burgess cutoff
def SB(l,r):
    return (1-np.exp(-10*r/(l*(l+1))))**(2*l+1)

#to test the new module compare with hydogenic orbitals for 1s, 4f and 6p
#r = np.linspace(1e-6,5*4**2,900)
#tot = np.zeros((900,4))
#for i in range(900):
#	tot[i] = [r[i],Unl(4,3,r[i]),W(4,3,r[i]),Wp(4,3,r[i])]
#np.save("1scom.npy",tot)
#print(nu)
#r = np.zeros(401)
#wtest = np.zeros((401,2))
#for i in range(401):	
#	wtest[i,0]=(Ucln(eau,l,i*5*8**2/400)*Wp(nu,l,i*5*8**2/400)-Uclpn(eau,l,i*5*8**2/400)*W(nu,l,i*5*8**2/400))-(-np.sin(np.pi*(nu-l))*2/np.pi)
#	wtest[i,1]=Ucln(eau,l,i*5*8**2/400)*Cclpn(eau,l,i*5*8**2/400)-Uclpn(eau,l,i*5*8**2/400)*Ccln(eau,l,i*5*8**2/400)-(2/np.pi)
#np.save("wtest.npy",np.array(wtest))
def Wtest():
	l = 5
	zion = 1.0
	eau = -0.5/(4**2)
	print(nuf(eau))
	nu = nuf(eau)
	r = np.linspace(1e-5,150,601)
	fogt = np.zeros((601,4))
	for i in range(601):
		[f,fp,g,gp] = [Uclno(eau,l,r[i]),Uclnop(eau,l,r[i]),W_scal(nuf(eau),l,r[i]),W_scalp(nuf(eau),l,r[i])]
		#[f,fp,g,gp,*args]
		fogt[i] = [r[i],f*gp-fp*g,-nu**l * 2/spc.factorial(l-nu),np.sqrt(nuf(eau))]
	np.savetxt("TestOutput/ftest.dat",fogt)



#Wtest()
#l=1
#eau = -0.5*(-1/1**2)
#np.sqrtnu = np.sqrt(-1/(2*eau))
#r = np.linspace(0,50,401)
#fogt = np.zeros((401,6))
#fgt = np.zeros((401,4))
#for i in range(401):
#	fogt[i] = gs.fogo(eau*2.0,l,r[i],1e-12)
#	fgt[i] = gs.seaton(l,eau*2.0,r[i],1)
#np.save("fogotbnu.npy",fogt)
#np.save("fgtbnu.npy",fgt)'''