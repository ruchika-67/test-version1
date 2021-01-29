import numpy as np
from scipy.integrate import quad
from scipy.interpolate import splrep, splev, UnivariateSpline
data=np.loadtxt("fs8z_rec.dat")
z=data[:,0]
fs8=data[:,1]
sigma=data[:,2]
print(len(z))

delp=[]
delta=[]
sigma8_0= 1
del_0 =1
om0= 0.3
print('ruchika')
delp_0 =1

deltap = -(del_0/sigma8_0)*(fs8/(1.+z))


def delp(z1):
    delp1 = UnivariateSpline(z,deltap, k=3, s=0)
    return delp1(z1)
print(delp(1))
############################################################chk after that

def funcfs8(z1):#making spline of integral for equatin 12
 fs81= UnivariateSpline(z, fs8/(1.+z), k=3, s=0)
 return fs81(z1)

print(funcfs8(1.1))
#simply integrating and getting equation 12
for i in range(len(z)):
 int1 = quad(funcfs8, 0, z[i])[0]
 delta1= del_0 -(del_0/sigma8_0)*int1
 delta.append(delta1)

print(len(delta))

def intfunc(z1):
    int1 = UnivariateSpline(z, (delta*(-deltap)/(1.+z)), k=3, s=0)
    int2 = quad(int1,0,z)[0]                        
    return int2(z1)

print(intfunc(1.1))

firstterm= (1.+z)**2.*((delp_0**2./(delp(z)**2.))
#print(firstterm[0])
                       
secondterm = -3.*om0*((1.+z)**2./(delp(z)**2.))*intfunc(z)
                       
def Hubble Constant(z):
 return = np.sqrt(first term +second term)                       
                       