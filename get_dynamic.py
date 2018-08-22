#this is to calculate the dynamic matrix from the force constants

import numpy as np
import math as mh
import cmath as ch
from numpy import linalg as LA

#information about the unit cell this will be modified later to read all the information from POSCAR
n_atom=2
n_mass=12
nx=5
ny=5
N_atom=n_atom*nx*ny
# constants needed
eV = 1.602E-19;
AA = 1E-10;
amu = 1.660538921E-27;  # atomic mass unit (in kg)
hbar = 1.054571726E-34; # Planck constant
mass_c= n_mass*amu;     # mass of C12 isotope atom
pi=mh.pi
#mass_C13 = 13.0*amu;    % mass of C13 isotope atom
#a_cc = 1.46*Angstrom;




# let us worry about the unit all_atoms_in_supercell
f1=open('POSCAR','r+')
f3=open('SPOSCAR','r+')
f2=open('FORCE_CONSTANTS_2ND','r+')


# get the information about the unit cells
lines1=f1.readlines()
a1=np.fromstring(lines1[2],sep=' ')*AA
a2=np.fromstring(lines1[3],sep=' ')*AA
a3=np.fromstring(lines1[4],sep=' ')*AA

V=np.dot(np.cross(a1,a2),a3)
# calculate the unit vector in reciprocal space
b1=2*pi*np.cross(a2,a3)/V
b2=2*pi*np.cross(a3,a1)/V


# calculate and prepare the FORCE_CONSTANTS   (!!worry about the unit latter)
Nc=N_atom*3
ss=(Nc,Nc)
FC=np.zeros(ss)
lines2=f2.readlines()

for i in range(N_atom):
    for j in range(N_atom):
        nl=4*(i*N_atom+j)+2
        FC[i*3:(i+1)*3,j*3:(j+1)*3]=np.fromstring(lines2[nl]+lines2[nl+1]+lines2[nl+2],sep=' ').reshape(3,3)

#print(FC-FC.T)
#print(FC.shape)
# output the force constant to file to check
#f=open('FC.dat','w+')
#f.writelines(FC)
np.savetxt('kapp.dat',FC[:9,:9],fmt='%.5e', delimiter = '\t')

##constructure a matrix

#read the atomic coordinates
lines3=f3.readlines()
a1=np.fromstring(lines3[2],sep=' ')*AA
a2=np.fromstring(lines3[3],sep=' ')*AA
a3=np.fromstring(lines3[4],sep=' ')*AA
aa=np.concatenate((a1,a2,a3),axis=0).reshape(3,3).T

#print(aa)

ss=(N_atom,3)
xxs=np.zeros(ss)

for i in range(N_atom):
    xxs[i,:]=np.fromstring(lines3[i+7],sep=' ')
xxs=aa.dot(xxs.T).T

#print(xxs.shape)
#print(xxs)

#Specify the q coordinate and get the dynamic matrix
ss=(n_atom*3,n_atom*3)
D=np.zeros(ss,dtype=complex)
nn=N_atom//2

#Specify the number of q qpoints

#print(nn)

for i in range(nn):
    D[:3,:3]=D[:3,:3]+FC[:3,i*3:(i+1)*3]*ch.exp(1j*(q*b2@(xxs[i,:]-xxs[0,:])))/mass_c*eV/(AA*AA)
    D[:3,3:6]=D[:3,3:6]+FC[:3,(i+nn)*3:(i+1+nn)*3]*ch.exp(1j*(q*b2@(xxs[i+nn,:]-xxs[0,:])))/mass_c*eV/(AA*AA)
    D[3:6,:3]=D[3:6,:3]+FC[3*nn:3*nn+3,i*3:(i+1)*3]*ch.exp((1j*q*b2@(xxs[i,:]-xxs[nn,:])))/mass_c*eV/(AA*AA)
    D[3:6,3:6]=D[3:6,3:6]+FC[3*nn:3*nn+3,(i+nn)*3:(i+1+nn)*3]*ch.exp(1j*(q*b2@(xxs[i+nn,:]-xxs[nn,:])))/mass_c*eV/(AA*AA)


w,v=LA.eig(D)
w0=0.1*np.amax(np.abs(w))

w,v=LA.eig(D+w0*np.identity(6))
print(np.sort(np.sqrt(w.real-w0))/(2*mh.pi))
#print(v)
