''' 
Program to calculate bias in 6dFGSv Bulk Flow from lognormal velocities
Morag I. Scrimgeour, 19/07/2015
'''

import numpy as np
import time
import cosmologycodes as cc
from usefulcodes import optbinsize
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import sys
speedlight = 299792.458

def MLE_BF(x_arr, y_arr, z_arr, r_arr, v_arr, sigma_n, sigmastar2):
    ngal = len(x_arr)
    Aij = np.zeros([3,3])
    rihat = np.array([x_arr/r_arr, y_arr/r_arr, z_arr/r_arr])
    for i in range(3):
        for j in range(3):
            Aij[i,j] = np.sum(rihat[i,n]*rihat[j,n] / (sigma_n[n]**2. + sigmastar2) for n in range(ngal))
    Aij_inv = np.linalg.pinv(Aij)
    rsum = np.array([np.sum(rihat[0,n]*v_arr[n]/(sigma_n[n]**2.+sigmastar2) for n in range(ngal)),
                     np.sum(rihat[1,n]*v_arr[n]/(sigma_n[n]**2.+sigmastar2) for n in range(ngal)),
                     np.sum(rihat[2,n]*v_arr[n]/(sigma_n[n]**2.+sigmastar2) for n in range(ngal))])
    u = np.dot(Aij_inv,rsum)
    return u

def readin6dFGSvData():

    # -------------------- Read in 6dFGSv data ----------------------
    Dz = []
    xerr = []
    zo = []
    xobs = []
    x_arr = []
    y_arr = []
    z_arr = []

    f = open('/Users/Morag/6dFGS/6dFGSv/data/6dFGSv_vpecs_meanx_sigmax.dat','r')
    #numgal:  8885
    #  meanx  MLx  stddevx
    line = f.readline()
    p = line.split()
    ngal = int(p[1])
    for line in f:
        if not line.startswith('#'):
            p = line.split()
            xobs.append(float(p[0]))
            xerr.append(float(p[2]))
    f.close()
    f = open('/Users/Morag/6dFGS/6dFGSv/data/6dFGSv_info.dat','r')
    #numgal:  8885
    #  RA  Dec  gl  gb  sgl  sgb  cz  zshift  Dz  x_Eq  y_Eq  z_Eq  x_Gal  y_Gal  z_Gal
    for line in f:
        if not line.startswith('#'):
            p = line.split()
            Dz.append(float(p[8]))
            x_arr.append(float(p[9]))
            y_arr.append(float(p[10]))
            z_arr.append(float(p[11]))
    f.close()

    print 'Number galaxies: %i' % (ngal)

    # Convert lists to numpy arrays
    Dz = np.array(Dz)
    xerr = np.array(xerr)
    x_arr = np.array(x_arr)
    y_arr = np.array(y_arr)
    z_arr = np.array(z_arr)
    xobs = np.array(xobs)
    r_arr = np.sqrt(x_arr**2 + y_arr**2 + z_arr**2)

	return Dz, xerr, x_arr, y_arr, z_arr, xobs, r_arr

def zeroVelocityMocks():

    sigmastar2 = 250.**2

    # -------------------- Read in 6dFGSv data ----------------------

    Dz, xerr, x_arr, y_arr, z_arr, xobs, r_arr = readin6dFGSvData()

    #---------------------- Convert Dz to z -------------------------   

    zfile = '/Users/Morag/Cosmology_codes/redshift_OM0.3175.dat'
    zo = cc.z_x_arr(Dz, zfile)

    #-------- Assign each galaxy zero x, perturbed by xerr ----------   
    # Give all galaxies v=0 but perturb x by xerr
    # Generate array of num random deltax values with std dev = xerr[gal]
    
    meanx = xerr*np.random.standard_normal(ngal)

    # Calculate derived Dr
    Dr_der = Dz / 10**meanx

    # Calculate new derived redshifts
    zr_der = cc.z_x_arr(Dr_der, zfile)

    vp_rand = speedlight * ((zo - zr_der)/(1.+zr_der))
    sigma_v = 0.324 * 100. * Dz

    #------------------------- Calculate BF ---------------------------

    u = MLE_BF(x_arr, y_arr, z_arr, r_arr, vp_rand, sigma_v, sigmastar2)

    print u

    return

def main():
    zeroVelocityMocks()
    return

if __name__ == "__main__":
    main()