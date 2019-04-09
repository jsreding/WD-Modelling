# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 09:05:56 2017

@author: Erik
"""

import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from ec_table_tools import ec_wdfitted_flux
from astropy.table import Table 
import glob,os

plt.ioff()

sig_sfb = 5.67040e-5 # stefan-boltzmann in erg cm^-2 K^-4 s^-1
solar_lum = 3.826e33 # in ergs/s
solar_radius = 6.9599e10 # in cm

# Disk relevant constants
clight = 2.997e10 # cm s^-1
hpl = 6.626e-27 #erg s
kB = 1.381e-16 #erg K^-1
parsec = 3.0857e18 #c
ec_table = Table.read('EC_Tablev3r3.csv')
ec_ir_excess = np.where(ec_table['excess flag'])[0]

T_IR = np.zeros(len(ec_table))
R_IR = np.zeros(len(ec_table))
R_wd = np.zeros(len(ec_table))
age_wd = np.zeros(len(ec_table))
m_wd = np.zeros(len(ec_table))

dir_name = 'bb_plots/'

files = glob.glob(dir_name+'/*')
for f in files:
    os.remove(f)

for test_ind in ec_ir_excess:
    ec_name = ec_table[test_ind]['ec name']
    fit_wavs,ec_fit_flux,ec_fit_eflux,mod_fit_flux_scaled,all_mod_wavs,all_mod_flux_scaled,scale,age,mass,R_wd_cm = ec_wdfitted_flux(ec_table,test_ind)
    ir_wavs = np.where(fit_wavs > 0.55)[0]
    fit_nu = clight*(10.0**4.0)/fit_wavs
    all_nu = clight*(10.0**4.0)/all_mod_wavs

    R_wd[test_ind] = R_wd_cm
    age_wd[test_ind] = age
    m_wd[test_ind] = mass
    
    dist = 10.*np.sqrt(1./scale) # Gives us distance in parsecs
    
    wd_radius = R_wd_cm # wd radius in cm
    dist_cm = dist*parsec # distance in cm
    dist_wd_rad = dist_cm/wd_radius # distance in wd radii

    @jit
    def bbody_only(nu,temp,radius,dist):
        return ((radius/(dist))**2.0)*((2.0*np.pi*hpl*(nu**3.0))/(clight**2.0))*(1.0/np.expm1(hpl*nu/(kB*temp)))*(10**23.0)
        # returns a blackbody monochromatic  flux at a distance in janskys, distance 
        # must be given in wd_radii
    
    bbody_Tarray = np.linspace(500,5500,401)
    bbody_Rarray_solar = np.linspace(0.01,1.2,239)
    bbody_Rarray = bbody_Rarray_solar/(wd_radius/solar_radius)

    def bbody_chi_sqr(Tarray,Rarray):
        sumflux = bbody_only(fit_nu[ir_wavs],Tarray,Rarray,dist_wd_rad) + mod_fit_flux_scaled[ir_wavs]
        diff = (sumflux - ec_fit_flux[ir_wavs])
        red_chi_sqr = np.sum((diff**2.)/(ec_fit_eflux[ir_wavs]**2.))#/(len(diff)-4.0)
        return red_chi_sqr    

    bbody_chi_array = np.zeros((len(bbody_Tarray)*len(bbody_Rarray),3))
    niter = 0
    for i in range(len(bbody_Tarray)):
        for j in range(len(bbody_Rarray)):
            chi_sqr = bbody_chi_sqr(bbody_Tarray[i],bbody_Rarray[j])
            bbody_chi_array[niter,0] = chi_sqr
            bbody_chi_array[niter,1] = bbody_Tarray[i]
            bbody_chi_array[niter,2] = bbody_Rarray[j]
            niter+=1
    
    best_ind = np.where(bbody_chi_array[:,0] == bbody_chi_array[:,0].min())[0]
    best_T = bbody_chi_array[best_ind,1]
    best_R = bbody_chi_array[best_ind,2]*wd_radius/solar_radius # gives radius in solar radii
    T_IR[test_ind] = best_T
    R_IR[test_ind] = best_R    
    print 'Now fitting: %s'%ec_name
    print 'Best fitting blackbody has T %s K and R %s solar' %(np.round(best_T),np.round(best_R,decimals=2))

    wd_Teff = ec_table[test_ind]['wd teff']
    wd_logg = ec_table[test_ind]['wd logg']
    sort_inds = np.argsort(all_mod_wavs)
    name = ec_table[test_ind]['ec name']
    otype = ec_table[test_ind]['sptype']

    plt.figure()
    plt.clf()
#    plt.title(name)#+' Type:%s'%otype + ' IR Teff: %s IR Radius %s' %(np.round(best_T[0]),np.round(best_R[0],decimals=2)))
    plt.xscale('log')
    plt.yscale('log')
    plt.errorbar(fit_wavs,ec_fit_flux,yerr=ec_fit_eflux,fmt='o',color='b',ecolor='b',markeredgecolor='b',markersize=2)
    plt.plot(all_mod_wavs[sort_inds],all_mod_flux_scaled[sort_inds],'r--',label='WD Model: Teff=%s$\,$(K) Logg=%s'%(np.round(wd_Teff,decimals=1),np.round(wd_logg,decimals=1)))  
    plt.plot(all_mod_wavs[sort_inds],all_mod_flux_scaled[sort_inds]+bbody_only(all_nu[sort_inds],best_T,best_R/(wd_radius/solar_radius),dist_wd_rad),'g-',label='WD Model + Blackbody: T=%s$\,$(K) R=%s$\,$R$_\odot$'%(best_T[0],best_R[0]))
    plt.xlabel('Wavelength ($\mu$m)')
    plt.xlim(0.1,11.)
    plt.ylabel('Flux (Jy)')
    plt.annotate(s='EC '+name,xy=(0.4,0.9),xycoords='figure fraction')
    plt.legend(loc=3)
    fig1 = plt.gcf()
    fname = dir_name +'/'+ name +'.jpg'
#    fig1.savefig(fname)
#    plt.close()
    fname = dir_name + name +'.eps'
#    fig1.savefig(fname,format='eps',bbox_inches='tight')
    plt.close()

#tir_col = Table.Column(T_IR,name='bb teff')
#rir_col = Table.Column(R_IR,name='bb radius')
#rwd_col = Table.Column(R_wd,name='wd radius')
#mwd_col = Table.Column(m_wd,name='wd mass')
#age_col = Table.Column(age_wd,name='wd age')
#ec_table['wd radius'] = rwd_col
#ec_table['wd age'] = age_col
#ec_table['wd mass'] = mwd_col
#ec_table['bb teff'] = tir_col
#ec_table['bb radius'] = rir_col
#ec_table.write('EC_Tablev3r4.csv',format='csv')