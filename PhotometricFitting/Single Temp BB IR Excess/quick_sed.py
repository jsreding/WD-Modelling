# -*- coding: utf-8 -*-
"""
Created on Tue Feb 09 11:49:43 2016
Modified on Mon Feb 18 :: 2019

@author: Erik, Josh R.
"""

import numpy as np
from bergeron_model_interp import bergeron_DA_phot, bergeron_DB_phot
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.interpolate import interp2d
from progress.bar import ChargingBar
import csv

from matplotlib import rc, rcParams
rcParams['font.family'] = 'afm'
rcParams['font.sans-serif'] = ['Helvetica']
rc('text', usetex=True)

def chisqr(obs, exp, error):
    chisqr = 0
    for n in range(len(obs)):
        # print("Deviation: ", obs[n]-exp[n])
        # print("Error: ", error[n])
        # print("Chi2: ", ((obs[n]-exp[n])**2)/(error[n]**2))
        chisqr = chisqr + ((obs[n]-exp[n])**2)/(error[n]**2)
        # print()
    return chisqr/(len(obs))

def find_nearest(array, value):
    array = np.asarray(array)
    min = array.argmin()
    idx1 = (np.abs(array[:min] - value)).argmin()
    idx2 = (np.abs(array[min:] - value)).argmin()+min
    print(idx1, idx2)
    return min, idx1, idx2

plt.close('all')

mass_mod = []
teff_mod = []
rad_mod = []

with open("/home/jsreding/Documents/UNC/Research/Projects/SpottedWDs/COModel_ThinH.csv", 'r') as f:
    reader = csv.DictReader(f, delimiter=',')
    for row in reader:
        mass_mod.append(float(row['Mass']))
        teff_mod.append(float(row['Teff']))
        rad_mod.append(float(row['rad']))

mrr = interp2d(mass_mod, teff_mod, rad_mod, kind='cubic')

lg = np.arange(7.5, 8.51, 0.01)
T = np.arange(8000., 9001., 1.)

Teff = 8369.
logg = 8.13
plx = 12.9437 #parallax in mas
# chisq = 9999999.
# chi = np.zeros(len(lg))
# rad = 0.
# bar = ChargingBar("Processing...", max=len(lg)*len(T))
# for n in range(len(T)):
#     Teff = np.round(T[n], 0)
#     for i in range(len(lg)):
#         logg = np.round(lg[i], 2)

mod_data = bergeron_DA_phot(logg,Teff)
mod_mags = mod_data[6::]

age = mod_data[5]/(10**6)
mass = mod_data[2]
# print("Mass: ", mass)

sig_sfb = 5.67040e-5 # stefan-boltzmann in erg cm^-2 K^-4 s^-1
Mbol = mod_data[3]
Mbol_sol = 4.75
L_sol = 3.826e33 # in ergs/s
R_sol = 6.9599e10 # in cm

# Disk relevant constants
clight = 2.997e10 # cm s^-1
hpl = 6.626e-27 #erg s
kB = 1.381e-16 #erg K^-1

L_wd = L_sol*(10**(0.4*(Mbol_sol-Mbol))) #luminosity of wd in ergs/s
R_wd = mrr(mass, Teff)[0] #returns R_wd in cm
# print("Radius: ", R_wd)

# mags    FUV, NUV,j_U,j_B,j_V,j_R,j_I,j_J   ,j_H   ,jon_K ,v_j,v_H,v_K,sds_u,sds_g,sds_r,sds_i,sds_z, W1, W2, W3, W4
EC_mags = [999,999,999,999,999,999,999,999,999,999,999,999,999,999,17.503,17.45,17.511,17.58,999,999,999,999]
EC_err =  [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.005,0.005,0.007,0.019,0.0,0.0,0.0,0.0]

# pivot_wav given as gal_FUV, gal_NUV, jon_U,jon_B,jon_V,jon_R,jon_I,2m_J,2m_H,2m_K,sdss_u,sdss_g,sdss_r,sdss_i,sdss_z, W1, W2, W3, W4
phot_pivot_wav = np.array([0.1528,0.2271,0.36,0.44,0.55,0.7,0.9,1.25,1.65,2.15,1.25,1.65,2.15,0.3543,0.477,0.6231,0.7625,0.9134,3.35,4.6,11.56,22.24])
# photometric zero pts given in janskys
#phot_zero_pts = np.array([3631,3631,1810,4260,3540,2870,2250,1593.,1024.,667.,1531.8,1024,667.43631,3631,3631,3631,3631,3631,306.681,170.663,29.0448,8.2839])
phot_zero_pts = np.array([3631.,3631.,1810.,4260.,3540.,2870.,2250.,1530.,1019.,631.,1531.8,1024,667.43631,3631.,3631.,3631.,3631.,3631.,306.681,170.663,29.0448,8.2839])
upper_zero_pts = np.array([29.0448,8.2839])

vista = False
if vista:
    mod_mags_add = np.array([mod_mags[7] - 0.077*(mod_mags[7] - mod_mags[8]),mod_mags[9] + 0.010*(mod_mags[7] - mod_mags[9])])
    mod_mags = np.concatenate((mod_mags,mod_mags_add))

D_wd = 1./(plx/1000.)

mod_iso_wav = np.array([0.1528,0.2271,0.3971,0.4481,0.5423,0.6441,0.8071,1.2350,1.6620,2.1590,1.2350,1.6620,2.1590,0.3146,0.4670,0.6156,0.7471,0.8918,3.3526,4.6028,11.5608,22.0883])
mod_corr = np.array([0,0,0.0915,0.0069,0.00,0.0018,-0.0014,-0.0140,0.0060,0.0080,-0.0140,0.0060,0.0080,0.0424,-0.0023,-0.0032,-0.0160,-0.0276,0,0,0,0])
mod_zero_pts = [3631.,3631.,1649.,4060.,3723.,3168.,2459.,1594.,1024.,666.8,1594.,1024.,666.8,3631.,3631.,3631.,3631.,3631.,306.681,170.663,29.0448,8.2839]

distmod = 5.*np.log10(D_wd)-5
EC_flux = (10**(-0.4*(np.array(EC_mags))))*phot_zero_pts
EC_flux_err = np.log(10)*0.4*np.array(EC_err)*EC_flux
mod_flux = (10**(-0.4*(np.array(mod_mags+distmod) + mod_corr)))*mod_zero_pts

scale = []
scale_err = []
for l in range(len(EC_mags)):
    if EC_err[l] == 0:
        pass
    else:
        if phot_pivot_wav[l] > 2.16 or phot_pivot_wav[l] < 0.25:
            pass
        else:
            scale.append(EC_flux[l]/mod_flux[l])
            scale_err.append(EC_flux_err[l])

flux_scale = np.average(scale,weights=1.0/(np.array(scale_err)**2.0))

sort_inds = np.argsort(mod_iso_wav)
EC_flux = EC_flux[sort_inds]
mod_iso_wav = mod_iso_wav[sort_inds]
EC_flux_err = EC_flux_err[sort_inds]
mod_flux = mod_flux[sort_inds]

plot_wav = np.linspace(0.1,25,5000)
plot_nubb = clight*(10.0**4.0)/plot_wav

phot_inds = np.where(EC_flux_err)[0]
EC_flux_pts = EC_flux[phot_inds]
EC_flux_err_pts = EC_flux_err[phot_inds]
wd_mod_pts = mod_flux[phot_inds]#*flux_scale
iso_wav_pts = mod_iso_wav[phot_inds]
iso_nu_pts = (clight*(10**4))/iso_wav_pts

D_wd_wdrad = (D_wd)*(3.08e18)/R_wd # Returns the distance to the wd in wd radii
def bbody_only(nu,T,R):
    return ((R/(D_wd_wdrad))**2.0)*((2.0*np.pi*hpl*(nu**3.0))/(clight**2.0))*(1.0/np.expm1(hpl*nu/(kB*T)))*(10**23.0)
bbody_scale = []
wd_bbody = bbody_only(iso_nu_pts,Teff,1.0)

for l in np.arange(2,len(wd_bbody)):
    b_scale = wd_mod_pts[l]/wd_bbody[l]
    bbody_scale.append(b_scale)
bbody_scale = np.mean(bbody_scale)

#         chi2test = chisqr(EC_flux_pts, wd_mod_pts, EC_flux_err_pts)
#         chi[i] = chi2test
#         bar.next()
#         if chi2test < chisq:
#             print()
#             chisq = chi2test
#             print(Teff, logg, mass, R_wd, chisq)
#             bestT = Teff
#             bestlg = logg
#             rad = R_wd
# print(bestT, bestlg, rad)
# bar.finish()
#
# plt.figure()
# plt.plot(lg, chi)
# plt.show()
#
# print("Chi square: ", chi.min())
#
# mind, ind1, ind2 = find_nearest(chi, chi.min()+np.sqrt(2./7.))
# print("Best: ", lg[mind], "+/-", np.mean([lg[mind]-lg[ind1], lg[ind2]-lg[mind]]))
# print(iso_wav_pts[3], iso_wav_pts[5])
# print(iso_nu_pts[3], iso_nu_pts[5])
# print(EC_flux_pts[3], EC_flux_pts[5])
# print("Spectral index = ", np.log(EC_flux_pts[3]/EC_flux_pts[5])/np.log(iso_nu_pts[5]/iso_nu_pts[3]))

plt.figure(1)
plt.clf()
plt.xlim(750., 320000.)
plt.xscale('log')
# plt.yscale('log')
# plt.title(r"DA No-NUV/u Best Fit ($\chi^2\approx1.9$)", fontsize=24)
plt.ylabel("Flux (mJy)", fontsize=20)
plt.xlabel(r"Wavelength (\AA)", fontsize=20)
plt.errorbar(iso_wav_pts*10000.,1000.0*EC_flux_pts,yerr=1000.0*EC_flux_err_pts,fmt='s',mfc='none',color='b',label='Observed Photometry', zorder=1)
#plt.errorbar(upper_wav_pts,1000.0*EC_flux_upper,yerr=EC_flux_upper_err,fmt='s',mfc='none',color='b',uplims=True)
plt.plot(iso_wav_pts*10000.,1000.0*wd_mod_pts,'k*',label='Bergeron Model Phot',linewidth=1.0, zorder=1)
plt.plot(plot_wav*10000.,bbody_only(plot_nubb,Teff,1.0)*1000.0,'k--',label=r'Blackbody Model ($\mathrm{T_{eff}}=%s$K, log g$=%s$)'%(Teff, logg), zorder=1)#*bbody_scale
# plt.axvspan(0, 2900., alpha=0.25, color='c')
# plt.axvspan(2900., 11500., alpha=0.25, color='y')
# plt.axvspan(11500., 320000., alpha=0.25, color='m')
plt.legend(fontsize=14, loc=7)
plt.show()

np.save('8369kbb.npy', [plot_wav*10000., bbody_only(plot_nubb,Teff,1.0)*1000.0])
