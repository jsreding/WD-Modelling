import numpy as np
from bergeron_model_interp import bergeron_DA_phot,bergeron_DB_phot
import matplotlib.pyplot as plt
from numba import jit
from scipy import interpolate
import matplotlib.patches as mpatches

from matplotlib import rc, rcParams
rcParams['font.family'] = 'afm'
rcParams['font.sans-serif'] = ['Helvetica']
rc('text', usetex=True)

"""Start with a simple SED calculation"""
plt.close('all')

Teff = 9104.
logg = 9.25
plx=12.9437 #parallax in mas
db=False
obj_name='J1252-0234'

if db:
    mod_data = bergeron_DB_phot(logg,Teff)
    mod_mags = mod_data[6::]
else:
    mod_data = bergeron_DA_phot(logg,Teff)
    mod_mags = mod_data[6::]

age = mod_data[5]/(10**6)
mass = mod_data[2]

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
R_wd = np.sqrt(L_wd/(4*np.pi*sig_sfb*(Teff**4.0))) #returns R_wd in cm

#Fuv,eFuv,NUV,eNUV
Galex=[22.682,0.218,19.232,0.024]

#U,B,V,R,I
Johnson_opt=[999,0.0,999,0.0,999,0.0,999,0.0,999,0.0]

#2mass J,H,K
Two_mass=[999,0.0,999,0.0,999,0.0]

#SDSS u,g,r,i,z
sdss=[18.036,0.012,17.503,0.005,17.450,0.005,17.511,0.007,17.580,0.019]

#AllWise
wise=[16.937,0.112,16.665,0.306,999,0.0,999,0.0]

# mags    FUV, NUV,j_U,j_B,j_V,j_R,j_I,j_J   ,j_H   ,jon_K ,v_j,v_H,v_K,sds_u,sds_g,sds_r,sds_i,sds_z, W1, W2, W3, W4
EC_mags = [Galex[0],Galex[2],Johnson_opt[0],Johnson_opt[2],Johnson_opt[4],Johnson_opt[6],Johnson_opt[8],Two_mass[0],Two_mass[2],Two_mass[4],999,999,999,sdss[0],sdss[2],sdss[4],sdss[6],sdss[8],wise[0],wise[2],wise[4],wise[6]]
EC_err =  [Galex[1],Galex[3],Johnson_opt[1],Johnson_opt[3],Johnson_opt[5],Johnson_opt[7],Johnson_opt[9],Two_mass[1],Two_mass[3],Two_mass[5],0.0,0.0,0.0,sdss[1],sdss[3],sdss[5],sdss[7],sdss[9],wise[1],wise[3],wise[5],wise[7]]

#EC_mags = [999,999,999,999,999,999,999,999,999,999,999,999,999,18.271,18.271,18.454,18.683,18.85,16.732,15.919,999,999]
#EC_err =  [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.05,0.0,0.0,0.0,0.0,0.016,0.007,0.008,0.012,0.05,0.060,0.098,0.0,0.0]

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

D_wdl =[]
D_wd_err = []

for i in range(len(EC_mags)):
    if EC_err[i] == 0:
        pass
    elif phot_pivot_wav[i] > 2.14:
        pass
    else:
        mu = EC_mags[i] - mod_mags[i]
        dist = 10**(mu/5. + 1)
        dist_err = 0.2*np.log(10)*(10**(0.2*mu + 1))*EC_err[i]
        D_wdl.append(dist)
        D_wd_err.append(dist_err)

D_wd = 1./(plx/1000.)#np.median(D_wdl)#np.average(D_wdl,weights=1.0/(np.array(D_wd_err)**2.0))

mod_iso_wav = np.array([0.1528,0.2271,0.3971,0.4481,0.5423,0.6441,0.8071,1.2350,1.6620,2.1590,1.2350,1.6620,2.1590,0.3146,0.4670,0.6156,0.7471,0.8918,3.3526,4.6028,11.5608,22.0883])
mod_corr = np.array([0,0,0.0915,0.0069,0.00,0.0018,-0.0014,-0.0140,0.0060,0.0080,-0.0140,0.0060,0.0080,0.0424,-0.0023,-0.0032,-0.0160,-0.0276,0,0,0,0])
mod_zero_pts = [3631.,3631.,1649.,4060.,3723.,3168.,2459.,1594.,1024.,666.8,1594.,1024.,666.8,3631.,3631.,3631.,3631.,3631.,306.681,170.663,29.0448,8.2839]

EC_flux = (10**(-0.4*(np.array(EC_mags))))*phot_zero_pts
EC_flux_err = np.log(10)*0.4*np.array(EC_err)*EC_flux
mod_flux = (10**(-0.4*(np.array(mod_mags) + mod_corr)))*mod_zero_pts

for l in range(len(EC_flux)):
    if EC_flux_err[l]>1e-10:
        if EC_flux_err[l]/EC_flux[l] < 0.05:
            EC_flux_err[l] = EC_flux[l]*0.05

scale = []
scale_err = []
for i in range(len(EC_mags)):
    if EC_err[i] == 0:
        pass
    else:
        if phot_pivot_wav[i] > 2.16 or phot_pivot_wav[i] < 0.25:
            pass
        else:
            scale.append(EC_flux[i]/mod_flux[i])
            scale_err.append(EC_flux_err[i])

flux_scale = np.average(scale,weights=1.0/(np.array(scale_err)**2.0))
flux_scale = np.median(scale)

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
wd_mod_pts = mod_flux[phot_inds]*flux_scale
iso_wav_pts = mod_iso_wav[phot_inds]
iso_nu_pts = (clight*(10**4))/iso_wav_pts

D_wd_wdrad = (D_wd)*(3.08e18)/R_wd # Returns the distance to the wd in wd radii

def bbody_only(nu,T,R):
    return ((R/(D_wd_wdrad))**2.0)*((2.0*np.pi*hpl*(nu**3.0))/(clight**2.0))*(1.0/np.expm1(hpl*nu/(kB*T)))*(10**23.0)

bbody_scale = []
wd_bbody = bbody_only(iso_nu_pts,Teff,1.0)

for i in np.arange(2,len(wd_bbody)):
    b_scale = wd_mod_pts[i]/wd_bbody[i]
    bbody_scale.append(b_scale)
bbody_scale = np.mean(bbody_scale)

plt.figure(1)
plt.clf()
plt.xscale('log')
plt.yscale('log')
plt.errorbar(iso_wav_pts,1000.0*EC_flux_pts,yerr=1000.0*EC_flux_err_pts,fmt='s',mfc='none',color='b',label='Observed Photometry',markersize=5.)
#plt.errorbar(upper_wav_pts,1000.0*EC_flux_upper,yerr=EC_flux_upper_err,fmt='s',mfc='none',color='b',uplims=True)
#plt.plot(iso_wav_pts,1000.0*wd_mod_pts,'k*',label='WD Model Phot',markersize=5.)
plt.plot(plot_wav,bbody_only(plot_nubb,Teff,1.0)*bbody_scale*1000.0,'k--',label='WD Blackbody Model: T=%s$\,$(K) log g=%s'%(Teff,logg))

"""Now we move on to single temp IR excess fitting"""

sig_sfb = 5.67040e-5 # stefan-boltzmann in erg cm^-2 K^-4 s^-1
solar_lum = 3.826e33 # in ergs/s
solar_radius = 6.9599e10 # in cm

# Disk relevant constants
clight = 2.997e10 # cm s^-1
hpl = 6.626e-27 #erg s
kB = 1.381e-16 #erg K^-1
parsec = 3.0857e18 #c

#  Quick filter to handle any left over weird magnitudes

#fit_wavs,ec_fit_flux,ec_fit_eflux,mod_fit_flux_scaled,all_mod_wavs,all_mod_flux_scaled,scale,age,mass,R_wd_cm = ec_wdfitted_flux(ec_table,test_ind)

fit_wavs=iso_wav_pts
ec_fit_flux=EC_flux_pts
ec_fit_eflux=EC_flux_err_pts
all_mod_wavs=mod_iso_wav
mod_fit_flux_scaled=wd_mod_pts
all_mod_flux_scaled=mod_flux*flux_scale

ir_wavs = np.where(fit_wavs > 0.55)[0]
fit_nu = clight*(10.0**4.0)/fit_wavs
all_nu = clight*(10.0**4.0)/all_mod_wavs

dist = D_wd# 10.*np.sqrt(1./scale) # Gives us distance in parsecs

wd_radius = R_wd # wd radius in cm
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
T_IR = best_T
R_IR = best_R
print('Best fitting blackbody has T %s K and R %s solar' %(np.round(best_T),np.round(best_R,decimals=2)))

sort_inds = np.argsort(all_mod_wavs)
plt.figure(1)
plt.plot(plot_wav,bbody_only(plot_nubb,Teff,1.0,dist_wd_rad)*bbody_scale*1000.0 + +bbody_only(plot_nubb,best_T,best_R/(wd_radius/solar_radius),dist_wd_rad)*1000,'g-',label='WD Model + Blackbody: T=%s$\,$(K) R=%s$\,$R$_\odot$'%(best_T[0],best_R[0]))
plt.xlabel('Wavelength ($\mu$m)', fontsize=20)
plt.xlim(0.1,11.)
plt.ylabel('Flux (mJy)', fontsize=20)
plt.legend(loc=3, fontsize=14)
plt.title(obj_name, fontsize=24)
plt.savefig('plots//'+obj_name+'_sed.png',format='png')
plt.show()

roc_wd,roc_tcool,roc_teff,roc_tir,roc_tau,roc_ref,gas,mdot = np.genfromtxt('roc_table_3.txt',unpack=True,delimiter=',')
rwdlogg8 = 0.0126 #r solar
roc_rir = ((roc_teff/roc_tir)**2.0)*rwdlogg8*np.sqrt(roc_tau/100.)
has_gas = np.where(gas==1)[0]
no_gas = np.where(gas==2)[0]

cols=['#809BC8', '#FF6666', '#FFCC66', '#64C204']

plt.figure(2)
plt.clf()
plt.plot(roc_tir,roc_rir/rwdlogg8,'d',mfc='None',mec='0.7',markersize=10.0,zorder=0)
#plt.plot(roc_tir,roc_rir,'+',mfc='None',mec='0.8',markersize=10.0,zorder=0)
plt.plot(roc_tir[has_gas],roc_rir[has_gas]/rwdlogg8,'d',mfc=cols[1],mec='k',markersize=10.0,markeredgewidth=2.0,zorder=1)
plt.plot(roc_tir[no_gas],roc_rir[no_gas]/rwdlogg8,'d',mfc='0.7',mec='0.7',markersize=10.0,zorder=1)
plt.plot(T_IR[0],R_IR[0]*solar_radius/R_wd,'*',mfc='g',mec='g',markersize=20.,zorder=2)
plt.xlabel('T$_\mathrm{eff}$ (K)',size='large')
plt.ylabel('R$_{BB}$ / R$_{WD}$',size='large')
plt.xlim(0,2100)
plt.ylim(0.,41)
plt.title(obj_name)
plt.savefig('plots//'+obj_name+'_dust_region.png',format='png')

"""Finally we produce chi-sqr IR excess fit plot"""

bbody_Rarraym,bbody_Tarraym = np.meshgrid(bbody_Rarray,bbody_Tarray)
vbbody_chi_sqr = np.vectorize(bbody_chi_sqr)
chi_sqr_grid = vbbody_chi_sqr(bbody_Tarraym,bbody_Rarraym)

chi_sqr_gridt = chi_sqr_grid[:,:]

fine_bbodyT_array = np.linspace(500,5500,1000)
fine_bbodyR_array = np.linspace(0.01,1.2,1000)

ifine_bbodyT,ifine_bbodyR = np.meshgrid(fine_bbodyT_array,fine_bbodyR_array/0.0126)

chi_sqr_spl = interpolate.RectBivariateSpline(bbody_Tarray,bbody_Rarray,chi_sqr_gridt)
fine_chi_sqr_grid = chi_sqr_spl.ev(ifine_bbodyT,ifine_bbodyR)

#print('Minimum Chi-sqr value is %s' % (chi_sqr_gridt.min())
levels = [0,2.3+chi_sqr_gridt.min(),6.18+chi_sqr_gridt.min(),11.8+chi_sqr_gridt.min(),chi_sqr_gridt.max()]

chab_mass,chab_teff,chab_lum,chab_logg,chab_rad,chab_li,chab_Vmag,chab_Rmag, \
        chab_Imag,chab_Jmag,chab_Kmag,chab_Lmag,chab_Mmag = np.genfromtxt('chab_tables.txt',unpack=True)
cb_mh,cb_Y,cb_mass,cb_age,cb_teff,cb_lum,cb_rad,cb_logTc,cb_logrhoc,cb_li,cb_be = np.genfromtxt('cb_97tables.txt',unpack=True)
cb_rad = (cb_rad/6.9599)

plt.figure(3)
plt.clf()
one_sig = mpatches.Patch(color='r',label='1-$\sigma$')
two_sig = mpatches.Patch(color='g',label='2-$\sigma$')
three_sig = mpatches.Patch(color='b',label='3-$\sigma$')
plt.plot(chab_teff,chab_rad,linestyle='none',marker='s',ms=10,mfc='0.85',mec='0.85',zorder=0)
plt.plot(cb_teff,cb_rad,linestyle='none',marker='s',ms=10,mfc='0.85',mec='0.85',zorder=0)
plt.plot(roc_tir,roc_rir,'d',mfc='0.7',mec='0.7',markersize=10.0,zorder=0)
#plt.plot(roc_tir,roc_rir,'+',mfc='None',mec='0.8',markersize=10.0,zorder=0)
plt.plot(roc_tir[has_gas],roc_rir[has_gas],'d',mfc=cols[1],mec='k',markersize=10.0,markeredgewidth=2.0,zorder=1)
plt.plot(roc_tir[no_gas],roc_rir[no_gas],'d',mfc='0.7',mec='0.7',markersize=10.0,zorder=1)
cont = plt.contour(fine_bbodyT_array,fine_bbodyR_array,fine_chi_sqr_grid,levels,colors='k',linewidths=0.5,linestyles='dashed',zorder=0)
contf = plt.contourf(fine_bbodyT_array,fine_bbodyR_array,fine_chi_sqr_grid,levels,colors=['r','g','b','w'],zorder=0)
plt.xlim(0,5000)
plt.ylim(0,0.7)
plt.legend(handles=[one_sig,two_sig,three_sig],loc=5)
plt.ylabel('Bbody Radius (Solar Radii)')
plt.xlabel('Bbody Temperature (K)')
plt.title(obj_name)
plt.savefig('plots//'+obj_name+'_bbody_chi_square.png',format='png')
