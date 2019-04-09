# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 12:10:40 2015

@author: Erik
"""

import numpy as np
import scipy.interpolate as interp

def bergeron_DA_phot(logg_in,Teff_in):
    Teff,logg,mass,mbol,bol_corr,jon_U,jon_B,jon_V,jon_R,jon_I,jon_J,jon_H,jon_K, \
    sdss_u,sdss_g,sdss_r,sdss_i,sdss_z,sdss_y,col_by,col_ub,col_vy,W1,W2,W3,W4 \
    ,age = np.genfromtxt('bergeron_Table_DA_WISE.txt',skip_header=4,unpack=True)
    
    Teff_gal,logg_gal,mass_gal,mbol_gal,jon_V_gal,FUV,NUV,age_gal = np.genfromtxt('bergeron_Table_DA_GALEX.txt',skip_header=2,unpack=True)
        
    all_data = [Teff,logg,mass,mbol,bol_corr,age,FUV,NUV,jon_U,jon_B,jon_V,jon_R,jon_I,jon_J,jon_H,jon_K,jon_J,jon_H,jon_K,sdss_u,sdss_g, \
    sdss_r,sdss_i,sdss_z,W1,W2,W3,W4]
    
    Teffi = Teff[0:60]
    loggi = np.arange(7.0,9.6,0.5)
    
    interp_data = []
    
    for data in all_data:
        datai = np.array([data[0:60],data[60:120],data[120:180],data[180:240],data[240:300],data[300:360]])  
        biv_spline = interp.RectBivariateSpline(loggi,Teffi,datai,kx=3,ky=3)
        interp_data.append(biv_spline(logg_in,Teff_in)[0][0])        
    return np.array(interp_data) 

def bergeron_DB_phot(logg_in,Teff_in):
    Teff,logg,mass,mbol,bol_corr,jon_U,jon_B,jon_V,jon_R,jon_I,jon_J,jon_H,jon_K, \
    sdss_u,sdss_g,sdss_r,sdss_i,sdss_z,sdss_y,col_by,col_ub,col_vy,W1,W2,W3,W4 \
    ,age = np.genfromtxt('bergeron_Table_DB_WISE.txt',skip_header=2,unpack=True)
    
    Teff_gal,logg_gal,mass_gal,mbol_gal,jon_V_gal,FUV,NUV,age_gal = np.genfromtxt('bergeron_Table_DB_GALEX.txt',skip_header=2,unpack=True)
        
    all_data = [Teff,logg,mass,mbol,bol_corr,age,FUV,NUV,jon_U,jon_B,jon_V,jon_R,jon_I,jon_J,jon_H,jon_K,jon_J,jon_H,jon_K,sdss_u,sdss_g, \
    sdss_r,sdss_i,sdss_z,W1,W2,W3,W4]
    
    Teffi = Teff[0:48]
    loggi = np.arange(7.0,9.1,0.5)
    
    interp_data = []
    
    for data in all_data:
        datai = np.array([data[0:48],data[48:96],data[96:144],data[144:192],data[192:240]])  
        biv_spline = interp.RectBivariateSpline(loggi,Teffi,datai,kx=3,ky=3)
        interp_data.append(biv_spline(logg_in,Teff_in)[0][0])        
    return np.array(interp_data) 
        
        
