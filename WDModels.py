import numpy as np
import csv
import sys
from scipy import interpolate as interp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from shapely.geometry import LineString
import time
from matplotlib import ticker
from sklearn.neighbors.kde import KernelDensity
import math

from matplotlib import rc, rcParams
rcParams['font.family'] = 'afm'
# rcParams['font.afm'] = ['Helvetica']
rc('text', usetex=True)

def find_nearest(array,value):
    d = np.shape(array)
    id = (np.abs(array-value)).argmin()
    idx, idy = id/d[1], id%d[1]
    return idx, idy, array[idx][idy]

class OOMFormatter(ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=False):
        self.oom = order
        self.fformat = fformat
        ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=False)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % ticker._mathdefault(self.format)

teff = []
logg = []
M = []
mbol = []
bc = []
u = []
g = []
r = []
i = []
z = []
U = []
jonB = []
V = []
jonR = []
I = []
J = []
H = []
K = []
age = []
with open(sys.argv[1], 'r') as f:
    # reader = csv.DictReader(f, delimiter=' ')
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        teff.append(float(row['Teff']))
        M.append(float(row['M/Mo']))
        mbol.append(float(row['Mbol']))
        bc.append(float(row['BC']))
        #Corrections from Holberg and Bergeron 2006
        u.append(float(row['u'])+0.0424)
        g.append(float(row['g'])-0.0023)
        r.append(float(row['r'])-0.0032)
        i.append(float(row['i'])-0.0160)
        z.append(float(row['z'])-0.0276)
        U.append(float(row['U'])+0.0915)
        jonB.append(float(row['B'])+0.0069)
        V.append(float(row['V'])+0.0000)
        jonR.append(float(row['R'])+0.0018)
        I.append(float(row['I'])-0.0014)
        J.append(float(row['J'])-0.0140)
        H.append(float(row['H'])+0.0060)
        K.append(float(row['K'])+0.0080)
        age.append(float(row['Age']))
# ug = np.asarray(u)-np.asarray(g)
# gr = np.asarray(g)-np.asarray(r)
# ri = np.asarray(r)-np.asarray(i)
# iz = np.asarray(i)-np.asarray(z)
# UB = np.asarray(U)-np.asarray(B)
# BV = np.asarray(B)-np.asarray(V)
# VR = np.asarray(V)-np.asarray(R)
# RI = np.asarray(R)-np.asarray(I)
# IJ = np.asarray(I)-np.asarray(J)
# JH = np.asarray(J)-np.asarray(H)
# HK = np.asarray(H)-np.asarray(K)

def colorfit(teff_in, logg_in):
    teffi = teff[:58]
    loggi = np.arange(7.0, 9.6, 0.5)
    color = [M]
    intdata = []
    for data in color:
        datai = np.array([data[:58], data[58:116], data[116:174], data[174:232], data[232:290], data[290:348]])
        bvspl_teff = interp.RectBivariateSpline(loggi, teffi, datai, kx=3, ky=3)
        intdata.append(bvspl_teff(logg_in, teff_in)[0][0])
    # teffi = teff[:48]
    # loggi = np.arange(7.0, 9.1, 0.5)
    # color = [ug, gr]
    # intdata = []
    # for data in color:
    #     datai = np.array([data[:48], data[48:96], data[96:144], data[144:192], data[192:240]])
    #     bvspl_teff = interp.RectBivariateSpline(loggi, teffi, datai, kx=3, ky=3)
    #     intdata.append(bvspl_teff(logg_in, teff_in)[0][0])
    return intdata

tf = input("teff: ")
lg = input("logg: ")
print(colorfit(tf, lg))

tr = np.linspace(1500, 120000, 1000)
lgr = np.linspace(7.0, 9.6, 1000)
U = np.zeros((1000, 1000))
G = np.zeros((1000, 1000))
R = np.zeros((1000, 1000))
I = np.zeros((1000, 1000))
Z = np.zeros((1000, 1000))
JONB = np.zeros((1000, 1000))
JONR = np.zeros((1000, 1000))
# MS = np.zeros((1000, 1000))
for l in range(len(lgr)):
    print(l)
    for t in range(len(tr)):
        U[l][t] = colorfit(tr[t], lgr[l])[0]
        G[l][t] = colorfit(tr[t], lgr[l])[1]
        R[l][t] = colorfit(tr[t], lgr[l])[2]
        I[l][t] = colorfit(tr[t], lgr[l])[3]
        Z[l][t] = colorfit(tr[t], lgr[l])[4]
        JONB[l][t] = colorfit(tr[t], lgr[l])[5]
        JONR[l][t] = colorfit(tr[t], lgr[l])[6]
        # MS[l][t] = colorfit(tr[t], lgr[l])[2]
np.save("usurf.npy", U)
np.save("gsurf.npy", G)
np.save("rsurf.npy", R)
np.save("isurf.npy", I)
np.save("zsurf.npy", Z)
np.save("JONBsurf.npy", JONB)
np.save("JONRsurf.npy", JONR)
U = np.load("usurf.npy")
G = np.load("gsurf.npy")
R = np.load("rsurf.npy")
I = np.load("isurf.npy")
Z = np.load("zsurf.npy")
JONB = np.load("JONBsurf.npy")
JONR = np.load("JONRsurf.npy")
plt.figure()
plt.contour(tr, lgr, U, [13.6], colors='c')
plt.contour(tr, lgr, JONB, [12.46], colors='b')
plt.contour(tr, lgr, G, [13.08], colors='g')
plt.contour(tr, lgr, R, [13.02], colors='y')
plt.contour(tr, lgr, I, [13.09], colors='m')
plt.contour(tr, lgr, JONR, [12.86], colors='r')
plt.contour(tr, lgr, Z, [13.15], colors='k')
plt.show()
# print("DONE")
# # # np.save("masssurf.npy", MS)
# #
# UG = np.load("ugsurf_DB.npy")
# GR = np.load("grsurf_DB.npy")
# # MS = np.load("masssurf.npy")
# # # # RI = np.load("risurf.npy")
# # surfs = [UG, GR]
# # nms = ['u-g', 'g-r']
# # for s in range(len(surfs)):
# # plt.figure()
# # plt.title("%s contours in $\log g$-$T_{eff}$ space"%(nms[s]), fontsize=18)
# # plt.xlabel("$T_{eff}$", fontsize=14)
# # plt.ylabel("$\log g$", fontsize=14)
# # plt.clabel(plt.contour(tr, lgr, MS, np.linspace(-1, 2, 20)))
# # plt.show()
# # # #
# trange = np.linspace(3500, 40000, 1000)
# lgrange = np.arange(7.0, 9.1, 0.5)
# plt.figure()
# # # ax = plt.axes(projection='3d')
# plt.title("K2 DB white dwarfs and $T_{eff}$, $\log g$ surface in $u-g$, $g-r$ space", fontsize=32)
# plt.ylabel("$u-g$", fontsize=24)
# plt.xlabel("$g-r$", fontsize=24)
# for l in lgrange:
#     ugr = []
#     grr = []
#     # rir = []
#     for t in trange:
#         ugr.append(colorfit(t, l)[0])
#         grr.append(colorfit(t, l)[1])
#         # rir.append(colorfit(t, l)[2])
#     plt.plot(grr, ugr, zorder=1)
#     # , label='$\log g = %s$'%(l)
# tcoarse = np.arange(3500, 40100, 500)
# for t in tcoarse:
#     ugr = []
#     grr = []
#     for l in lgrange:
#         ugr.append(colorfit(t, l)[0])
#         grr.append(colorfit(t, l)[1])
#         # rir.append(colorfit(t, l)[2])
#     if t%5000. == 0:
#         plt.annotate('$%s K$'%(t), xy=(grr[-1],ugr[-1]), xycoords='data', fontsize=18)
#     plt.plot(grr, ugr, color='black', zorder=1)
# plt.annotate('$log g = 7.0$', xy=(-0.35,0.55), xycoords='data', fontsize=18)
# unug = np.load('unug.npy').astype(np.float)
# ungr = np.load('ungr.npy').astype(np.float)
# unug_da = unug[:982]
# ungr_da = ungr[:982]
# unug_db = unug[982:1075]
# ungr_db = ungr[982:1075]
# unug_dc = unug[1075:1140]
# ungr_dc = ungr[1075:1140]
# unug_other = unug[1140:]
# ungr_other = ungr[1140:]
# spug = np.load('spug.npy').astype(np.float)
# spgr = np.load('spgr.npy').astype(np.float)
# spug_da = spug[:22]
# spgr_da = spgr[:22]
# spug_db = spug[22:25]
# spgr_db = spgr[22:25]
# spug_dc = spug[25:29]
# spgr_dc = spgr[25:29]
# spug_other = spug[29:]
# spgr_other = spgr[29:]
# plt.scatter(ungr_da, unug_da, label='Unspotted DA', s=20, color='c', alpha=.25, zorder=2)
# plt.scatter(ungr_db, unug_db, label='Unspotted DB', s=20, color='m', alpha=.25, zorder=2)
# plt.scatter(ungr_dc, unug_dc, label='Unspotted DC', s=20, color='y', alpha=.25, zorder=2)
# plt.scatter(ungr_other, unug_other, label='Unspotted Other', s=20, color='k', alpha=.25, zorder=2)
# plt.scatter(spgr_da, spug_da, label='Spotted DA', marker='*', s=400, color='c', edgecolors='Black', zorder=3)
# plt.scatter(spgr_db, spug_db, label='Spotted DB', marker='*', s=400, color='m', edgecolors='Black', zorder=3)
# plt.scatter(spgr_dc, spug_dc, label='Spotted DC', marker='*', s=400, color='y', edgecolors='Black', zorder=3)
# plt.scatter(spgr_other, spug_other, label='Spotted Other/Unknown', marker='*', s=400, color='k', edgecolors='Black', zorder=3)
# plt.legend(fontsize=18)
# plt.xlim(-0.6, 0.3)
# plt.ylim(-0.6, 0.8)
# plt.tick_params(axis='both', which='major', labelsize=18)
# plt.gca().invert_yaxis()
# plt.show()
#
tr = np.linspace(1500, 120000, 1000)
lgr = np.linspace(7.0, 9.5, 1000)
UG = np.load("ugsurf.npy").astype(np.float)
GR = np.load("grsurf.npy").astype(np.float)
RI = np.load("risurf.npy").astype(np.float)
IZ = np.load("izsurf.npy").astype(np.float)

# un_temps = []
# un_lg = []
# test = 'y'
# while test != 'y':
t1 = input("u-g = ")
t2 = input("g-r = ")
t3 = input("r-i = ")
t4 = input("i-z = ")
inp = [t1, t2, t3, t4]
# count = 0
# success = 0
# for u in range(len(spug)):
#     if count == 35:
#         print("***")
#         print(success)
#         print("***END OF DA***")
#         print("***")
#     elif count == 39:
#         print("***")
#         print(success)
#         print("***END OF DB***")
#         print("***")
#     elif count == 41:
#         print("***")
#         print(success)
#         print("***END OF DC***")
#         print("***")
    # inp = [spug[u], spgr[u]]
col = [UG, GR, RI, IZ]
try:
    for c in range(len(col)):
        cont = plt.contour(tr, lgr, col[c], [inp[c]])
        # for n in range(2):
        #     print(n)
        if c == 0:
            # ugt = [p[0] for p in cont.collections[0].get_paths()[0].vertices]
            # print(ugt)
            # uglg = [p[1] for p in cont.collections[0].get_paths()[0].vertices]
            path1 = LineString(cont.collections[0].get_paths()[0].vertices)
            print("ug")
        elif c == 1:
            # grt = [p[0] for p in cont.collections[0].get_paths()[0].vertices]
            # grlg = [p[1] for p in cont.collections[0].get_paths()[0].vertices]
            path2 = LineString(cont.collections[0].get_paths()[0].vertices)
            print("gr")
        elif c == 2:
            # rit = [p[0] for p in cont.collections[0].get_paths()[0].vertices]
            # rilg = [p[1] for p in cont.collections[0].get_paths()[0].vertices]
            path3 = LineString(cont.collections[0].get_paths()[0].vertices)
            print("ri")
        elif c == 3:
            # izt = [p[0] for p in cont.collections[0].get_paths()[0].vertices]
            # izlg = [p[1] for p in cont.collections[0].get_paths()[0].vertices]
            path4 = LineString(cont.collections[0].get_paths()[0].vertices)
            print("iz")
        plt.show()
    # Tfinal, LGfinal = path1.intersection(path2).x, path1.intersection(path2).y
    # print(Tfinal, LGfinal)
    # un_temps.append(Tfinal)
    # un_lg.append(LGfinal)
    # count += 1
    # success += 1
except:
    print("Outside the range of the Bergeron Models")
    # count += 1
#
# print(count)


#
# un_temps = np.load("un_teff3d.npy").astype(np.float)
# un_lgs = np.load("un_logg3d.npy").astype(np.float)
# sp_temps = np.load("sp_teff3d.npy").astype(np.float)
# sp_lgs = np.load("sp_logg3d.npy").astype(np.float)
# #
#
# sphist, spbins = np.histogram(sp_temps, [0, 7337, 8521, 9494, 10280, 11363, 12580, 13837, 15166, 16380, 18239, 20589, 23810, 28987, 55967, 120000])
# unhist, unbins = np.histogram(un_temps, [0, 7337, 8521, 9494, 10280, 11363, 12580, 13837, 15166, 16380, 18239, 20589, 23810, 28987, 55967, 120000])
# print(len(sphist))
#
# spfrac = np.zeros(len(sphist))
# unfrac = np.zeros(len(unhist))
# frac = np.zeros(len(sphist))
# for b in range(len(sphist)):
#     print(sphist[b], unhist[b])
#     frac[b] = sphist[b]/(sphist[b]+unhist[b])
#     spfrac[b] = sphist[b]/len(sp_temps)
#     unfrac[b] = unhist[b]/len(un_temps)
#
# print(frac)
# print(unhist[-1], sphist[-1])
#
# plt.figure()
# plt.rcParams["patch.force_edgecolor"] = True
# plt.minorticks_on()
# ax = plt.gca()
# plt.title("K2 White Dwarf Spot Fraction by 100 WD Temperature Bin", fontsize=32)
# plt.ylabel("Spotted (%)", fontsize=24)
# plt.xlabel("$T_{eff}$ (K)", fontsize=24)
# # plt.xscale('log')
# # plt.bar([2150, 7179, 8100, 9149, 9974, 10842, 11852, 13080, 14244, 15304, 16438, 18240, 20224, 22997, 27292, 37222], spfrac*100, align='edge', width=[7179-2151, 8100-7179, 9149-8100, 9974-9149, 10842-9974, 11852-10842, 13080-11852, 14244-13080, 15304-14244, 16438-15304, 18240-16438, 20224-18240, 22997-20224, 27292-22997, 37222-27292, 40000-37222], alpha=0.33, label='Fraction of total Spotted')
# # plt.bar([2150, 7179, 8100, 9149, 9974, 10842, 11852, 13080, 14244, 15304, 16438, 18240, 20224, 22997, 27292, 37222], unfrac*100, align='edge', width=[7179-2151, 8100-7179, 9149-8100, 9974-9149, 10842-9974, 11852-10842, 13080-11852, 14244-13080, 15304-14244, 16438-15304, 18240-16438, 20224-18240, 22997-20224, 27292-22997, 37222-27292, 40000-37222], alpha=0.33, label='Fraction of total Unspotted')
# barlist = plt.bar([3648, 7337, 8521, 9494, 10280, 11363, 12580, 13837, 15166, 16380, 18239, 20589, 23810, 28987, 55967], frac*100, align='edge', width=[7337-3648, 8521-7337, 9494-8521, 10280-9494, 11363-10280, 12580-11363, 13837-12580, 15166-13837, 16380-15166, 18239-16380, 20589-18239, 23810-20589, 28987-23810, 55967-28987, 120000-55967], color='c', alpha=0.5)
# barlist[-1].set_color('m')
# # [8100-2151, 9974-8100, 11852-9974, 14243-11852, 16473-14244, 20223-16438, 27292-20224, 120000-27292]
# # [7179-2151, 8100-7179, 9149-8100, 9974-9149, 10842-9974, 11852-10842, 13080-11852, 14244-13080, 15304-14244, 16438-15304, 18240-16438, 20224-18240, 22997-20224, 27292-22997, 37222-27292, 120000-37222]
# plt.tick_params(axis='x', which='minor')
# plt.tick_params(axis='both', which='major', labelsize=18)
# # plt.xticks([2151, 7179, 8100, 9149, 9974, 10842, 11852, 13080, 14244, 15304, 16438, 18240, 20224, 22997, 27292, 37222])
# plt.xscale('log')
# plt.xlim(4000, 100000)
# plt.grid(axis='x', which='both', linestyle='--')
# plt.tick_params(axis='x', which='major', labelsize=18, pad=10)
# plt.tick_params(axis='x', which='minor')
# ax.xaxis.set_major_formatter(OOMFormatter(3, "%0.0f"))
# ax.xaxis.set_minor_formatter(OOMFormatter(3, "%0.0f"))
# ax.ticklabel_format(axis='x', style='sci', scilimits=(3,3))
# plt.show()
#
# # fig  = plt.figure()
# # ax1 = fig.add_subplot(211)
# # ax1.minorticks_on()
# # ax1.set_title("K2 white dwarf mass histogram", fontsize=18)
# # ax1.set_ylabel("Density (# of stars)", fontsize=14)
# # ax1.set_xlabel("Mass ($M_{\odot}$)", fontsize=14)
# # ax1.hist(un_mass, bins=np.arange(0, 1.45, 0.05), label='Unspotted')
# # ax1.set_xlim(0, 1.4)
# # plt.legend()
# # ax2 = fig.add_subplot(212)
# # ax2.minorticks_on()
# # ax2.set_ylabel("Density (# of stars)", fontsize=14)
# # ax2.set_xlabel("Mass ($M_{\odot}$)", fontsize=14)
# # ax2.hist(sp_mass, bins=np.arange(0, 1.45, 0.05), color='Orange', label='Spotted')
# # ax2.set_xlim(0, 1.4)
# # plt.legend()
# # plt.show()
#
# # xplot = np.linspace(0, 120000, 1000)[:, np.newaxis]
# # un_kde = KernelDensity(kernel='gaussian', bandwidth=1000.0).fit(un_temps.reshape(-1, 1))
# # sp_kde = KernelDensity(kernel='gaussian', bandwidth=1000.0).fit(sp_temps.reshape(-1, 1))
# # un_log_dens = un_kde.score_samples(xplot)
# # sp_log_dens = sp_kde.score_samples(xplot)
#
# # plt.figure()
# # plt.minorticks_on()
# # plt.title("K2 White Dwarf Density vs $T_{eff}$", fontsize=32)
# # plt.ylabel("Density", fontsize=24)
# # plt.xlabel("Temperature (K)", fontsize=24)
# # plt.plot(xplot[:, 0], np.exp(un_log_dens), color='c', alpha=0.5, label='Unspotted')
# # plt.fill_between(xplot[:, 0], np.exp(un_log_dens), color='c', alpha=0.5)
# # plt.plot(xplot[:, 0], np.exp(sp_log_dens), color='m', alpha=0.5, label='Spotted')
# # plt.fill_between(xplot[:, 0], np.exp(sp_log_dens), color='m', alpha=0.5)
# # plt.xlim(0, 120000)
# # plt.tick_params(axis='both', which='major', labelsize=18)
# # plt.legend(fontsize=18)
# # plt.show()
#
# for u in range(len(un_lgs)):
#     if un_lgs[u] > 9.5:
#         un_lgs[u] = 9.5
#
# for s in range(len(sp_lgs)):
#     if sp_lgs[s] > 9.5:
#         sp_lgs[s] = 9.5
#
# plt.figure()
# plt.minorticks_on()
# ax = plt.gca()
# plt.title(r'K2 White Dwarfs, $\log g, T_\mathrm{eff}$ Estimates', fontsize=32)
# plt.xlim(4000, 100000)
# plt.ylim(6.5, 9.6)
# plt.xlabel(r'$T_\mathrm{eff}$ (kK)', fontsize=24)
# plt.ylabel(r'$\log g$', fontsize=24)
# plt.xscale('log')
# ax.axvspan(10500, 12500, alpha=0.2, zorder=0, color='c')
# ax.axvspan(20000, 32000, alpha=0.2, zorder=0, color='m')
# plt.scatter(un_temps[:887], un_lgs[:887], marker='o', color='c', zorder=1, alpha=0.2, label='Unspotted, DA')
# plt.scatter(un_temps[887:967], un_lgs[887:967], marker='o', color='m', zorder=1, alpha=0.2, label='Unspotted, DB')
# plt.scatter(un_temps[967:1016], un_lgs[967:1016], marker='o', color='y', zorder=1, alpha=0.2, label='Unspotted, DC')
# plt.scatter(un_temps[1016:], un_lgs[1016:], marker='o', color='k', zorder=1, alpha=0.2, label='Unspotted, Other/Unk.')
# plt.scatter(sp_temps[:32], sp_lgs[:32], s=400, marker='*', color="c", edgecolors='Black', zorder=2, label="Spotted, DA")
# plt.scatter(sp_temps[32:35], sp_lgs[32:35], s=400, marker='*', color="m", edgecolors='Black', zorder=2, label="Spotted, DB")
# plt.scatter(sp_temps[35:37], sp_lgs[35:37], s=400, marker='*', color="y", edgecolors='Black', zorder=2, label="Spotted, DC")
# plt.scatter(sp_temps[37:], sp_lgs[37:], s=400, marker='*', color="k", edgecolors='Black', zorder=2, label="Spotted, Other/Unk.")
# plt.scatter(un_temps[827:829], un_lgs[827:829], marker='o', facecolors='None', edgecolors='Black', s=400, label='Magnetic ($B \geq 2MG$)')
# plt.scatter(un_temps[831:834], un_lgs[831:834], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_temps[835], un_lgs[835], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_temps[839:843], un_lgs[839:843], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_temps[846], un_lgs[846], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_temps[848], un_lgs[848], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_temps[850:853], un_lgs[850:853], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_temps[854:866], un_lgs[854:866], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_temps[868], un_lgs[868], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_temps[870:872], un_lgs[870:872], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_temps[956:964], un_lgs[956:964], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_temps[1038], un_lgs[1038], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(sp_temps[21:24], sp_lgs[21:24], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(sp_temps[27:30], sp_lgs[27:30], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(sp_temps[34], sp_lgs[34], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.tick_params(axis='both', which='major', labelsize=22)
# plt.tick_params(axis='x', which='minor', labelsize=20)
# ax.xaxis.set_major_formatter(OOMFormatter(3, "%0.0f"))
# ax.xaxis.set_minor_formatter(OOMFormatter(3, "%0.0f"))
# ax.ticklabel_format(axis='x', style='plain')#, scilimits=(0,0))
# # plt.xticks(np.arange(4000, 10000, 1000), ('4', '5', '6', '7', '8', '9', '10'))
# # plt.xticks(np.arange(10000, 100000, 10000), ('20', '30', '40', '50', '60', '70', '80', '90', '100'))
# plt.grid()
# plt.grid(axis='x', which='minor', linestyle='--')
# plt.legend(fontsize=18)
# plt.gca().invert_yaxis()
# plt.show()
#
# un_G = np.load('un_GaiaG.npy').astype(float)
# un_BR = np.load('un_GaiaBR.npy').astype(float)
# sp_G = np.load('sp_GaiaG.npy').astype(float)
# sp_BR = np.load('sp_GaiaBR.npy').astype(float)
#
# plt.figure()
# plt.title("Gaia CMD of K2 White Dwarfs", fontsize=24)
# plt.xlabel("$G_\mathrm{BP} - G_\mathrm{RP}$", fontsize=20)
# plt.xlim(-0.75, 1.25)
# plt.ylabel("$G_\mathrm{abs}$", fontsize=20)
# plt.ylim(5, 15)
# plt.scatter(un_BR[:945], un_G[:945], marker='o', color='c', zorder=1, alpha=0.2, label='DA')
# plt.scatter(un_BR[945:1037], un_G[945:1037], marker='o', color='m', zorder=1, alpha=0.2, label='DB')
# plt.scatter(un_BR[1037:1089], un_G[1037:1089], marker='o', color='y', zorder=1, alpha=0.2, label='DC')
# plt.scatter(un_BR[1089:], un_G[1089:], marker='o', color='k', zorder=0, alpha=0.2, label='Other/Unk.')
# # plt.scatter(sp_BR[:29], sp_G[:29], s=400, marker='*', color='c', edgecolors='k', zorder=2, label='Spotted, DA')
# # plt.scatter(sp_BR[29:33], sp_G[29:33], s=400, marker='*', color='m', edgecolors='k', zorder=2, label='Spotted, DB')
# plt.scatter(sp_BR[35], sp_G[35], s=400, marker='*', color='k', edgecolors='k', zorder=2, label='J1252')
# # plt.scatter(sp_BR[36:], sp_G[36:], s=400, marker='*', color='k', edgecolors='k', zorder=2, label='Spotted, Other/Unk.')
# plt.scatter(un_BR[887], un_G[887], marker='o', facecolors='None', edgecolors='Black', s=400, label='Magnetic ($B \geq 2MG$)')
# plt.scatter(un_BR[889:892], un_G[889:892], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_BR[897:901], un_G[897:901], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_BR[904], un_G[904], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_BR[907:910], un_G[907:910], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_BR[912:921], un_G[912:921], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_BR[923], un_G[923], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_BR[925:928], un_G[925:928], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_BR[930:932], un_G[930:932], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_BR[1023:1031], un_G[1023:1031], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_BR[1036], un_G[1036], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(un_BR[1113], un_G[1113], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(sp_BR[19], sp_G[19], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(sp_BR[21:23], sp_G[21:23], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(sp_BR[26], sp_G[26], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.scatter(sp_BR[32], sp_G[32], marker='o', facecolors='None', edgecolors='Black', s=400)
# plt.grid()
# plt.legend(fontsize=14)
# plt.gca().invert_yaxis()
# plt.tick_params(axis='both', which='major', labelsize=18)
# plt.show()
#
# un = open("/home/jsreding/Documents/UNC/Research/Projects/SpottedWDs/UnspottedTeffLoggSpec", "r")
# unlines = un.readlines()
# un_temps = [float(l.strip().split('\t')[0]) for l in unlines]
# un_lgs = [float(l.strip().split('\t')[1]) for l in unlines]
# sp = open("/home/jsreding/Documents/UNC/Research/Projects/SpottedWDs/SpottedTeffLoggSpec", "r")
# splines = sp.readlines()
# sp_temps = [float(l.strip().split('\t')[0]) for l in splines]
# sp_lgs = [float(l.strip().split('\t')[1]) for l in splines]
#
# un_mass = np.zeros(len(un_temps))
# sp_mass = np.zeros(len(sp_temps))
# for i in range(len(un_mass)):
#     un_mass[i] = colorfit(un_temps[i], un_lgs[i])[-1]
# for i in range(len(sp_mass)):
#     sp_mass[i] = colorfit(sp_temps[i], sp_lgs[i])[-1]
#
# plt.figure()
# plt.hist(un_mass, bins=20, alpha=0.5, color='c')
# plt.hist(sp_mass, bins=20, alpha=0.5, color='m')
# plt.show()
