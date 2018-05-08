import numpy as np
import csv
import sys
from scipy import interpolate as interp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

def find_nearest(array,value):
    d = np.shape(array)
    id = (np.abs(array-value)).argmin()
    idx, idy = id/d[1], id%d[1]
    return idx, idy, array[idx][idy]

teff = []
logg = []
M = []
mbol = []
bc = []
u = []
g = []
r = []
age = []
ug = []
gr = []
ri = []
i = 0
with open(sys.argv[1], 'rb') as f:
    reader = csv.DictReader(f, delimiter=' ')
    for row in reader:
        teff.append(float(row['Teff']))
        M.append(float(row['M/Mo']))
        mbol.append(float(row['Mbol']))
        bc.append(float(row['BC']))
        u.append(float(row['u']))
        g.append(float(row['g']))
        r.append(float(row['r']))
        age.append(float(row['Age']))
        ug.append((float(row['u'])+0.0424) - (float(row['g'])-0.0023))
        gr.append((float(row['g'])-0.0023) - (float(row['r'])-0.0032))
        ri.append((float(row['r'])-0.0032) - (float(row['i'])-0.016))

def colorfit(teff_in, logg_in):
    teffi = teff[:58]
    loggi = np.arange(7.0, 9.6, 0.5)
    color = [ug, gr, ri]
    intdata = []
    for data in color:
        datai = np.array([data[:58], data[58:116], data[116:174], data[174:232], data[232:290], data[290:348]])
        bvspl_teff = interp.RectBivariateSpline(loggi, teffi, datai, kx=1, ky=1)
        intdata.append(bvspl_teff(logg_in, teff_in)[0][0])
    return intdata

tr = np.linspace(1500, 120000, 1000)
lgr = np.linspace(7.0, 9.5, 1000)
# UG = np.zeros((1000, 1000))
# GR = np.zeros((1000, 1000))
# RI = np.zeros((1000, 1000))
# for l in range(len(lgr)):
#     print l
#     for t in range(len(tr)):
#         UG[l][t] = colorfit(tr[t], lgr[l])[0]
#         GR[l][t] = colorfit(tr[t], lgr[l])[1]
#         RI[l][t] = colorfit(tr[t], lgr[l])[2]
# np.save("ugsurf", UG)
# np.save("grsurf", GR)
# np.save("risurf", RI)

# UG = np.load("ugsurf.npy")
# GR = np.load("grsurf.npy")
# RI = np.load("risurf.npy")
# surfs = [UG, GR, RI]
# nms = ['u-g', 'g-r', 'r-i']
# for s in range(len(surfs)):
#     plt.figure()
#     plt.title("%s contours in $\log g$-$T_{eff}$ space"%(nms[s]), fontsize=18)
#     plt.xlabel("$T_{eff}$", fontsize=14)
#     plt.ylabel("$\log g$", fontsize=14)
#     plt.clabel(plt.contour(tr, lgr, surfs[s], np.linspace(-1, 2, 100)))
#     plt.show()

# trange = np.linspace(1500, 120000, 1000)
# lgrange = np.arange(7.0, 9.6, 0.5)
# plt.figure()
# # ax = plt.axes(projection='3d')
# plt.title("$T_{eff}$, $\log g$ surface in $u-g$, $g-r$ parameter space", fontsize=18)
# plt.ylabel("$u-g$", fontsize=14)
# plt.xlabel("$g-r$", fontsize=14)
# for l in lgrange:
#     ugr = []
#     grr = []
#     rir = []
#     for t in trange:
#         ugr.append(colorfit(t, l)[0])
#         grr.append(colorfit(t, l)[1])
#         rir.append(colorfit(t, l)[2])
#     plt.plot(grr, ugr, label='$\log g = %s$'%(l))
# tcoarse = np.arange(1500, 120100, 500)
# for t in tcoarse:
#     ugr = []
#     grr = []
#     for l in lgrange:
#         ugr.append(colorfit(t, l)[0])
#         grr.append(colorfit(t, l)[1])
#         rir.append(colorfit(t, l)[2])
#     if t%5000. == 0:
#         plt.annotate('$%s K$'%(t), xy=(grr[-1],ugr[-1]), xycoords='data', fontsize=14)
#     plt.plot(grr, ugr, color='black')
# plt.legend()
# plt.xlim(-0.6, 0.35)
# plt.ylim(-0.6, 1.0)
# plt.gca().invert_yaxis()
# plt.show()

tr = np.linspace(1500, 120000, 1000)
lgr = np.linspace(7.0, 9.5, 1000)
UG = np.load("ugsurf.npy")
GR = np.load("grsurf.npy")
RI = np.load("risurf.npy")
test = 'n'
while test != 'y':
    t1 = raw_input("u-g = ")
    t2 = raw_input("g-r = ")
    # t3 = raw_input("r-i = ")
    inp = [t1, t2]
    col = [UG, GR]
    # inp = [t1, t2, t3]
    # col = [UG, GR, RI]
    from shapely.geometry import LineString
    try:
        for c in range(len(col)):
            cont = plt.contour(tr, lgr, col[c], [inp[c]]);
            if c == 0:
                ugt = [p[0] for p in cont.collections[0].get_paths()[0].vertices]
                uglg = [p[1] for p in cont.collections[0].get_paths()[0].vertices]
                path1 = LineString(cont.collections[0].get_paths()[0].vertices)
            elif c == 1:
                grt = [p[0] for p in cont.collections[0].get_paths()[0].vertices]
                grlg = [p[1] for p in cont.collections[0].get_paths()[0].vertices]
                path2 = LineString(cont.collections[0].get_paths()[0].vertices)
            # elif c == 2:
            #     rit = [p[0] for p in cont.collections[0].get_paths()[0].vertices]
            #     rilg = [p[1] for p in cont.collections[0].get_paths()[0].vertices]
            #     path3 = LineString(cont.collections[0].get_paths()[0].vertices)
    except:
        print "Outside the range of the Bergeron Models"
    # Tfinal = []
    # LGfinal = []
    # for p in path1.intersection(path2):
    #     Tfinal.append(p.x)
    #     LGfinal.append(p.y)
    Tfinal, LGfinal = path1.intersection(path2).x, path1.intersection(path2).y
    print Tfinal, LGfinal

    # plt.figure()
    # plt.title("Color Contour Intersections", fontsize=20)
    # plt.xlabel("$T_{eff}$", fontsize=14)
    # plt.ylabel("$\log g$", fontsize=14)
    # plt.plot(ugt, uglg, label='$u-g=%s$'%(t1))
    # plt.plot(grt, grlg, label='$g-r=%s$'%(t2))
    # # plt.plot(rit, rilg, label='$r-i=%s$'%(t3))
    # plt.scatter(Tfinal, LGfinal, label='$T_{eff}=%s, \log g=%s$'%(Tfinal, LGfinal))
    # plt.legend(fontsize=14)
    # plt.show()
    test = raw_input("Done? ")

# dahteff = [14078., 15254., 17326., 20335., 11140., 16750., 11351., 22775., 19134., 17904., 17936., 22490., 25636., 8753., 6886., 16175., 18764., 18529., 15613., 15980., 10604., 18090., 22510., 12284., 10182., 19142., 18694., 9871., 7525., 20160., 13458., 12275., 14592., 16700., 28315., 47547., 66229., 28125., 21696., 59983., 61801., 10620.]
# dahlogg = [8.45, 9.02, 9.02, 9.0, 8.738, 9.12, 8.33, 9.05, 7.81, 7.79, 8.21, 8.69, 7.44, 9.06, 9.0, 8.37, 9.16, 8.09, 8.15, 8.16, 7.97, 7.91, 10.0, 8.4, 6.84, 7.86, 8.07, 8.04, 9.02, 9.0, 9.0, 9.0, 5.65, 8.32, 10.0, 8.84, 7.58, 9.22, 8.08, 8.0, 6.799, 9.189]
# ugr = []
# grr = []
# for i in range(len(dahteff)):
#     ugr.append(colorfit(dahteff[i], dahlogg[i])[0])
#     grr.append(colorfit(dahteff[i], dahlogg[i])[1])
# plt.scatter(grr, ugr, label='DAH', marker='o', zorder=2)
#
# dbahteff = [18254., 19992., 17772., 18421., 16005., 23629., 15413.]
# dbahlogg = [8.07, 8.02, 8.17, 7.85, 8.13, 8.32, 8.02]
# ugr = []
# grr = []
# for i in range(len(dbahteff)):
#     ugr.append(colorfit(dbahteff[i], dbahlogg[i])[0])
#     grr.append(colorfit(dbahteff[i], dbahlogg[i])[1])
# plt.scatter(grr, ugr, label='DBAH', marker='^', zorder=2)
#
# dbhteff = [34521., 40000.]
# dbhlogg = [8.44, 9.61]
# ugr = []
# grr = []
# for i in range(len(dbhteff)):
#     ugr.append(colorfit(dbhteff[i], dbhlogg[i])[0])
#     grr.append(colorfit(dbhteff[i], dbhlogg[i])[1])
# plt.scatter(grr, ugr, label='DBH', marker='s', zorder=2)
#
# dbzhmteff = [16503.]
# dbzhmlogg = [8.33]
# ugr = []
# grr = []
# for i in range(len(dbzhmteff)):
#     ugr.append(colorfit(dbzhmteff[i], dbzhmlogg[i])[0])
#     grr.append(colorfit(dbzhmteff[i], dbzhmlogg[i])[1])
# plt.scatter(grr, ugr, label='DBZHM', marker='x', zorder=2)
#
# dqhteff = [11700.]
# dqhlogg = [8.5]
# ugr = []
# grr = []
# for i in range(len(dqhteff)):
#     ugr.append(colorfit(dqhteff[i], dqhlogg[i])[0])
#     grr.append(colorfit(dqhteff[i], dqhlogg[i])[1])
# plt.scatter(grr, ugr, label='DQH', marker='*', zorder=2)
# plt.xlim(-0.6, 0.5)
# plt.ylim(1.0, -0.6)
# plt.legend()
# plt.show()
