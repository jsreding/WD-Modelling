from astropy import time, coordinates, units as u
import numpy as np

names = ["EPIC 201647223", "EPIC 248433650", "EPIC 201741620", "EPIC 248717409"]
coord = [("11:15:27.321", "02:46:22.02"), ("10:34:48.927", "00:52:01.43"), ("11:41:32.984", "04:20:28.89"), ("10:16:42.003", "08:44:00.04")]
bjd = [(2456850.231986, 0.000171), (2457945.46526, 0.00035), (2456850.223495, 0.000095), (2457946.2467, 0.0019)]
per = [(0.121426608, 0.00000090), (0.0915579, 0.0000014), (0.062289147, 0.00000026), (0.261939, 0.000022)]

offset = np.zeros(len(bjd))
for n in range(len(bjd)):
    offset[n] = int((time.Time.now().jd - bjd[n][0])/per[n][0])

name = raw_input('Which object? [2016, 2484, 2017, 2487] ')
if name == "2016":
    x = 0
elif name == "2484":
    x = 1
elif name == "2017":
    x = 2
elif name == "2487":
    x = 3

tar = []
for n in range(20):
    tar.append((bjd[x][0]+(offset[x]+n)*per[x][0], np.sqrt(bjd[x][1]**2+(offset[x]+n)*per[x][1]**2)))

obj = coordinates.SkyCoord(coord[x][0], coord[x][1], unit=(u.hourangle, u.deg))
mcdonald = coordinates.EarthLocation.of_site('mcdonald')
times = time.Time(tar, format='jd', scale='tdb', location=mcdonald)
barycorr = times.light_travel_time(obj)

print names[x]
for i in range(len(tar)):
    h = int(tar[i][1]*24)
    m = int((tar[i][1]*24-h)*60)
    s = ((tar[i][1]*24-h)*60-m)*60
    print "Center:", (times[i][0]-barycorr[i][0]).utc.datetime, "+/-", h, ':', m, ':', s
