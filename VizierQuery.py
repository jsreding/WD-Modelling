import numpy as np
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
import csv

v = Vizier(columns=['RA_ICRS', 'DE_ICRS', 'umag', 'e_umag', 'gmag', 'e_gmag', 'rmag', 'e_rmag'])
# l = [('RA_ICRS', 'DE_ICRS', 'Plx', 'e_Plx', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Gmag', 'e_Gmag', 'BPmag', 'e_BPmag', 'RPmag', 'e_RPmag')]

with open('/home/jsreding/Documents/UNC/Research/Projects/MagWDs/UnspottedList.csv', 'r') as list:
    read = csv.reader(list, delimiter='\t')
    i = 1
    for row in read:
        i += 1
        ra, dec = row[0], row[1]
        result = v.query_region(coord.SkyCoord(ra=ra, dec=dec, frame='icrs', unit=(u.hourangle, u.deg)), radius=0.001*u.deg, catalog='V/139/sdss9')
        if len(result) == 0:
            print(ra, dec, "NA")
        else:
            print(ra, dec, result['V/139/sdss9']['umag'][0], result['V/139/sdss9']['e_umag'][0], result['V/139/sdss9']['gmag'][0], result['V/139/sdss9']['e_gmag'][0], result['V/139/sdss9']['rmag'][0], result['V/139/sdss9']['e_rmag'][0])
            # l.append((ra, dec, result['V/139/sdss9']['Plx'][0], result['V/139/sdss9']['e_Plx'][0], result['V/139/sdss9']['pmRA'][0], result['V/139/sdss9']['e_pmRA'][0], result['V/139/sdss9']['pmDE'][0], result['V/139/sdss9']['e_pmDE'][0], result['V/139/sdss9']['Gmag'][0], result['V/139/sdss9']['e_Gmag'][0], result['V/139/sdss9']['BPmag'][0], result['V/139/sdss9']['e_BPmag'][0], result['V/139/sdss9']['RPmag'][0], result['V/139/sdss9']['e_RPmag'][0]))
# for x in l:
#     print(x)
# print result['umag'][0] - result['gmag'][0]
# print result['gmag'][0] - result['rmag'][0]
