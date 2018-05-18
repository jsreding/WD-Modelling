import numpy as np
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
import csv

v = Vizier(columns=['_RAJ2000', '_DEJ2000', 'umag', 'e_umag', 'gmag', 'e_gmag', 'rmag', 'e_rmag'])
l = [('RA', 'Dec', 'umag', 'e_umag', 'gmag', 'e_gmag', 'rmag', 'e_rmag')]

with open('/home/jsreding/Documents/UNC/Research/Projects/MagWDs/SpotList.csv', 'rb') as list:
    read = csv.reader(list, delimiter=' ')
    for row in read:
        ra, dec = row[0], row[1]
        result = v.query_region(coord.SkyCoord(ra=ra, dec=dec, frame='icrs', unit=(u.hourangle, u.deg)), radius=0.001*u.deg, catalog='V/139/sdss9')
        if len(result) == 0:
            print ra, dec, "NA", "NA", "NA", "NA", "NA", "NA"
            l.append((ra, dec, "NA", "NA", "NA", "NA", "NA", "NA"))
        else:
            print ra, dec, result['V/139/sdss9']['umag'][0], result['V/139/sdss9']['e_umag'][0], result['V/139/sdss9']['gmag'][0], result['V/139/sdss9']['e_gmag'][0], result['V/139/sdss9']['rmag'][0], result['V/139/sdss9']['e_rmag'][0]
            l.append((ra, dec, result['V/139/sdss9']['umag'][0], result['V/139/sdss9']['e_umag'][0], result['V/139/sdss9']['gmag'][0], result['V/139/sdss9']['e_gmag'][0], result['V/139/sdss9']['rmag'][0], result['V/139/sdss9']['e_rmag'][0]))
    # print zip(*l)
# print result['V/139/sdss9']['umag'][0] - result['V/139/sdss9']['gmag'][0]
# print result['V/139/sdss9']['gmag'][0] - result['V/139/sdss9']['rmag'][0]
