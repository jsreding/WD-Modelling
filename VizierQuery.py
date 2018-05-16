import numpy as np
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
import csv

v = Vizier(columns=['_RAJ2000', '_DEJ2000', 'umag', 'gmag', 'rmag'])
l = [("RA", "Dec", "umag", "gmag", "rmag")]

with open('/home/jsreding/Documents/UNC/Research/Projects/MagWDs/SpotList.csv', 'rb') as list:
    read = csv.reader(list, delimiter='\t')
    for row in read:
        ra, dec = row[0], row[1]
        result = v.query_region(coord.SkyCoord(ra=ra, dec=dec, frame='icrs', unit=(u.hourangle, u.deg)), radius=0.1*u.deg, catalog='V/147/sdss12')
        if len(result) == 0:
            print ra, dec, "NA", "NA", "NA"
            l.append((ra, dec, "NA", "NA", "NA"))
        else:
            print ra, dec, result['V/147/sdss12']['umag'][0], result['V/147/sdss12']['e_umag'][0], result['V/147/sdss12']['gmag'][0], result['V/147/sdss12']['e_gmag'][0], result['V/147/sdss12']['rmag'][0],
            result['V/147/sdss12']['e_rmag'][0]
            l.append((ra, dec, result['V/147/sdss12']['umag'][0], result['V/147/sdss12']['gmag'][0], result['V/147/sdss12']['rmag'][0]))
# print result['V/147/sdss12']['umag'][0] - result['V/147/sdss12']['gmag'][0]
# print result['V/147/sdss12']['gmag'][0] - result['V/147/sdss12']['rmag'][0]
