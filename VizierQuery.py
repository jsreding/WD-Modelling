import numpy as np
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u

v = Vizier(columns=['_RAJ2000', '_DEJ2000', 'umag', 'gmag', 'rmag'])
result = v.query_region(coord.SkyCoord(ra='1:12:43.2', dec='+20:52:16.53', frame='icrs', unit=(u.hourangle, u.deg)), radius=0.1*u.deg, catalog='V/147/sdss12')
print result['V/147/sdss12']['umag'][0], result['V/147/sdss12']['gmag'][0], result['V/147/sdss12']['rmag'][0]
print result['V/147/sdss12']['umag'][0] - result['V/147/sdss12']['gmag'][0]
print result['V/147/sdss12']['gmag'][0] - result['V/147/sdss12']['rmag'][0]
