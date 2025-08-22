from astropy.io import ascii

CENTRAL_WAVE_URL = 'https://github.com/splus-collab/splus_filters/blob/master/central_wavelengths.csv'

try:
    __filters_table__ = ascii.read(CENTRAL_WAVE_URL)
except:
    __filters_table__ = None