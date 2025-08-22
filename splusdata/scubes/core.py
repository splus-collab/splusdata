import sys
import numpy as np
import astropy.units as u
from tqdm.auto import tqdm
from astropy.io import fits
from astropy.wcs import WCS
from os import makedirs, getcwd
from astropy.table import Table
from splusdata.core import Core
import astropy.constants as const
from os.path import join, exists, isfile

from splusdata.scripts.args import create_parser
from splusdata.features.io import print_level, convert_coord_to_degrees
from splusdata.vars import BANDS, BANDWAVEINFO, get_band_info

__script_author__ = 'Eduardo A. D. Lacerda <dhubax@gmail.com>'
__script_version__ = '0.1.idr6-beta'

SPLUS_MOTD_TOP = '┌─┐   ┌─┐┬ ┬┌┐ ┌─┐┌─┐ '
SPLUS_MOTD_MID = '└─┐───│  │ │├┴┐├┤ └─┐ '
SPLUS_MOTD_BOT = '└─┘   └─┘└─┘└─┘└─┘└─┘ '
SPLUS_MOTD_SEP = '----------------------'

SCUBES_PROG_DESC = f'''
{SPLUS_MOTD_TOP} | Create S-PLUS galaxies data cubes, a.k.a. S-CUBES. 
{SPLUS_MOTD_MID} | S-CUBES is an organized FITS file with data, errors, 
{SPLUS_MOTD_BOT} | mask and metadata about some galaxy present on any 
{SPLUS_MOTD_SEP} + S-PLUS observed tile. Any problem contact:

   {__script_author__}

The input values of RA and DEC will be converted to degrees using the 
scubes.utilities.io.convert_coord_to_degrees(). All scripts with RA 
and DEC inputs parse angles in two different units:

- **hourangle**: using *hms* divisors; Ex: *10h37m2.5s*
- **degrees**: using *:* or *dms*  divisors; Ex: *10:37:2.5* or *10d37m2.5s*

Note that *10h37m2.5s* is a totally different angle from *10:37:2.5* 
(*159.26 deg* and *10.62 deg* respectively).

'''

SCUBES_ARGS = {
    # optional arguments
    'redo': ['r', dict(action='store_true', default=False, help='Enable redo mode to overwrite final cubes.')],
    'force': ['f', dict(action='store_true', default=False, help='Force overwrite of existing files.')],
    'size': ['l', dict(default=500, type=int, help='Size of the cube in pixels. If size is a odd number, the program will choose the closest even integer.')],
    'workdir': ['w', dict(default=getcwd(), help='Working directory.')],
    'verbose': ['v', dict(action='count', default=0, help='Verbosity level.')],
    'username': ['U', dict(default=None, help='S-PLUS Cloud username.')],
    'password': ['P', dict(default=None, help='S-PLUS Cloud password.')],
     #'data_release': ['d', dict(default='dr4', type=str, help='Select S-PLUS Data Release')],

    # positional arguments
    'field': ['pos', dict(metavar='SPLUS_TILE', help='Name of the S-PLUS field')],
    'ra': ['pos', dict(metavar='RA', help="Object's right ascension")],
    'dec': ['pos', dict(metavar='DEC', help="Object's declination")],
    'object': ['pos', dict(metavar='OBJECT_NAME', help="Object's name")],
}

def _getval_array(pathlist, key, ext):
    return np.array([fits.getval(img, key, ext) for img in pathlist])

def _getdata_array(pathlist, ext):
    return np.array([fits.getdata(img, ext=ext) for img in pathlist])

def _getheader_array(pathlist, ext):
    return np.array([fits.getheader(img, ext=ext) for img in pathlist])

def _getval_array_mem(hdull, key):
    return np.array([hdul.header.get(key) for hdul in hdull])

def _getdata_array_mem(hdull):
    return np.array([hdul.data for hdul in hdull])

def _getheader_array_mem(hdull):
    return np.array([hdul.header for hdul in hdull])

def _get_band_info_array(prop):
    return np.array([get_band_info(band)[prop] for band in BANDS])

class SCubes:
    def __init__(self, ra, dec, field, size=None, username=None, password=None, verbose=1):
        self.conn = Core(username, password, verbose=verbose)
        self.field = field
        self.verbose = verbose
        self.cubepath = None
        self._stamp_config(ra, dec, size)

    def _stamp_config(self, ra, dec, size):
        self.ra, self.dec = convert_coord_to_degrees(ra, dec)
        self.size = size

    def _getval(self, obj, key, ext, mem=False):
        return _getval_array_mem(obj, key) if mem else _getval_array(obj, key, ext)

    def _getdata(self, obj, ext, mem=False):
        return _getdata_array_mem(obj) if mem else _getdata_array(obj, ext)

    def _getheader(self, obj, ext, mem=False):
        return _getheader_array_mem(obj) if mem else _getheader_array(obj, ext)

    def _download_calibrated_stamps(self, objname, outpath=None, force=False):
        images = []
        wimages = []
        _kw_args = dict(ra=self.ra, dec=self.dec, size=self.size, field_name=self.field)
        for b in tqdm(BANDS, desc=f'{objname} @ {self.field} - Downloading', leave=True, position=0):
            b = b.upper().replace('J0', 'F')
            kw_args = _kw_args.copy()
            kw_args.update(band=b, weight=False)            
            if outpath is not None:
                filename = f'{objname}_{self.field}_{b}_{self.size}x{self.size}_swp.fits.fz'
                kw_args.update(outfile=join(outpath, filename))
                _ = self.conn.calibrated_stamp(**kw_args)
                images.append(kw_args['outfile'])
            else:
                images.append(self.conn.calibrated_stamp(**kw_args))
            # wei
            kw_args['weight'] = True
            if outpath is not None:
                kw_args['outfile'] = join(outpath, filename.replace('swp', 'swpweight'))
                _ = self.conn.stamp(**kw_args)
                wimages.append(kw_args['outfile'])
            else:
                wimages.append(self.conn.stamp(**kw_args)[1])
        self.images = images
        self.wimages = wimages

    def _photospectra(self, flam_scale=None, ext=1):
        flam_scale = 1e-19 if flam_scale is None else flam_scale
        _c = const.c
        scale = (1/flam_scale)
        #Jy to to erg/s/cm/cm/Hz
        #Jy2fnu = - 2.5*(np.log10(3631) - 23)  # 48.5999343777177...
        Jy2fnu = 3631e-23

        self.wl__b = _get_band_info_array('pivot_wave')*u.Angstrom

        self.flam_unit = u.erg / u.s / u.cm / u.cm / u.AA
        self.fnu_unit = u.erg / u.s / u.cm / u.cm / u.Hz

        mem = True if self.cubepath is None else False
        calib_data__byx = self._getdata(self.images, ext, mem=mem)
        fnu__byx = calib_data__byx*Jy2fnu*self.fnu_unit
        flam__byx = scale*(fnu__byx*_c/self.wl__b[:, None, None]**2).to(self.flam_unit)

        gain__b = self._getval(self.images, 'GAIN', ext, mem=mem)
        gain__byx = gain__b[:, None, None]
        weidata__byx = np.abs(self._getdata(self.wimages, ext, mem=mem))
        calib_data_err__byx = np.sqrt(1/weidata__byx + np.abs(calib_data__byx)/gain__byx)
        efnu__byx = calib_data_err__byx*Jy2fnu*self.fnu_unit
        eflam__byx = scale*(efnu__byx*_c/self.wl__b[:, None, None]**2).to(self.flam_unit)

        self.flam__byx = flam__byx
        self.eflam__byx = eflam__byx
        self.weidata__byx = weidata__byx

    def _stamp_WCS_to_cube_header(self, header):
        '''
        Convert WCS information from stamp to cube header.

        Parameters
        ----------
        header : :class:`~astropy.io.fits.Header`
            FITS header containing WCS information.

        Returns
        -------
        :class:`~astropy.io.fits.Header`
            Cube header with updated WCS information.
        '''        
        w = WCS(header)
        nw = WCS(naxis=3)
        nw.wcs.cdelt[:2] = w.wcs.cdelt
        nw.wcs.crval[:2] = w.wcs.crval
        nw.wcs.crpix[:2] = w.wcs.crpix
        nw.wcs.ctype[0] = w.wcs.ctype[0]
        nw.wcs.ctype[1] = w.wcs.ctype[1]
        try:
            nw.wcs.pc[:2, :2] = w.wcs.pc
        except:
            pass
        return nw.to_header()
    
    def _weights_mask_hdu(self):
        # WEIGHTS MASK HDU
        w__byx = self.weidata__byx
        wmask__byx = np.where(w__byx < 0, 1, 0)
        wmask__yx = wmask__byx.sum(axis=0)
        wmask_hdu = fits.ImageHDU(wmask__yx)
        wmask_hdu.header['EXTNAME'] = ('WEIMASK', 'Sum of negative weight pixels (from 1 to 12)')
        return wmask_hdu
    
    def _metadata_hdu(self, ext=1):
        # METADATA HDU
        tab = [BANDS]
        names = ['filter', 'central_wave', 'pivot_wave', 'PSFFWHM']
        for item in names[1:-1]:
            tab.append(_get_band_info_array(item))      
        psffwhm__b = []
        for img in self.images:
            if self.cubepath is not None:
                hdr = fits.getheader(img, ext=ext)
            else:
                hdr = img.header
            key = [k for k in hdr.keys() if 'FWHMMEAN' in k]
            if len(key) == 1:
                psffwhm__b.append(hdr.get(key[0]))
        tab.append(psffwhm__b)       
        meta_hdu = fits.BinTableHDU(Table(tab, names=names))
        meta_hdu.header['EXTNAME'] = 'METADATA'
        return meta_hdu

    def _create_cube_hdulist(self, objname, ext=1):
        cube_prim_hdu = fits.PrimaryHDU()
        cube_prim_hdu.header['TILE'] = self.field
        cube_prim_hdu.header['OBJECT'] = objname
        cube_prim_hdu.header['SIZE'] = (self.size, 'Side of the stamp in pixels')
        cube_prim_hdu.header['RA'] = self.ra
        cube_prim_hdu.header['DEC'] = self.dec
        cube_prim_hdu.header.update(
            self._stamp_WCS_to_cube_header(
                self.images[0].header if self.cubepath is None else fits.getheader(self.images[0], ext=0)
            )
        )
        # DATA HDU
        flam_hdu = fits.ImageHDU(self.flam__byx.value, cube_prim_hdu.header)
        flam_hdu.header['EXTNAME'] = ('DATA', 'Name of the extension')       
        # ERRORS HDU
        eflam_hdu = fits.ImageHDU(self.eflam__byx.value, cube_prim_hdu.header)
        eflam_hdu.header['EXTNAME'] = ('ERRORS', 'Name of the extension')
        hdul = [cube_prim_hdu, flam_hdu, eflam_hdu]
        # INFO TO HEADERS
        for hdu in hdul[1:]:
            hdu.header['BSCALE'] = (self.flam_scale, 'Linear factor in scaling equation')
            hdu.header['BZERO'] = (0, 'Zero point in scaling equation') 
            hdu.header['BUNIT'] = (f'{self.flam_unit}', 'Physical units of the array values')
        hdul.append(self._weights_mask_hdu())
        hdul.append(self._metadata_hdu(ext))
        return hdul

    def write(self, cubepath, overwrite=False):
        print_level(f'writting cube {cubepath}', 1, self.verbose)
        fits.HDUList(self.cube).writeto(cubepath, overwrite=overwrite)
        print_level(f'Cube successfully created!', 1, self.verbose)    

    def create_cube(self, flam_scale=None, objname=None, outpath=None, redo=False, data_ext=1, write_fits=False):
        self.flam_scale = 1e-19 if flam_scale is None else flam_scale
        self.objname = 'myobj' if objname is None else objname
        self.outpath = '.' if outpath is None else outpath

        if outpath is not None:
            try: 
                makedirs(self.outpath)
            except FileExistsError:
                print_level(f'{self.outpath}: directory already exists', 1, self.verbose)    

            self.cubepath = join(self.outpath, f'{self.objname}_cube.fits')

            if exists(self.cubepath) and not redo:
                raise OSError('SCube exists!')

        self._download_calibrated_stamps(objname, outpath, force=redo)
        self._photospectra(flam_scale, ext=data_ext)

        self.cube = self._create_cube_hdulist(objname, ext=data_ext)

        if (self.cubepath is not None) and write_fits:
            self.write(self.cubepath, redo)
           
def scubes_argparse(args):
    '''
    A particular parser of the command-line arguments for `scubes` entry-point script.

    Parameters
    ----------
    args : :class:`argparse.Namespace`
        Command-line arguments parsed by :meth:`argparse.ArgumentParser.parse_args`

    Returns
    -------
    :class:`argparse.Namespace`
        Command-line arguments parsed.
    '''
    # closest even
    args.size = round(args.size/2)*2
    return args

def scubes():
    '''
    Entry-point function for creating S-PLUS galaxy data cubes (S-CUBES).

    Raises
    ------
    SystemExit
        If SExtractor is not found.

    Returns
    -------
    None
    '''
    from splusdata.scubes.core import SCubes

    parser = create_parser(args_dict=SCUBES_ARGS, program_description=SCUBES_PROG_DESC)
    # ADD VERSION OPTION
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__script_version__))
    args = scubes_argparse(parser.parse_args(args=sys.argv[1:]))

    if args.verbose == 0:
        # dactivate warnings without verbosity
        import warnings
        from astropy.wcs import FITSFixedWarning
        from astropy.io.fits.verify import VerifyWarning

        warnings.simplefilter('ignore', category=VerifyWarning)
        warnings.simplefilter('ignore', category=FITSFixedWarning)

    scube = SCubes(
        ra=args.ra, dec=args.dec, field=args.field, 
        size=args.size, username=args.username, password=args.password,
        verbose=args.verbose,
    )

    outpath = join(args.workdir, args.object)
    scube.create_cube(objname=args.object, outpath=outpath, redo=args.redo, data_ext=1, write_fits=True)