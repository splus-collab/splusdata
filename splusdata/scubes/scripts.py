import sys
from os import getcwd
from os.path import join

from splusdata.scripts.args import create_parser

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

    creator = SCubes(
        ra=args.ra, dec=args.dec, field=args.field, 
        size=args.size, username=args.username, password=args.password,
        verbose=args.verbose,
    )

    outpath = join(args.workdir, args.object)
    _ = creator.create_cube(objname=args.object, outpath=outpath, force=args.force, data_ext=1, write_fits=True)