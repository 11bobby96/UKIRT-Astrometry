import ssl
from glob import glob
from tqdm import tqdm
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord


def find_gaia_stars(fits_dir, output_dir):
    """
    Query Gaia sources overlapping each FITS image and save results to CSV.

    This function iterates over semester directories in `fits_dir`, loads each
    FITS image and its WCS, estimates the sky footprint of the image, queries
    the Gaia archive for sources within that region, and writes the resulting
    catalog to a per-image CSV file under `output_dir`.

    Parameters
    ----------
    fits_dir
        Base directory containing semester subdirectories with FITS images.
    output_dir
        Destination root directory where Gaia query CSV files will be written.

    Returns
    -------
    None
    """
    semesters = sorted(glob(f'{fits_dir}/U*'))
    output_root = Path(output_dir)
    for s in semesters:
        files = glob(f'{s}/*.fits')
        files.sort()

        for file in tqdm(files):
            file = Path(file)

            semester = file.parts[-3]
            target = file.parts[-2]

            data = fits.getdata(file, hdu=0)
            wcs = WCS(fits.getheader(file, hdu=0))

            central_ra, central_dec, width, height = query_region(data, wcs)
            results = gaia_query(central_ra, central_dec, width, height)

            gaia_save_path = output_root / semester / target
            gaia_save_path.mkdir(parents=True, exist_ok=True)

            output_csv = gaia_save_path / f'{file.stem}.csv'

            df_gaia = results.to_pandas()
            df_gaia.to_csv(output_csv, index=True)


def query_region(image, wcs_header):
    """
    Estimate the central sky position and angular size of a FITS image.

    Computes the sky coords of the image center and the four corners using the
    given WCS, then derives an approximate rectangular query region in RA/Dec.

    Parameters
    ----------
    image
        2D image array used to determine the image center in pixel coordinates.
    wcs_header
        WCS object used to convert pixel coordinates to sky coordinates.

    Returns
    -------
    central_ra
        Right ascension of the image center (degrees).
    central_dec
        Declination of the image center (degrees).
    width
        Approximate angular width of the image footprint in RA (degrees).
    height
        Approximate angular height of the image footprint in Dec (degrees).
    """
    center = wcs_header.pixel_to_world(np.shape(image)[0] / 2, np.shape(image)[1] / 2)
    top_right = wcs_header.pixel_to_world(0, 0)
    top_left = wcs_header.pixel_to_world(0, 2047)
    bottom_right = wcs_header.pixel_to_world(2047, 0)
    bottom_left = wcs_header.pixel_to_world(2047, 2047)

    min_ra = min(top_right.ra.value, top_left.ra.value, bottom_right.ra.value, bottom_left.ra.value)
    max_ra = max(top_right.ra.value, top_left.ra.value, bottom_right.ra.value, bottom_left.ra.value)
    min_dec = min(top_right.dec.value, top_left.dec.value, bottom_right.dec.value, bottom_left.dec.value)
    max_dec = max(top_right.dec.value, top_left.dec.value, bottom_right.dec.value, bottom_left.dec.value)

    width = abs(max_ra - min_ra)
    height = abs(max_dec - min_dec)
    return center.ra.value, center.dec.value, width, height


def gaia_query(central_ra, central_dec, width, height):
    """
    Query the Gaia archive for sources in a rectangular region on the sky.

    Builds a sky coordinate from the provided center position and performs a
    rectangular region query with the given width and height. A small padding
    is added to both dimensions to reduce edge losses.

    Parameters
    ----------
    central_ra
        Right ascension of the query center (degrees).
    central_dec
        Declination of the query center (degrees).
    width
        Angular width of the query region (degrees).
    height
        Angular height of the query region (degrees).

    Returns
    -------
    results
        Gaia query results table returned by the Gaia archive query.
    """
    Gaia.ROW_LIMIT = -1
    coord = SkyCoord(ra=central_ra, dec=central_dec, unit=(u.degree, u.degree), frame='icrs')
    width = u.Quantity(width + 0.02, u.deg)
    height = u.Quantity(height + 0.02, u.deg)
    ssl._create_default_https_context = ssl._create_unverified_context
    r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
    return r
