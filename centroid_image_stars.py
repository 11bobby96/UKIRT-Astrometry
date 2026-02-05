import os
import subprocess
import pandas as pd
from glob import glob
from tqdm import tqdm
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS


def centroid_image_stars(fits_dir, output_dir):
    """
    Run source centroiding on a directory of FITS images.

    This function searches for WCS-calibrated FITS images, runs CASU's
    imcore source extraction on each file, parses the resulting ellipse
    catalog, converts centroids to world coordinates using the image WCS,
    and saves the results to CSV files. Temporary .ell and .cat files
    produced by imcore are deleted after conversion to CSV. Failures are
    reported but do not stop execution.

    Parameters
    ----------
    fits_dir
        Base directory containing FITS images organized by target and semester.
    output_dir
        Destination directory where centroid catalogs will be written.

    Returns
    -------
    None
    """
    fits_files = sorted(glob(f'{fits_dir}/*/*_hdu*.fits'))
    output_dir = Path(output_dir)

    for fits_path in tqdm(fits_files):
        try:
            fits_path = Path(fits_path)
            target = fits_path.parts[-3]
            semester = fits_path.parts[-2]

            output_path = output_dir / target / semester
            output_path.mkdir(parents=True, exist_ok=True)

            ell_path, cat_path = run_imcore_on_fits(fits_path, output_path)
            df = parse_casu_ell_file(ell_path, fits_path)

            csv_path = ell_path.replace('.ell', '.csv')
            df.to_csv(csv_path, index=False)

            os.remove(ell_path)
            os.remove(cat_path)

        except Exception as e:
            print(f'Failed on {fits_path}: {e}')


def run_imcore_on_fits(fits_path, output_dir, ipix=2, thresh=3, cattype=2):
    """
    Run CASU imcore source extraction on a FITS image.

    Parameters
    ----------
    fits_path
        Path to the FITS image to process.
    output_dir
        Directory where output catalog files will be written.
    ipix
        Minimum number of pixels above background for an object to be
        considered a detection. I have not put much thought into this
        parameter. This should be optimized in the future.
    thresh
        The isophotal analysis threshold, specified in terms of the number of
        standard deviations above the local background level. The global sigma
        is computed as part of the background estimation. The minimum allowed
        value is 1 and the recommended value for general purpose deep image
        detection is 1.25.
    cattype
        Output catalogue type: 1 == INT WFC, 2 == WFCAM,
        3 == Basic, 4 == Object Mask, 6 == VIRCAM or VST

    Returns
    -------
    ell_path
        Path to the generated ellipse (.ell) catalog.
    cat_path
        Path to the generated source catalog (.cat).

    Notes
    -----
    http://casu.ast.cam.ac.uk/surveys-projects/software-release/imcore
    Uses the confidence map from confidence_map in the analysis. If a filename
    of noconf is specified no confidence map is used in the analysis.
    """
    base = os.path.basename(fits_path).replace('.fits', '').replace('.fit', '')
    cat_path = os.path.join(output_dir, f'{base}.cat')
    ell_path = os.path.join(output_dir, f'{base}.ell')

    cmd = [
        'imcore',
        fits_path,
        'noconf',
        cat_path,
        str(ipix),
        str(thresh),
        f'--cattype={cattype}'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    return ell_path, cat_path


def parse_casu_ell_file(ell_path, fits_path):
    """
    Parse a CASU imcore ellipse file and compute sky coordinates.

    Reads centroid and shape parameters from a CASU .ell file and converts
    pixel coordinates into right ascension and declination using the FITS
    WCS solution.

    Parameters
    ----------
    ell_path
        Path to the CASU ellipse (.ell) file.
    fits_path
        Path to the corresponding FITS image containing a valid WCS.

    Returns
    -------
    dataframe
        Table containing x, y, a, b, theta, raWCS, and decWCS for each source.
    """
    with fits.open(fits_path) as hdul:
        wcs_header = WCS(hdul[0].header)

    x_list = []
    y_list = []
    a_list = []
    b_list = []
    theta_list = []
    ra_list = []
    dec_list = []
    with open(ell_path, 'r') as file:
        for line in file:
            if line.startswith('image; ellipse'):
                parts = line.split()
                x = float(parts[2]) - 1
                y = float(parts[3]) - 1
                a = float(parts[4])
                b = float(parts[5])
                theta = float(parts[6])

                world_coords = wcs_header.pixel_to_world(x, y)
                ra = world_coords.ra.value
                dec = world_coords.dec.value

                x_list.append(x)
                y_list.append(y)
                a_list.append(a)
                b_list.append(b)
                theta_list.append(theta)
                ra_list.append(ra)
                dec_list.append(dec)

    return pd.DataFrame({
        'x': x_list,
        'y': y_list,
        'a': a_list,
        'b': b_list,
        'theta': theta_list,
        'raWCS': ra_list,
        'decWCS': dec_list
    })
