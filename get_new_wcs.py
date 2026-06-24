import time
from glob import glob
from tqdm import tqdm
from pathlib import Path
from astropy.io import fits
from astroquery.astrometry_net import AstrometryNet


def get_new_wcs(data_dir, output_dir, frame):
    """
    Generate new WCS solutions for FITS images and write updated files.
    
    This function iterates through semester directories within `data_dir`,
    extracts a specified HDU from each FITS file, writes it to `output_dir`,
    solves for a new astrometric plate solution using Astrometry.net,
    and updates the output FITS header with the resulting WCS.

    Parameters
    ----------
    data_dir : str or pathlib.Path
        Base directory containing semester subdirectories, such as U* folders.
    output_dir : str or pathlib.Path
        Destination directory for WCS-updated FITS images.
    frame : int
        HDU index to extract from each FITS file, such as 3 for a WFCAM
        science frame.

    Returns
    -------
    None
    """
    semesters = sorted(glob(f'{data_dir}/U*'))

    for i, semester in enumerate(tqdm(semesters)):
        files = sorted(glob(f'{semester}/*.fit'))

        for file in tqdm(files):
            hdu = _save_target_frame(file, output_dir, frame)
            wcs = _plate_solution(hdu)
            _update_wcs_header(hdu, wcs)


def _save_target_frame(file, output_dir, frame):
    """
    Extract a specific HDU from a FITS file and save it as a new Primary HDU.

    The primary header, HDU 0, and the target extension header are merged.
    The resulting image is written to a directory structure organized by
    target and semester. The header keyword 'CC_PRES' is removed if present
    to avoid conflicts.

    Parameters
    ----------
    file : str or pathlib.Path
        Path to the input FITS file.
    output_dir : str or pathlib.Path
        Base directory where extracted frames will be written.
    frame : int
        HDU index to extract from the input FITS file.

    Returns
    -------
    out_path : pathlib.Path
        Path to the newly written FITS file.
    """
    hdulist = fits.open(file)

    image_data = hdulist[frame].data
    header_hdu0 = hdulist[0].header
    header_hdu3 = hdulist[frame].header

    merged_header = header_hdu0.copy()
    merged_header.update(header_hdu3)
    if 'CC_PRES' in merged_header:
        del merged_header['CC_PRES']

    hdu = fits.PrimaryHDU(data=image_data, header=merged_header)

    file = Path(file)
    target = file.parts[-3]
    semester = file.parts[-2]
    image_name = file.stem

    out_dir = Path(output_dir) / target / semester
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / f'{image_name}_hdu{frame}.fits'
    hdu.writeto(out_path, overwrite=True)
    hdulist.close()
    return out_path


def _plate_solution(image_path, max_attempts=3):
    """
    Solve for an astrometric WCS solution using Astrometry.net.

    This function submits the FITS file at `image_path` to Astrometry.net and
    returns the resulting WCS header. If the plate solve fails, the function
    retries up to `max_attempts` times before raising the final exception.

    Parameters
    ----------
    image_path : str or pathlib.Path
        Path to the input FITS image to be solved.
    max_attempts : int, optional
        Maximum number of plate-solve attempts. Default is 3.

    Returns
    -------
    wcs_header : astropy.io.fits.Header
        WCS header solution returned by Astrometry.net.
    """
    ast = AstrometryNet()
    ast.api_key = 'zemvfxdxbnimplin'

    for attempt in range(max_attempts):

        try:
            wcs_header = ast.solve_from_image(image_path, solve_timeout=1000)
            return wcs_header

        except Exception as e:

            print(f"\nPlate solve failed (attempt {attempt+1}/{max_attempts})")
            print(e)

            if attempt < max_attempts - 1:
                print("Retrying in 10 seconds...")
                time.sleep(10)
            else:
                raise


def _update_wcs_header(image_path, wcs_header):
    """

    Update a FITS file header with a new WCS solution.
    The FITS file at `image_path` is modified in place. The new WCS keywords
    from `wcs_header` are added to the primary HDU header.

    Parameters
    ----------
    image_path : str or pathlib.Path
        Path to the FITS file whose header will be updated.
    wcs_header : astropy.io.fits.Header
        WCS header returned by `plate_solution`.

    Returns
    -------
    None
    """
    with fits.open(image_path, mode='update') as hdulist:
        hdulist[0].header.update(wcs_header)
        hdulist.flush()
