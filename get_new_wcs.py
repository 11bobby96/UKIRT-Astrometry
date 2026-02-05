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
    data_dir
        Base directory containing semester subdirectories (e.g., U* folders).
    output_dir
        Destination directory for WCS-updated FITS images.
    frame
        HDU index to extract from each FITS file (e.g., 3 for WFCAM science frame).

    Returns
    -------
    None
    """
    semesters = sorted(glob(f'{data_dir}/U*'))
    for semester in tqdm(semesters):
        files = sorted(glob(f'{semester}/*.fit'))

        for file in tqdm(files):
            hdu = save_target_frame(file, output_dir, frame)
            wcs = plate_solution(file)
            update_wcs_header(hdu, wcs)
            break


def save_target_frame(file, output_dir, frame):
    """
    Extract a specific HDU from a FITS file and save it as a new Primary HDU.

    The primary header (HDU 0) and the target extension header are merged.
    The resulting image is written to a directory structure organized by
    target and semester. The header keyword 'CC_PRES' is removed if present
    to avoid conflicts.

    Parameters
    ----------
    file
        Path to the input FITS file.
    output_dir
        Base directory where extracted frames will be written.
    frame
        HDU index to extract.

    Returns
    -------
    out_path
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


def plate_solution(image_path):
    """
    Extract a specific HDU from a FITS file and save it as a new Primary HDU.

    The primary header (HDU 0) and the target extension header are merged.
    The resulting image is written to a directory structure organized by
    target and semester.

    Parameters
    ----------
    file
        Path to the input FITS file.
    output_dir
        Base directory where extracted frames will be written.
    frame
        HDU index to extract.

    Returns
    -------
    out_path
        Path to the newly written FITS file.
    """
    ast = AstrometryNet()
    ast.api_key = 'zemvfxdxbnimplin'
    wcs_header = ast.solve_from_image(image_path, solve_timeout=1000)
    return wcs_header


def update_wcs_header(image_path, wcs_header):
    """
    Update a FITS file header with a new WCS solution. The recently created
    FITS file is modified in place.

    Parameters
    ----------
    image_path
        Path to the FITS file whose header will be updated.
    wcs_header
        WCS header returned by `plate_solution`.

    Returns
    -------
    None
    """
    with fits.open(image_path, mode='update') as hdulist:
        hdulist[0].header.update(wcs_header)
        hdulist.flush()
