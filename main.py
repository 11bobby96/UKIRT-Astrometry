"""
Run the full UKIRT astrometric preprocessing pipeline for a single target.

This script:
1) Generates new WCS solutions for raw FITS images,
2) Centroids sources using CASU imcore,
3) Queries Gaia sources overlapping each image.
"""

from get_new_wcs import get_new_wcs
from centroid_image_stars import centroid_image_stars
from find_gaia_stars import find_gaia_stars

target = 'W0436+1901'
base_dir = '/Users/bobbystiller/Documents/UKIRT_Bobby/UKIRT_data'
save_dir = '/Users/bobbystiller/Desktop/test'

if __name__ == '__main__':

    get_new_wcs(
        data_dir=f'{base_dir}/{target}',
        output_dir=f'{save_dir}/data',
        frame=3
    )

    centroid_image_stars(
        fits_dir=f'{save_dir}/data/{target}',
        output_dir=f'{save_dir}/centroids'
    )

    find_gaia_stars(
        fits_dir=f'{save_dir}/data/{target}',
        output_dir=f'{save_dir}/gaia'
    )
