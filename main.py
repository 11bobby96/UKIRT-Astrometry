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

test_target = 'W0436+1901'

if __name__ == '__main__':

    get_new_wcs(data_dir=f'/Users/bobbystiller/Documents/UKIRT_Bobby/UKIRT_data/{test_target}',
                       output_dir='/Users/bobbystiller/Desktop/test/data',
                       frame=3)

    centroid_image_stars(fits_dir=f'/Users/bobbystiller/Desktop/test/data/{test_target}',
                         output_dir='/Users/bobbystiller/Desktop/test/centroids')

    find_gaia_stars(fits_dir=f'/Users/bobbystiller/Desktop/test/data/{test_target}',
                    output_dir='/Users/bobbystiller/Desktop/test/gaia')
