"""
Run the full UKIRT astrometric preprocessing pipeline for a single target.

This script:
1) Generates new WCS solutions for raw FITS images,
2) Centroids sources using CASU imcore,
3) Queries Gaia sources overlapping each image.
"""

import logging
import warnings
from astropy import log

warnings.filterwarnings("ignore")
log.setLevel(logging.CRITICAL)

from get_new_wcs import get_new_wcs
from centroid_image_stars import centroid_image_stars
from find_gaia_stars import find_gaia_stars

# target = 'W0310+2046'  # frame: 3, status: 0/3
# target = 'W0338+1717'  # frame: 3, status: 0/3
# target = 'W0353+0418'  # frame: 2, status: 3/3
# target = 'W0412+1044'  # frame: 3, status: 0/3
# target = 'W0419+2036'  # frame: 3, status: 0/3
# target = 'W0427+0743'  # frame: 3, status: 0/3
# target = 'W0430+1058'  # frame: 3, status: 0/3
target = 'W0436+1901'  # frame: 3, status: 3/3
# target = 'W0439+2025'  # frame: 3, status: 0/3
# target = 'W0441+2130'  # frame: 3, status: 0/3
# target = 'W0532+1119'  # frame: 3, status: 0/3
base_dir = '/Users/bobbystiller/Documents/UKIRT_Bobby/UKIRT_data'
save_dir = '/Users/bobbystiller/Documents/sandbox'

if __name__ == '__main__':

    # get_new_wcs(
    #     data_dir=f'{base_dir}/{target}',
    #     output_dir=f'{save_dir}/data',
    #     frame=3
    # )
    #
    # centroid_image_stars(
    #     fits_dir=f'{save_dir}/data/{target}',
    #     output_dir=f'{save_dir}/centroids'
    # )
    #
    # find_gaia_stars(
    #     fits_dir=f'{save_dir}/data/{target}',
    #     output_dir=f'{save_dir}/gaia_2016'
    # )
    print('hi')
