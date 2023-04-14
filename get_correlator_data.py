import sys
import os
import numpy as np
from astropy.time import Time

import sprite_utils as su

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Get data from SPRITE correlator')
    parser.add_argument('--sched_file', '-s', type=str,
                        help='Schedule file detailing observation')
    parser.add_argument('--catalog_file', '-c', type=str,
                        help='Catalog file with all source names, RAs, and Decs')
    parser.add_argument('--no_load_switch', action='store_true',
                        help='Use this to NOT move ambient load into place periodically during observation')
    parser.add_argument('--no_alt_check', action='store_true',
                        help='Do not check whether or not sources are always above 30 degrees altitude')

    args = parser.parse_args()

    if 'no_load_switch':
        switch_load = False
    else:
        switch_load = True

    if 'no_alt_check':
        check_src_alt = False
    else:
        check_src_alt = True


    obs_info_list = su.parse_schedule(args.sched_file, args.catalog_file, check_src_alt)
    print(obs_info_list)

    config_lowband = '/home/sprite/SPRITE/sprite_correlator_sw/config/sprite_lowband.yaml'
    config_highband = '/home/sprite/SPRITE/sprite_correlator_sw/config/sprite_highband.yaml'

    grab_h5_file = '/home/sprite/SPRITE/sprite_correlator_sw/src/corr_grab_h5.py'
    move_cal_load_file = '/home/sprite/SPRITE/move_cal_load.py'

    for obs_info in obs_info_list:

        source_name, start_time, ra, dec, obslen = obs_info

        cal_load_iter = int(np.ceil(obslen/(float(60*5))))
        print(cal_load_iter)

        lines = ['#!/bin/bash','\n\n']
        lines.append('python %s -- %s %s %s %s %s rp &'%(grab_h5_file, config_lowband, source_name, ra, dec, obslen))
        lines.append('\n\n')
        lines.append('python %s -- %s %s %s %s %s &' % (grab_h5_file, config_highband, source_name, ra, dec, obslen))
        lines.append('\n\n')
        if switch_load:
            lines.append('python %s %d &' %(move_cal_load_file, cal_load_iter))
            lines.append('\n\n')

        bash_file = '/home/sprite/_get_data.bash'
        with open('%s' %(bash_file), 'w') as f:
            f.writelines(lines)

        cmd = 'bash %s' %bash_file

        while start_time > Time.now():
            continue

        os.system(cmd)

