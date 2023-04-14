#!/bin/bash

# arguments are source name, RA, DEC, obslen

/home/sprite/anaconda3/envs/r2_casper/bin/python /home/sprite/installs/sprite_correlator_sw/src/corr_grab_h5_updated_v4.py -- /home/sprite/installs/sprite_correlator_sw/config/sprite_lowband.yaml ${1} ${2} ${3} ${4} rp  &

/home/sprite/anaconda3/envs/r2_casper/bin/python /home/sprite/installs/sprite_correlator_sw/src/corr_grab_h5_updated_v4.py -- /home/sprite/installs/sprite_correlator_sw/config/sprite_highband.yaml ${1} ${2} ${3} ${4} &


