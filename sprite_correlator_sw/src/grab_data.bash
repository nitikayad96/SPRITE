#!/bin/bash

# arguments are source name, RA, DEC, obslen

/home/sprite/anaconda3/envs/r2_casper/bin/python /home/sprite/SPRITE/sprite_correlator_sw/src/corr_grab_h5.py -- /home/sprite/SPRITE/sprite_correlator_sw/config/sprite_lowband.yaml ${1} ${2} ${3} ${4} ${1} science /home/sprite/vikram/scratch/HB rp  &

/home/sprite/anaconda3/envs/r2_casper/bin/python /home/sprite/SPRITE/sprite_correlator_sw/src/corr_grab_h5.py -- /home/sprite/SPRITE/sprite_correlator_sw/config/sprite_highband.yaml ${1} ${2} ${3} ${4} ${1} science /home/sprite/vikram/scratch/LB &


