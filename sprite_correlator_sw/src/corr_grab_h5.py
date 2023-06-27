import sys
import os
import time
import numpy as np
import adc5g as adc
import pylab
import socket
import ami as AMI
import helpers as helpers
import control as control
import file_writer as fw
import pylab
import signal
import logging
import casperfpga
import datetime
from astropy.time import Time
import pause

from red_pitaya_utils import calc_phi_step, calc_rp_error

logger = helpers.add_default_log_handlers(logging.getLogger(__name__))

def write_data(writer, d, timestamp, mytime):
    writer.append_data('xeng_raw0', d.shape, d, np.int64)
    writer.append_data('timestamp0', [1], timestamp, float)
    writer.append_data('timestamp_UTC', [1], mytime, float)

def signal_handler(signum, frame):
    """
    Run when kill signals are caught
    """
    print("Received kill signal %d. Closing files and exiting"%signum)
    writer.close_file()
    try:
        ctrl.close_sockets()
    except:
       pass #this is poor form
    exit()


if __name__ == '__main__':

    from optparse import OptionParser

    # set clock
    now = datetime.datetime.now()
    micr = now.microsecond
    delt = datetime.timedelta(seconds=3,microseconds=1000000-micr)
    
    
    ################ Parse Arguments #########################
    p = OptionParser()
    p.set_usage('%prog [options] [CONFIG_FILE]')
    p.set_description(__doc__)
    opts, args = p.parse_args(sys.argv[1:])

    if args == []:
        config_file = None
    else:
        config_file = args[0]

    src_name = args[1]
    ra = args[2]
    dec = args[3]
    obslen = float(args[4])
    track = args[5] # name of target source in track
    purpose = args[6] # one of [science,gain,bandpass,flux]
    data_path = args[7]


    if len(args) == 9 and args[8] == 'rp':
        switch_rp = True

    else:
        switch_rp = False

    if 'lowband' in config_file:
        band = 'low'

    elif 'highband' in config_file:
        band = 'high'
        


    ################# Initialize observation file ######################
    writer = fw.H5Writer(config_file=config_file)
    bl_order = [[0,0], [1,1], [0,1]]
    writer.set_bl_order(bl_order)
    writer.set_source_details(src_name, ra, dec)
    writer.set_metadata(band, track, purpose)
    
    ############## Initialize red pitaya values ########################

    t = Time.now()
    t.format='isot'
    start_utc = t.value
    phi_step_list = calc_phi_step(ra, dec, start_utc, 2, 1)
    actual_phi_rate, rp_phi_rate = calc_rp_error(phi_step_list)
    writer.set_fringe_rates(actual_phi_rate[0], rp_phi_rate[0])

    if switch_rp:
        rp = casperfpga.CasperFpga('rp-f07179.ovro.pvt')
        phi_step_32b = round(phi_step_list[0] * (2 ** 31))
        rp.write_int('phi_step', round(0))
        #rp.write_int('phi_init', 0 * (2 ** 32))

    ################ Inilatize correlator values #######################
    corr = AMI.spriteSbl(config_file=config_file) #passive=True)
    time.sleep(0.1)

    xeng = corr.xengs[0]

    # some initial values for the loop
    datavec = np.zeros([corr.n_chans*corr.n_bands,corr.n_bls,corr.n_pols,2],dtype=np.int64)
    print('Shape of datavec: ', datavec.shape)

    current_obs = None
    receiver_enable = False
    last_meta_timestamp = time.time()
    cnt=0
    ps = [0, 0.25, 0.5, 0.75]

    # Catch keyboard interrupt and kill signals (which are initiated by amisa over ssh)
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    ################### run  data capture loop for given observation time ########################

    if current_obs is None:
        fname = 'corr_TEST_%d.h5'%(time.time())
        writer.start_new_file(fname,data_path)
        current_obs = 'test'

    # start at reasonable time
    pause.until(now+delt)
    if switch_rp:
        rp.write_int('phi_init', 0 * (2 ** 32))

    mcnt_old = xeng.read_uint('mcnt_lsb')

    while cnt<int(obslen/(50000.*8.192e-6)):

        mcnt = xeng.read_uint('mcnt_lsb')
        if mcnt != mcnt_old:
            mytime = time.time()
            mcnt_old = mcnt            
            ps_idx = (cnt+1)%4
            if switch_rp:
                rp.write_int('phi_init', ps[ps_idx] * (2 ** 32))
                if cnt==0:
                    rp.write_int('reg_cntrl', 1)
                    rp.write_int('reg_cntrl', 0)
                    rp.write_int('phi_step', phi_step_32b)                
            d = corr.snap_corr_wide(wait=False,combine_complex=False)

            cnt += 1

            if d is not None:
                datavec[:,0,0,1] = d['corr00']
                datavec[:,1,0,1] = d['corr11']
                datavec[:,2,0,1] = d['corr01'][0:2048] #datavec[:,:,:,1] should be real
                datavec[:,2,0,0] = d['corr01'][2048:] #datavec[:,:,:,0] should be imag
                print("got new correlator data from %s band with timestamp %.4f at time %.4f"%(band, mytime, time.time()))
                write_data(writer, datavec, d['timestamp'], mytime)
                d_old = d
                t_old = mytime
                datavec_old = datavec

            else:
                write_data(writer, datavec_old, d_old['timestamp'], t_old)
                print("Failed to send because MCNT changed during snap")

        time.sleep(0.01)
    


        
