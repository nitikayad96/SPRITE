import numpy as np
import adc5g as adc
import corr
import time
import sys
import struct
import pylab
import socket

def snap(r,name,format='L',man_trig=True):
    n_bytes = struct.calcsize('=%s'%format)
    d = r.snapshot_get(name, man_trig=man_trig)
    return np.array(struct.unpack('>%d%s'%(d['length']/n_bytes,format),d['data']))

def uint2int(d,bits,bp):
    dout = np.array(d,dtype=float)
    dout[dout>(2**(bits-1))] -= 2**bits
    dout /= 2**bp
    return dout

def dbs(x):
    return 10*np.log10(x)



if __name__ == '__main__':
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('%prog [options]')
    p.set_description(__doc__)
    p.add_option('-p', '--skip_prog', dest='prog_fpga',action='store_false', default=True, 
        help='Skip FPGA programming (assumes already programmed).  Default: program the FPGAs')
    #p.add_option('-e', '--skip_eq', dest='prog_eq',action='store_false', default=True, 
    #    help='Skip configuration of the equaliser in the F engines.  Default: set the EQ according to config file.')
    p.add_option('-v', '--verbosity', dest='verbosity',type='int', default=0, 
        help='Verbosity level. Default: 0')
    p.add_option('-r', '--roach', dest='roach',type='str', default='192.168.0.111', 
        help='ROACH IP address or hostname. Default: 192.168.0.111')
    p.add_option('-s', '--secondroach', dest='secondroach',type='str', default='192.168.0.111', 
        help='second ROACH IP address or hostname. Default: 192.168.0.111')
    p.add_option('-b', '--boffile', dest='boffile',type='str', default='ami_fx_sbl_wide.bof', 
        help='Boffile to program. Default: ami_fx_sbl_wide.bof')
    p.add_option('-a', '--acc_len', dest='acc_len',type='int', default='1024', 
        help='Number of spectra to accumulate. Default: 1024')
    p.add_option('-f', '--fft_shift', dest='fft_shift',type='str', default='111111111111', 
        help='FFT shift schedule. Enter as a 12-bit binary string. Default: 111111111111 (i.e. shift every stage)')
    p.add_option('-t', '--tvg', dest='tvg',action='store_true', default=False, 
        help='Use corner turn tvg. Default:False')
    p.add_option('-m', '--manual_sync', dest='manual_sync',action='store_true', default=False, 
        help='Use this flag to issue a manual sync (useful when no PPS is connected). Default: Do not issue sync')
    p.add_option('-n', '--network', dest='network',action='store_true', default=False, 
        help='Send data out over tcp')
    p.add_option('-l', '--plottting', dest='plotting',action='store_true', default=False,
                 help='Do plotting')

    opts, args = p.parse_args(sys.argv[1:])


    if opts.network:
        #set up the socket
        TCP_IP = '127.0.0.1'
        TCP_PORT = 10000
        BUFFER_SIZE = 1024
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((TCP_IP,TCP_PORT))

    print 'Connecting to %s'%opts.roach
    r = corr.katcp_wrapper.FpgaClient(opts.roach)
    time.sleep(0.2)
    print 'ROACH is connected?', r.is_connected()

    print 'Connecting to %s'%opts.secondroach
    r2 = corr.katcp_wrapper.FpgaClient(opts.secondroach)
    time.sleep(0.2)
    print 'Second ROACH is connected?', r2.is_connected()

    
    if opts.prog_fpga:
        print 'Programming ROACHes with boffile %s'%opts.boffile
        r.progdev(opts.boffile)
        time.sleep(0.5)
        r2.progdev(opts.boffile)
        time.sleep(0.5)
        print 'Estimating clock speed...'
        print 'Clock speed is %d MHz'%r.est_brd_clk()
        print 'Second clock speed is %d MHz'%r2.est_brd_clk()
        print 'Calibrating ADCs'

        adc.calibrate_all_delays(r,0,snaps=['snapshot_adc0'],verbosity=opts.verbosity)
        adc.calibrate_all_delays(r,1,snaps=['snapshot_adc1'],verbosity=opts.verbosity)

        adc.calibrate_all_delays(r2,0,snaps=['snapshot_adc0'],verbosity=opts.verbosity)
        adc.calibrate_all_delays(r2,1,snaps=['snapshot_adc1'],verbosity=opts.verbosity)

        
    print 'reseting control sw'
    r.write_int('control',0)
    print 'Setting accumulation length to %d'%opts.acc_len
    r.write_int('acc_len',opts.acc_len)
    print 'reseting second control sw'
    r2.write_int('control',0)
    print 'Setting second accumulation length to %d'%opts.acc_len
    r2.write_int('acc_len',opts.acc_len)

    print 'Setting fft-shift to %s'%opts.fft_shift
    r.write_int('fft_shift0',int(opts.fft_shift,2))
    r.write_int('fft_shift1',int(opts.fft_shift,2))
    print 'Setting second fft-shift to %s'%opts.fft_shift
    r2.write_int('fft_shift0',int(opts.fft_shift,2))
    r2.write_int('fft_shift1',int(opts.fft_shift,2))

    COARSE_DELAY=10
    print 'Setting coarse delays to %d'%COARSE_DELAY
    r.write_int('coarse_delay0',COARSE_DELAY)
    r.write_int('coarse_delay1',COARSE_DELAY)
    print 'Setting second coarse delays to %d'%COARSE_DELAY
    r2.write_int('coarse_delay0',COARSE_DELAY)
    r2.write_int('coarse_delay1',COARSE_DELAY)

    # critical step to get both ROACHes in sync
    
    print 'Arming pps'
    ctrl = r.read_uint('control')
    ctrl2 = r2.read_uint('control')
    ctrl = ctrl | (1<<2)
    ctrl2 = ctrl2 | (1<<2)
    r.write_int('control',ctrl)
    r2.write_int('control',ctrl2)
    time.sleep(1.5)
    ctrl = ctrl & ((2**32-1)-(1<<2))
    ctrl2 = ctrl2 & ((2**32-1)-(1<<2))
    r.write_int('control',ctrl)
    r2.write_int('control',ctrl2)

