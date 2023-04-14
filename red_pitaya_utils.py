import time
import datetime
import numpy as np
import struct
import casperfpga
import sys

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord, AltAz

#################################################################################################################
# INITIALIZE RED PITAYA (USEFUL FOR CASE OF POWER CYCLE)
#################################################################################################################
def init_rp(bof = 'firmware/loberotator_final_2022-07-14_1419.fpg'):
    """
    Initialize BRAM and required registers for Red Pitaya - need to do this every time board is power cycled
    :param bof: firmware file used to program FPGA - should be either .bof or .fpg
    """
    # connect to board
    rp = casperfpga.CasperFpga('rp-f07179.ovro.pvt')

    # program board with bof file
    rp.upload_to_ram_and_program(bof)

    # check that board is running
    print(rp.is_running())
    print(rp.is_connected())

    # program bram with values of sin and cos
    N = 2 ** 11
    n = np.arange(N) / float(N - 1)

    bram_dict = {'sin': np.sin(2 * np.pi * n), 'cos': np.cos(2 * np.pi * n)}
    bram_32b = {}

    for key in bram_dict.keys():

        scaled_sig_32b = np.empty(N)

        sig = bram_dict[key]
        for i in range(len(sig)):
            scaled_sig_32b[i] = sig[i] * (2 ** 31) - 1

        vals_32b = scaled_sig_32b.astype(int)

        # iterate through all values, pack them as hex data and write each word to correct bram address
        bram_name = 'bram_%s_phi' % key

        for j in range(2048):
            buf_val = struct.pack('>1i', vals_32b[j])
            rp.write(bram_name, buf_val, 4 * j)

        # read back from bram to ensure this worked properly
        bram_32b[key] = np.array(struct.unpack('>2048i', rp.read(bram_name, 2048 * 4, 0)))

    rp.write_int('phi_init', 0 * (2 ** 32))
    rp.write_int('phi_step', int(4e-9 * (2 ** 31)))

    delay = 0.46  # rad
    c1 = -1 / np.sin(delay)
    c2 = np.cos(delay) / np.sin(delay)

    rp.write_int('c1', int(c1 * (2 ** 14)))
    rp.write_int('c2', int(c2 * (2 ** 14)))

    rp.write_int('reg_cntrl', 1)
    rp.write_int('reg_cntrl', 0)


#################################################################################################################
# CALCULATE LOBE ROTATION RATE FOR ADJUSTING VALUES IN RED PITAYA
#################################################################################################################
C = 3e10  # speed of light
B = 24e2  # baseline length
B_VEC = [0, B, 0]  # baseline vector
NU = 90e9  # observing frequency
RP_CLK_T = 1 / (125.22e6)  # clock rate
OVRO = EarthLocation.of_site('ovro')


def s_hat(alt_list, az_list):
    """
    Calculate source direction vectors given an altitude and azimuth list
    :param alt_list: list of altitudes in degrees
    :param az_list: list of azimuths in degrees
    :return: list of direction unit vectors
    :rtype: numpy array
    """

    alt_list = np.deg2rad(alt_list)
    az_list = np.deg2rad(az_list)

    s_hat_list = []
    for i in range(len(alt_list)):
        alt = alt_list[i]
        az = az_list[i]
        s_hat_list.append([np.cos(alt) * np.cos(az), np.cos(alt) * np.sin(az), np.sin(alt)])

    return (np.array(s_hat_list))

def calc_phi_step(ra, dec, start_UTC, t_obs, t_update):
    """
    Calculate phase rate updates in units of fringe cycles per clock cycle to be sent to Red Pitaya for each time \
    interval requested
    :param ra: string of RA of source in HMS
    :param dec: string of Dec of source in DM
    :param start_UTC: string in isot format for beginning of observation
    :param t_obs: length of observation in minutes
    :param t_update: update interval in minutes
    :return: list of phase rates
    :rtype: numpy array
    """
    # create an astropy time object with all times to calculate dphi/dt at
    t_start = Time(start_UTC, format='fits', scale='utc')
    t_intervals = np.arange(0, t_obs, t_update) * u.min
    t_steps = t_start + t_intervals

    # create a sky coordinates object given source RA and Dec and
    # convert to altaz frame
    src = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    altaz_frame = AltAz(obstime=t_steps, location=OVRO)
    altaz = src.transform_to(altaz_frame)

    # using the altaz info, generate source direction vectors and use that to
    # calculate the phase rate in rad/s
    s = s_hat(altaz.alt, altaz.az)
    theta = np.arccos(np.dot(B_VEC, s.T) / float(B))
    dtau_dt = B * np.gradient(theta, t_intervals.to(u.s).value) * np.sin(theta) / float(C)
    dphi_dt = 2 * np.pi * NU * dtau_dt

    # convert dphi_dt from rad/s to fringe cycles/red pitaya clock cycle
    phi_step_list = dphi_dt * RP_CLK_T / float((2 * np.pi))

    phi_step_list = np.array(phi_step_list)

    return phi_step_list

def calc_rp_error(phi_step_list):
    """
    Calculate both actual phase rate and phase rate loaded into Red Pitaya in Hz
    :param phi_step_list: list of phase rates in units of fringe cycle per clock cycle
    :return: list of actual phase rates and list of loaded phase rates
    :rtype: numpy array, numpy array
    """

    actual_phi_rate = phi_step_list / RP_CLK_T

    phi_step_list_32b = np.round(phi_step_list * (2 ** 31))
    rp_phi_rate = ( phi_step_list_32b / (2**31) ) / RP_CLK_T

    return actual_phi_rate, rp_phi_rate