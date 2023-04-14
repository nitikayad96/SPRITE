import time
import datetime
import numpy as np

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord, AltAz

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