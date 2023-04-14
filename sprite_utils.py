import sys
import os
import struct
import time
import datetime
import socket

import pandas as pd
import numpy as np

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord, AltAz

#################################################################################################################
# PARSING CATALOG AND SCHEDULE FILES NECESSARY FOR RUNNING CORRELATOR CODE
#################################################################################################################
def cat_to_df(catalog_file):
    """
    Takes a calibrator file used with the COMAP control system and creates a pandas dataframe containing RA and Dec \
    information for each source
    :param catalog_file: catalog file ending in ".cat"
    :return: source details (source name, RA, Dec)
    :rtype: pandas DataFrame
    """
    with open(catalog_file) as f:
        lines = f.readlines()

    columns = ['Source Name', 'RA', 'Dec']

    data = []
    for line in lines:

        if line[0] == '#':
            continue

        else:

            if line[-1:] == '\n':
                line = line[:-1]

            if line[-1] == ' ':
                line = line[:-1]

            epoch, source, ra, dec = line.split(' ')
            data.append([source, ra, dec])

    df = pd.DataFrame(data, columns=columns)

    return df


def get_date():
    """
    Get the current date
    :return: current dat
    :rtype: string
    """
    t = Time.now()
    t.format = 'iso'
    date = t.value.split(' ')[0]
    return date


def check_alt(start_time, ra, dec, obslen):
    """
    Check that the altitude of a source is above 30 degrees for the entire observation time
    :param start_time: astropy time object for date and time in UTC that observation begins
    :param ra: RA of source
    :param dec: Dec of source
    :param obslen: length of observation in seconds
    :return: True if source is always above 30 degrees, False if not
    :rtype: Bool
    """
    obstime = start_time + np.linspace(0, obslen, 100) * u.s
    src = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    altaz_frame = AltAz(obstime=obstime, location=OVRO)
    altaz = src.transform_to(altaz_frame)

    return sum(1 * (altaz.alt.value < 30)) == 0


def parse_schedule(sched_file, catalog_file = 'calibrator_catalog.cat', check_src_alt = True):
    """
    Given a SPRITE schedule file and a source catalog file, create a list of observation details for each source to be \
    passed to the correlator code. Each row in the output array will contain the source name, start time, RA, Dec, and \
    observation length (in that order)
    :param sched_file: schedule file ending in ".sch"
    :param catalog_file: catalog file ending in ".cat", default set to "calibrator_catalog.cat"
    :param check_src_alt: Bool, default True - if True, raise exception for sources below 30 degrees
    :return: list of details from schedule file
    :rtype: numpy array
    """
    with open(sched_file) as f:
        lines = f.readlines()

    obs_info = []
    bad_sources = []
    ref_time = Time.now()
    for line in lines:

        if line[-1] == '\n':
            line = line[:-1]

        if 'observe_source' in line:

            _, source_name, start_time, duration = line.split(' ')

            if duration[-1] == 's':
                obslen = float(duration[:-1])
            elif duration[-1] == 'm':
                obslen = 60 * float(duration[:-1])
            elif duration[-1] == 'h':
                obslen = 3600 * float(duration[:-1])

        elif 'pointing_scan' in line:
            _, source_name, start_time, offset, step_size = line.split(' ')

            offset = offset[:-1]
            deg_off, min_off, sec_off = offset.split(':')
            # arcsec_off = full range for both x and y directions
            arcsec_off = 2 * 2 * (float(deg_off) * 60 * 60 + float(min_off) * 60 + float(sec_off))

            deg_step, min_step, sec_step = step_size.split(':')
            arcsec_step = (float(deg_step) * 60 * 60 + float(min_step) * 60 + float(sec_step))

            n_step = arcsec_off / arcsec_step
            obslen = 10 * n_step + 10  # in seconds

        else:
            continue

        source_name = source_name[:-1]
        start_time = start_time[:-1]

        df = cat_to_df(catalog_file)
        src_row = df[df['Source Name'] == source_name]
        ra = src_row['RA'].values[0]
        dec = src_row['Dec'].values[0]

        start_time = Time('%s %s' % (get_date(), start_time))
        if start_time < ref_time:
            start_time = start_time + 1 * u.day

        source_row = [source_name, start_time, ra, dec, obslen]
        obs_info.append(source_row)

        if not check_alt(*source_row[1:]):
            bad_sources.append(source_name)

    if check_src_alt:
        if len(bad_sources) > 0:
            except_msg = 'The following sources are below an altitude 30 degrees for part/all of this observation:'
            for bad_src in bad_sources:
                except_msg+='\n%s'%bad_src
            raise Exception(except_msg)

    return np.array(obs_info)

#################################################################################################################
# SOCKET CONNECTION FOR TCP/IP COMMUNICATION TO CAN MODULES FOR CONTROLLING HOT LOAD
#################################################################################################################

# define connection details
HOST1 = 'c1.ovro.pvt'
HOST2 = 'c2.ovro.pvt'
PORT1 = 15000
PORT2 = 15001

# create byte arrays for switching between sky and ambient load
hex_array_sky = ['02', '01', '02', '90', '01', 'ff', 'ff', '01', '00', '00', '00', '00', '00', '00', '00', '00']
packet_sky = np.zeros(16, dtype='uint8')
for i, h in enumerate(hex_array_sky):
    packet_sky[i] = np.cast['uint8'](int(h, 16))
psky = packet_sky.tobytes()

hex_array_load = ['02', '01', '02', '90', '01', 'ff', 'ff', '01', '01', '00', '00', '00', '00', '00', '00', '00']
packet_load = np.zeros(16, dtype='uint8')
for i, h in enumerate(hex_array_load):
    packet_load[i] = np.cast['uint8'](int(h, 16))
pload = packet_load.tobytes()


def switch_cal_load():
    """
    Establishes socket connection to writing port for CAN modules on C1 and C2 to control calibration loads.
    Move the load into place, wait 15 seconds, and then move back to sky
    """
    # establish socket connection for C1
    s_listen_C1 = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s_listen_C1.connect((HOST1, PORT1))

    s_write_C1 = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s_write_C1.connect((HOST1, PORT2))

    # establish socket connection for C2
    s_listen_C2 = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s_listen_C2.connect((HOST2, PORT1))

    s_write_C2 = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s_write_C2.connect((HOST2, PORT2))

    try:

        s_write_C1.send(pload)
        s_write_C2.send(pload)

        time.sleep(15)

        s_write_C1.send(psky)
        s_write_C2.send(psky)

    except Exception as e:
        print(e)

    finally:
        s_listen_C1.close()
        s_write_C1.close()

        s_listen_C2.close()
        s_write_C2.close()
