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
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, solar_system_ephemeris, get_body

OVRO = EarthLocation.of_site('ovro')

#################################################################################################################
# OBSERVATION PLANNING TOOLS
#################################################################################################################

def dqtau_nextper(nper=1):
    """
    Find next periastron date and time for binary system DQ Tau
    :param nper: number of future periastrons to return (default = 1)
    :return: list of dates in ISO format of next nper periastron
    :rtype: list
    """
    # get periastron reference time from Czekala et al (2016)
    # https://iopscience.iop.org/article/10.3847/0004-637X/818/2/156
    ref = 2400000
    T_per = Time(ref + 47433.507, format='jd')
    T_now = Time.now()
    P = 15.80158

    # number of full cycles that have already occurred
    n_cyc = int((T_now.jd - T_per.value) / P)
    next_per_list = []
    for n in np.arange(nper):
        # get date of next cycle
        next_per = T_per + ((n_cyc + n + 1) * P) * u.day
        next_per_list.append(next_per.iso)

    return next_per_list

def louis2014_Tb(nu):
    """
`   Get the brightness temperature of Uranus based on the model presented in Louis et al 2014
    :param nu: frequency in GHz
    :return: brightness temperature in Kelvin
    """
    a = np.log10(nu/100.0)
    Tb = 120 - 80.9*a + 22.2*(a**2)
    return Tb


# from Archinal et al. (2011)
# https://ui.adsabs.harvard.edu/abs/2011CeMDA.109..101A/abstract
r_eq = 25559e5 # +/- 4e5 cm
dr_eq = 4e5

r_pol = 24973e5 # +/- 20e5 cm
dr_pol = 20e5


Tb = louis2014_Tb(90) # K at 100 GHz -- get value for 90 GHz
dTb = 0.1*Tb # assume around 10% error
#unable to find a good source for this anywhere
# eyeballed off Fig 12 in https://www.aanda.org/articles/aa/full_html/2017/11/aa30311-16/aa30311-16.html

dD = 712000e5 # ephemeris errors as reported on the astropy website
# https://docs.astropy.org/en/stable/coordinates/solarsystem.html#precision-of-the-built-in-ephemeris

nu = 90e9 # GHz
k = 1.38e-16 # boltzmann
c = 3e10 # cm/s

Jy = 1e-23
def get_uranus_flux_density(obs_time):
    """
    Get the flux density of Uranus at the time it was observed
    :param obs_time: astropy time object of the start time of the observation
    :return: flux density and flux density error of Uranus in Jy
    :rtype: float, float
    """

    # get distance using astropy planet ephemeris
    with solar_system_ephemeris.set('builtin'):
        uranus = get_body('uranus', obs_time, OVRO)

    distance = uranus.distance.to(u.cm).value

    # get solid angle of Uranus given distance
    theta_eq = r_eq / distance
    theta_pol = r_pol / distance

    # get flux density and flux density error
    solid_angle = np.pi * theta_eq * theta_pol
    flux_density = 2 * (nu ** 2) * k * Tb * solid_angle / (c ** 2)

    flux_density_error = flux_density * np.sqrt(
        (dTb / Tb) ** 2 + (dr_eq / r_eq) ** 2 + (dr_pol / r_pol) ** 2 + (-2 * dD / distance) ** 2)

    return flux_density/Jy, flux_density_error/Jy

#################################################################################################################
# PARSING CATALOG AND SCHEDULE FILES NECESSARY FOR RUNNING CORRELATOR CODE
#################################################################################################################

# get all ALMA calibrators in sprite source catalog
def cat_to_dict(catalog_file):
    """
    Takes SPRITE catalog file and creates two dictionaries containing calibrator sources and target sources. \
    The source names are the keys and astropy SkyCoords objects for them are the data
    :param catalog_file: catalog file ending in ".cat"
    :return: calibrator dict, target dict
    :rtype: dict
    """

    with open(catalog_file, 'r') as f:
        lines = np.array(f.readlines())


    almacal_dict = {}
    ovrocal_dict = {}
    target_dict = {}

    idx1 = np.argwhere(lines == '# ALMA CALIBRATOR SOURCES\n').ravel()[0]
    idx2 = np.argwhere(lines == '# OVRO CALIBRATOR SOURCES\n').ravel()[0]
    idx3 = np.argwhere(lines == '# SCIENCE TARGETS\n').ravel()[0]

    for line in lines[idx1:idx2]:

        if line[0] == '#' or line[0] == '\n':
            continue

        else:

            if line[-1:] == '\n':
                line = line[:-1]

            if line[-1] == ' ':
                line = line[:-1]

            epoch, source, ra, dec = line.split(' ')

            coords = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))

            almacal_dict[source] = coords

    for line in lines[idx2:idx3]:

        if line[0] == '#' or line[0] == '\n':
            continue

        else:

            if line[-1:] == '\n':
                line = line[:-1]

            if line[-1] == ' ':
                line = line[:-1]

            epoch, source, ra, dec = line.split(' ')

            coords = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))

            ovrocal_dict[source] = coords

    for line in lines[idx3:]:

        if line[0] == '#' or line[0] == '\n':
            continue

        else:

            if line[-1:] == '\n':
                line = line[:-1]

            if line[-1] == ' ':
                line = line[:-1]

            epoch, source, ra, dec = line.split(' ')

            coords = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))

            target_dict[source] = coords

    cal_dict = {}
    cal_dict.update(almacal_dict)
    cal_dict.update(ovrocal_dict)

    return cal_dict, target_dict


def find_closest_cal(target_coords, cal_dict, print_cals=False):
    """
    Takes an input astropy SkyCoords object and looks through a given calibrator dictionary to identify \
    the nearest source to use as a calibrator. You can use the "cat_to_dict" function to create the dict \
    In the event that the source itself if in the calibrator dictionary, it will select the next nearest \
    source.
    :param target_coords: astropy SkyCoords object for target source
    :param cal_dict: dict of all calibrators
    :param print_cals: if True, the function will print the source name and separation distance
    :return: name of source located at smallest separation from target
    :rtype: str
    """

    seps = []

    for key in cal_dict:

        cal_coords = cal_dict[key]

        sep = target_coords.separation(cal_coords).value

        # allow observations with check sources to have a calibrator source
        # other than themselves
        if sep < 0.01:
            seps.append(np.inf)
        else:
            seps.append(sep)

    idx = np.argmin(seps)
    cal_source = list(cal_dict.keys())[idx]
    if print_cals:
        print('Best Calibrator: %s' % cal_source)
        print('Separation: %f degrees' % seps[idx])

    return cal_source

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

        if line[0] == '#' or line[0] == '\n':
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

def check_alt(sch):
    """
    Make sure for each source being observed, the source is
    :param sch: schedule list
    """

    for sch_row in sch:
        source_name, start_time, ra, dec, duration = sch_row
        obstime = start_time + np.linspace(0, duration, 100) * u.s
        src = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        altaz_frame = AltAz(obstime=obstime, location=OVRO)
        altaz = src.transform_to(altaz_frame)

        is_up = sum(1 * (altaz.alt.value < 30)) == 0

        if not is_up:
            raise Exception('Observation of %s at %s is not above 30 degrees altitude for entire duration'
                            % (source_name, start_time))


def check_timing(sch):
    """
    Make sure subsequent observations are at least 5 minutes apart from each other to allow sufficient \
    time for long slews
    :param sch: schedule list
    """

    for i in range(len(sch) - 1):
        obs0 = sch[i]
        obs1 = sch[i + 1]

        t0_end = obs0[1] + obs0[4] * u.s
        t0_end_plus5 = t0_end + 4.99 * u.min
        t1_start = obs1[1]

        if t0_end_plus5 > t1_start:
            raise Exception('Observation of %s at %s not at least 5 minutes after end of previous observation of'
                            % (obs1[0], str(obs1[1])) +
                            '%s ending at %s'% (obs0[0], t0_end))


def get_reference_date(sched_file):
    """
    Get the reference date for this schedule file from the corresponding line in the given schedule file. This would \
    be used for automatically advancing the times in a given schedule such that the source altitudes match those \
    on the reference day. To indicate a reference date in your schedule file, add the following line with the \
    approriate date filled in: "# reference date: YYYY-MM-DD." This function will raise an exception if no reference \
    date is given in the file
    :param sched_file: schedule file ending in .sch
    :return: the reference date
    :rtype: astropy Time
    """
    with open(sched_file) as f:
        lines = f.readlines()

    for line in lines:
        if 'reference date' in line:
            date = line.split(' ')[-1][:-1]

            date = date + " 00:00:00"

            return date

    # if you get to this line, then no ref date was given
    raise Exception('No reference date is given in this file.\n' +
                    'Please include the following line somewhere in this schedule file with the correct date:\n' +
                    '# reference date: YYYY-MM-DD')


def parse_schedule(sched_file, catalog_file):
    """
    Given a SPRITE schedule file and a source catalog file, create a list of observation details for each source to be \
    passed to the correlator code. Each row in the output array will contain the source name, start time, RA, Dec, and \
    observation length (in that order)
    :param sched_file: schedule file ending in ".sch"
    :param catalog_file: catalog file ending in ".cat"
    :param check_src_alt: Bool, default True - if True, raise exception for sources below 30 degrees
    :return: list of details from schedule file
    :rtype: numpy array
    """

    df = cat_to_df(catalog_file)

    with open(sched_file) as f:
        lines = f.readlines()

    obs_info = []
    bad_sources = []
    ref_time = Time.now()
    for line in lines:

        if line[0] == '#':
            continue

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
            obslen = 20 * n_step + 10  # in seconds

        elif 'five_point_scan' in line:
            _, source_name, start_time, offset = line.split(' ')

            obslen = 5 + (90*5) + 5  # in seconds

        else:
            continue

        source_name = source_name[:-1]
        start_time = start_time[:-1]

        start_time = Time('%s %s' % (get_date(), start_time))
        if start_time < ref_time:
            start_time = start_time + 1 * u.day


        if source_name in ['jupiter', 'saturn', 'uranus', 'moon', 'sun']:
            with solar_system_ephemeris.set('builtin'):
                p = get_body(source_name, start_time, OVRO)

            ra, dec = p.to_string('hmsdms').split(' ')
            ra = ra[:-1]
            dec = dec[:-1]
            for sym in ['h', 'd', 'm']:
                ra = ra.replace(sym, ':')
                dec = dec.replace(sym, ':')

        else:
            try:
                src_row = df[df['Source Name'] == source_name]
                ra = src_row['RA'].values[0]
                dec = src_row['Dec'].values[0]

            except:
                print('%s not in catalog'%source_name)

        source_row = [source_name, start_time, ra, dec, obslen]
        obs_info.append(source_row)

    return np.array(obs_info)


def sort_schedule(sched):
    """
    Take a schedule list and sort so that observation scheduled to occur soonest is scheduled first

    :param sched: list of schedule details created from "parse_schedule" function
    :return: list of details from schedule file sorted so most recent item appears first
    :rtype: numpy array
    """
    t_now = Time.now()

    row_start = 0
    dt_compare = 1 * u.day
    for i, sched_row in enumerate(sched):
        t_obs = sched_row[1]

        dt = t_obs - t_now

        if dt.value > 0 and dt.value < dt_compare.value:
            row_start = i
            dt_compare = dt

    new_sched = np.concatenate((sched[row_start:], sched[:row_start]))

    return new_sched


def advance_sched(sched, ref_date):
    """
    Given a schedule list and a reference date, adjust all schedule times such that the altitudes of the sources \
    will match those on the reference date. This is done by subtracting
    :param sched: list of schedule details created from "parse_schedule" function
    :param ref_date: reference date of schedule file
    :return: list of details from schedule file with adjusted time
    :rtype: numpy array
    """
    adv_sch = []

    diff = Time.now() - Time(ref_date, format='iso')
    n_days = int(diff.value)

    for sch_row in sched:
        source_name, start_time, ra, dec, duration = sch_row

        start_time_new = start_time - n_days * 3.93 * u.min

        while start_time_new < Time.now():
            start_time_new = start_time_new + 1 * u.day

        adv_sch.append([source_name, start_time_new, ra, dec, duration])

    return np.array(adv_sch)


def add_cal_scans(sch_row, cal_dict, target_dict):
    """
    For each source in schedule, find the near calibrator source and add routine scans of it

    :param sch_row: one row of schedule, created from "parse_schedule" or "sort_schedule" function
    :param cal_dict: dictionary containing all calibration sources and coordinates
    :return: list of details from schedule file, with all necessary calibrator scans added
    :rtype: numpy array
    """
    new_sched = []

    source_name, start_time, ra, dec, duration = sch_row
    # each observation should be up to 10 minutes on source, 1 minute on cal, plus 1 minute to slew each way
    # therefore each cycle is up to 13 minutes

    track = source_name
    if track == 'uranus':
        purpose = 'flux'
        new_sched.append([source_name, start_time, ra, dec, duration, track, purpose])

    elif track in ['3C84', '3C273']:
        purpose = 'bandpass'
        new_sched.append([source_name, start_time, ra, dec, duration, track, purpose])

    else:

        if duration <= 60 * 5:
            raise Exception('Science or check source %s must be observed for at least 5 minutes to leave enough time '+
                            'for calibrator scans'%source_name)

        try:
            cal_source = find_closest_cal(target_dict[source_name], cal_dict, print_cals=False)

        except:
            cal_source = find_closest_cal(cal_dict[source_name], cal_dict, print_cals=False)

        cal_ra, cal_dec = cal_dict[cal_source].to_string('hmsdms').split(' ')


        cal_ra = cal_ra[:-1]
        cal_dec = cal_dec[:-1]
        for sym in ['h', 'd', 'm']:
           cal_ra = cal_ra.replace(sym, ':')
           cal_dec = cal_dec.replace(sym, ':')

        if duration < 60 * 13:

            t = start_time

            duration_target = duration - (2 * 60)

            new_sched.append([source_name, t, ra, dec, duration_target, track, 'science'])
            t = t + duration_target * u.second + 1 * u.min
            new_sched.append([cal_source, t, cal_ra, cal_dec, 1 * 60, track, 'gain'])


        else:

            ncyc = int(duration / (60 * 13))

            t = start_time
            for i in range(ncyc):
                new_sched.append([source_name, t, ra, dec, 10 * 60, track, 'science'])
                t = t + 11 * u.min
                new_sched.append([cal_source, t, cal_ra, cal_dec, 1 * 60, track, 'gain'])
                t = t + 2 * u.min

            remain = duration - ncyc * (60 * 13)
            new_sched.append([source_name, t, cal_ra, cal_dec, remain, track, 'science'])

    return np.array(new_sched)

def get_full_schedule(sched_file, catalog_file, check_src_alt=True):
    """
    Given a SPRITE schedule file and a source catalog file, create a list of observation details for each source to be \
    passed to the correlator code. Each row in the output array will contain the source name, start time, RA, Dec, and \
    observation length (in that order). The rows will be sorted in order of most imminently occuring observation times \
    and will have all necessary calibration scans inserted into it. It will also write out a new .sch with the same name \
    and _full appended to it that contains all the scans
    :param sched_file: schedule file ending in ".sch"
    :param catalog_file: catalog file ending in ".cat"
    :param check_src_alt: Bool, default True - if True, raise exception for sources below 30 degrees
    :param uploadtoC2: Bool, default True - if True, scp the new schedule file into \
    comapc2:/home/comap/nyadlapa/schedules/.
    :return: list of details from schedule file
    :rtype: numpy array
    """

    sch = parse_schedule(sched_file, catalog_file)
    ref_date = get_reference_date(sched_file)
    adv_sched = advance_sched(sch, ref_date)
    new_sched = sort_schedule(adv_sched)

    if check_src_alt:
        check_alt(new_sched)

    check_timing(new_sched)

    cal_dict, target_dict = cat_to_dict(catalog_file)
    full_sched = np.array([]).reshape(0, 7)

    for sch_row in new_sched:
        sch_row_full = add_cal_scans(sch_row, cal_dict, target_dict)
        full_sched = np.concatenate((full_sched, sch_row_full), axis=0)

    split_path = sched_file.split(os.sep)
    filebase, ext = split_path[-1].split('.')
    split_path[-1] = filebase + '_full.' + ext
    full_sched_file = os.sep.join(split_path)

    lines = []
    lines.append('import /home/comap/nyadlapa/schedules/spriteSchedLib.sch\n')
    lines.append('catalog /home/comap/sprite_source_catalog.cat\n\n')

    command = 'observe_source'
    for row in full_sched:
        source, time, ra, dec, duration, track, purpose = row
        start_time = str(time).split(' ')[1][:8]
        duration = int(duration)
        lines.append('%s %s, %s, %ds\n' % (command, source, start_time, duration))

    lines.append('stow\n')

    with open(full_sched_file, 'w') as f:
        f.writelines(lines)

    cmd = 'scp %s comapc2:nyadlapa/schedules/.' %(full_sched_file)

    print('Upload full schedule to comap-centos machine with this command:\n'+
          '%s'%cmd+
          '\n')

    return full_sched





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
