'''
radiono

authors | James Aguirre, Immanuel Washington, Saul Kohn

purpose | Module used to gather information from IONEX files

Functions
---------
std_hour | converts hour into consistent string representation
get_results | writes ionospheric RM to files for future use
'''
from __future__ import print_function
import os
from astropy import constants as c

rad_path = os.path.dirname(os.path.realpath(__file__))
base_path = os.path.abspath(os.path.join(rad_path, '..'))
tec_dir = os.path.join(base_path, 'TEC')
TECU = 1e16
TEC2m2 = 0.1 * TECU
earth_radius = c.R_earth.value #6371000.0 # in meters
tesla_to_gauss = 1e4

def std_hour(UT, verbose=True):
    '''
    converts hour into standardized string

    Parameters
    ----------
    UT | int: numbered hour of time
    verbose | Optional[bool]: whether to print values or not

    Returns
    -------
    str: standardized string hour
    '''
    if verbose:
        print(int(UT))
    if UT < 10:
        hour = '0{hour}'.format(hour=int(UT))
    else:
        hour = '{hour}'.format(hour=int(UT))

    return hour

def get_results(hour, new_file, B_para, TEC_path, RMS_TEC_path):
    '''
    writes ionospheric RM to file

    Parameters
    ----------
    hour | str: standardized hour
    new_file | str: filename to write RM to
    B_para | array: B field parallel to zenith
    TEC_path | array: line of sight TEC
    RMS_TEC_path | array: line of sight RMS TEC
    '''
    # Saving the Ionosheric RM and its corresponding
    # rms value to a file for the given 'hour' value
    IFR = 2.6e-17 * B_para * TEC_path
    RMS_IFR = 2.6e-17 * B_para * RMS_TEC_path

    with open(new_file, 'w') as f:
        for tp, tf, ifr, rms_ifr in zip(TEC_path, B_para, IFR, RMS_IFR):
            f.write(('{hour} {TEC_path} '
                     '{B_para} {IFR} '
                     '{RMS_IFR}\n').format(hour=hour,
                                           TEC_path=tp,
                                           B_para=tf,
                                           IFR=ifr,
                                           RMS_IFR=rms_ifr))

if __name__ == '__main__':
    print('This is not a script anymore')
