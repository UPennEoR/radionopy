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
import healpy as hp
import numpy as np
from astropy import constants as c
from radiono.scripts import ion_altaz as ia

rad_path = os.path.dirname(os.path.realpath(__file__))
base_path = os.path.abspath(os.path.join(rad_path, '..'))
rm_dir = os.path.join(base_path, 'RM_files')
ionex_dir = os.path.join(base_path, 'TEC')
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

def get_results(hour, new_file, B_para, TEC_path, RMS_TEC_path, write_to_file=True):
    '''
    writes ionospheric RM to file

    Parameters
    ----------
    hour | str: standardized hour
    new_file | str: filename to write RM to
    B_para | array: B field parallel to zenith
    TEC_path | array: line of sight TEC
    RMS_TEC_path | array: line of sight RMS TEC
    write_to_file | bool: whether or not to write to file
    '''
    # Saving the Ionosheric RM and its corresponding
    # rms value to a file for the given 'hour' value
    IFR = 2.6e-17 * B_para * TEC_path
    RMS_IFR = 2.6e-17 * B_para * RMS_TEC_path

    if write_to_file:
        with open(new_file, 'w') as f:
            for tp, tf, ifr, rms_ifr in zip(TEC_path, B_para, IFR, RMS_IFR):
                f.write(('{hour} {TEC_path} '
                         '{B_para} {IFR} '
                         '{RMS_IFR}\n').format(hour=hour,
                                               TEC_path=tp,
                                               B_para=tf,
                                               IFR=ifr,
                                               RMS_IFR=rms_ifr))
    return B_para, IFR, RMS_IFR

class RM(object):
    '''
    superobject to do altaz or radec RM value generation
    gets the RM based on input

    Methods
    -------
    '''
    def __init__(self, lat_str, lon_str, time_strs, height=0, ionex_dir=ionex_dir, rm_dir=rm_dir):
        '''
        initalizes RM object

        Parameters
        ----------
        lat_str | str: latitude that the observation was taken at
        lon_str | str: longitude that the observation was taken at
        time_strs | list[str]: list of dates for observation
        height | Optional[float]: height observation taken at in meters
        ionex_dir | Optional[str]: directory in which ionex files are / will be located
        rm_dir | Optional[str]: directory in which RM data files are / will be located
        '''
        self.lat_str = lat_str
        self.lon_str = lon_str
        self.time_strs = time_strs
        self.height = height
        self.ionex_dir = ionex_dir
        self.rm_dir = rm_dir
        self.nside = 16

    def _radec(self, ra_strs, dec_strs):
        '''
        outputs RM data from radec calculation

        Parameters
        ----------
        ra_strs | list[str]: list of RAs for observation, corresponds by element to dec_strs
        dec_strs | list[str]: list of DECs for observation, corresponds by element for ra_strs
        '''
        #self.ra_strs = ra_strs
        #self.dec_strs = dec_strs
        #B_para, RMs, dRMs = 
        pass

    def _hp_arr(self):
        '''
        generates array of altitudes and azimuths for healpix map

        Returns
        -------
        tuple:
            array[float]: array of altitudes
            array[float]: array of azimuths
        '''
        npix = hp.nside2npix(self.nside)
        ipix = np.arange(npix)
        theta, phi = hp.pix2ang(self.nside, ipix)

        alt_src = 90. - np.degrees(np.array(theta))
        az_src = np.degrees(np.array(phi))

        return alt_src, az_src

    def _altaz(self, time_str):
        '''
        outputs RM data from altaz calculation

        Parameters
        ----------
        time_str | str: string representation of date and time

        Returns
        -------
        tuple:
            array[float]: parallel B field array
            array[float]: RM data
            array[float]: RM error data
        '''
        alt_src, az_src = self._hp_arr()
        B_para, RMs, dRMs = ia.ion_RM(time_str, self.lat_str, self.lon_str, alt_src, az_src, verbose=False)

        return B_para, RMs, dRMs

if __name__ == '__main__':
    print('This is not a script anymore')
